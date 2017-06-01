#Author Tristan VanFossen
#v1.3
#5-31-2017

##Changelog
##Changed to RHR function instead of HR. Performs reliability score and time decay weighting on the 10 PtP times stored to attempt to find more reliable data.
        ##RHR Function built from detailed writeup by Sean.
        ##Not implemented is the time decay on the reliability score before weighting, but directly reducing the weight of the HR by the time decay instead


import redis
import math
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time
import json
from collections import OrderedDict
from threading import Timer
from statsmodels import robust
from scipy.interpolate import interp1d




def redisConnect():
    return redis.StrictRedis(host='localhost',port=6379,db=0)

def findRelMax(buffer, threshold): ##Finds the first relative max
    for i in range (0,len(buffer)):
        if buffer[i]>threshold and i+20<len(buffer): ##Checks the relative max against the next 60 samples (63 samples per 250ms, there should be no peak higher than the selected value for atleast 60 samples
            for j in range (i, i+20):
                if buffer[i]<buffer[j]:
                    break
            else:
                return buffer[i],i
    else:
        return min(buffer), 63 #If no relative max is found, then return the minimum (For new low threshold) and the length of the buffer (To clear all data)


def clearBuffer(buffer, count): ##Clears the given buffer to the given point
    temp = np.array([])
    for i in range(0,count):
            temp = np.append(temp,i)
    return np.delete(buffer, temp)

def sigmoid(x):
    return 1/(1 + math.exp(-x))

def calcRHR(peaks):
    ibiBuffer = np.array([])
    ibiTimes = {}
    rIbiDict = {}
    avg = 0
    numReliable = 0

    for i in range (1, len(peaks)):
        if peaks[i] - peaks[i-1]>290 and peaks[i] - peaks[i-1]<1500: #Filter out any IBIs that are outside the threshold
            ibiBuffer = np.append(ibiBuffer, peaks[i] - peaks[i-1])
            ibiTimes[peaks[i]-peaks[i-1]] = (peaks[i] + peaks[i-1])/2 #Take the average of the times of the two peaks to determine a timestamp for later use
            numReliable+=1 #Number of IBIs (out of the 10) to be used for threshold decay change over time
    if (len(ibiBuffer)>0):
        ibiMedian = np.median(ibiBuffer) #Find the median of the IBIs
        ibiMad = robust.mad(ibiBuffer) #Find the median absolute deviation of the IBIs to form the actual cutoff

        uniformBuffer = np.array([]) #Creation of a normalized IBI set around the median of the actual set

        for i in range (1, 40):
            uniformBuffer = np.append(uniformBuffer, ibiMedian - i/38 *76* math.log10(76))
        for i in range (1,40):
            uniformBuffer = np.append(uniformBuffer, ibiMedian + i/38 *76* math.log10(76))

        uniformMad = robust.mad(uniformBuffer) #Find the median absolute deviation of the normalized set

        H = ibiMad/uniformMad #Ratio of the actual and normalized deviation to create a cutoff

        rDist = 1-min(0.99,H) # Reliability of the normalized distribution

        for i in ibiBuffer:
            rDelta = sigmoid(abs(i-ibiMedian)/(1.5*ibiMad)) #Score for each IBI as compared to the median
            reliability = math.sqrt(rDelta*rDist) #The score for each IBI is modified by the normalized score (rDist) and then square rooted to scale the drop off

            rIbiDict[i] = [reliability, ibiTimes[i]] #Creation of a dictionary for use in the reliability and time decay weighting

        curTime = (datetime.now().hour * 60 * 60 + datetime.now().minute * 60 + datetime.now().second) * 1000 + datetime.now().microsecond / 1000 #Time of "now"

        maxT = 0
        for i in rIbiDict: #Find the oldest peak in the dictionary, to be used as the scalar for the time decay
            if (curTime-rIbiDict[i][1])>maxT:
                maxT = curTime - rIbiDict[i][1]

        for i in rIbiDict: # generate the scaled sum of the IBIs
            avg+=(i*rIbiDict[i][0]*math.exp(-(curTime-rIbiDict[i][1])/(maxT)))

        w = 0
        for i in rIbiDict: #Find the scalar of the IBIs
            w+=(rIbiDict[i][0]*math.exp(-(curTime-rIbiDict[i][1])/(maxT)))

        if not math.isnan(60/(avg/w/1000)):
            return 60/(avg/w/1000), numReliable
        else: #If the weighted HR returns nan, check for weihted HR with ONLY reliability
            for i in rIbiDict:  # generate the scaled sum of the IBIs
                avg += (i * rIbiDict[i][0])

            w = 0
            for i in rIbiDict:  # Find the scalar of the IBIs
                w += (rIbiDict[i][0])

            if not math.isnan(60/(avg/w/1000)):
                return 60/(avg/w/1000), numReliable
            else: #If still giving a bad value, then calculate from old HR style
                return calcHR(peaks)




def calcHR(peaks): ##OLD HR CALCULATION
    avg = 0
    numReliable = 0

    for i in range (1, len(peaks)):
        if peaks[i] - peaks[i-1]>290 and peaks[i] - peaks[i-1]<1500: ##Min/max IBI
            #print(peaks[i], peaks[i-1])
            avg += peaks[i]-peaks[i-1]
            numReliable+=1
    if numReliable>len(peaks)/3:
        return 60/(avg/numReliable/1000), numReliable
    else:
        return -1, numReliable

def pushMsg(r, ZephyrID):

    for i in ZephyrID: ##Push both HRs to the DcsLogger feed every second
        t = time.mktime(datetime.now().timetuple())*1000 + int(datetime.now().microsecond/1000)

        msg = OrderedDict(
            [("timestamp", int(t)), ("device", "LPHR"), ("deviceId", str(i)), ("messageType", "refHR"),
             ("data", str(int(ZephyrID[i][6])))])

        json_msg = json.dumps(msg)

        r.lpush('DcsLogger', json_msg)  ##Push the data set to redis 'Dcslogger' with timestamp, ID, and HR\

    Timer(1, pushMsg, (r, ZephyrID)).start()



def ECG():

    print("LPHR starting")

    redisKey = 'DcsECG' ##Key to pull data from
    r.flushdb()

    numZephyr = 210

    maxVal = 0
    numReliable = 0

    Timer(1, pushMsg, (r, ZephyrID)).start()

    #plt.ion()
    #plt.plot()

    while True:
        temp = r.blpop(redisKey, timeout=5)

        if temp != None:
            temp = temp[1].decode('utf8')
            temp = json.loads(temp)

            for i in ZephyrID:
                if i == str(temp["deviceId"]) or str(temp["deviceId"]) == "--":
                    break
            else:
                numZephyr+=1
                ZephyrID[str(temp["deviceId"])] = [np.array([]),np.array([]),np.array([]),500, 0.95,pow(10, math.log10(0.95) * 0.1), -1,0, numZephyr]
                #1st - ECG buffer
                #2nd times buffer
                #2rd peaks buffer
                #4th threshold - initial threshold to 500
                #5th decay per sec - Initial to 0.9
                #6th thresholddecay - base from initial decay per sec
                #7th HR - init to -1
                #8th maxValTime

            if temp["device"] == 'BH3' and temp["messageType"] == 'EDP':

                currentZephyr = str(temp["deviceId"])  # Saves off the current Zephyr ID

                temp2 = temp['data']

                for i in range(0,len(temp2)-1): ##takes the current set of samples and appends them to the current buffer
                    if (int(temp2[i])!=None):
                        curTime = (datetime.now().hour*60*60 + datetime.now().minute*60+datetime.now().second)*1000+datetime.now().microsecond/1000

                        ZephyrID[currentZephyr][0] = np.append(ZephyrID[currentZephyr][0], int(temp2[i]))
                        ZephyrID[currentZephyr][1] = np.append(ZephyrID[currentZephyr][1], curTime)

                #plt.clf()
                #plt.plot(ZephyrID[currentZephyr][0])
                #plt.axhline(ZephyrID[currentZephyr][3])
                #plt.pause(0.01)


                if len(ZephyrID[currentZephyr][0])>10*63: ##When 10 sets of 63 samples are in the buffer, begin calculations

                    maxVal, ZephyrID[currentZephyr][7] = findRelMax(ZephyrID[currentZephyr][0], ZephyrID[currentZephyr][3])

                    ZephyrID[currentZephyr][2] = np.append(ZephyrID[currentZephyr][2], ZephyrID[currentZephyr][1][ZephyrID[currentZephyr][7]])

                    if len(ZephyrID[currentZephyr][2])>=11: ##Adjust the number of peaks to detect from 11 to 2 for rough PtP calculation - Further refinement of PtP for HRV
                        ZephyrID[currentZephyr][6], numReliable = calcRHR(ZephyrID[currentZephyr][2])

                        print(currentZephyr + ": " + str(ZephyrID[currentZephyr][6]))
                        
                        if ZephyrID[currentZephyr][6] != -1:
                            ZephyrID[currentZephyr][2] = clearBuffer(ZephyrID[currentZephyr][2], 1) ##Clear the last value from the peaks
                        elif len(ZephyrID[currentZephyr][2])>=15:
                            ZephyrID[currentZephyr][2] = clearBuffer(ZephyrID[currentZephyr][2],10)

                        if numReliable == 10 or numReliable == 9: # number of reliable peaks (IBI between 290 and 1500ms) is 9 or 10, further decrease the decay rate
                            if ZephyrID[currentZephyr][4]<0.95:
                                ZephyrID[currentZephyr][4] +=0.001
                            ZephyrID[currentZephyr][5] = pow(10, math.log10(ZephyrID[currentZephyr][4]) * 0.1)
                        elif numReliable <9  and numReliable>5: # Not enough reliable peaks, increase decay rate to catch more variation
                            if ZephyrID[currentZephyr][4] > 0.80:
                                ZephyrID[currentZephyr][4] -=0.001
                            ZephyrID[currentZephyr][5] = pow(10, math.log10(ZephyrID[currentZephyr][4]) * 0.1)

                    if ZephyrID[currentZephyr][7] != 0: ##If there is a clear peak found, clear buffer to that peak
                        ZephyrID[currentZephyr][0] = clearBuffer(ZephyrID[currentZephyr][0], ZephyrID[currentZephyr][7]+20)
                        ZephyrID[currentZephyr][1] = clearBuffer(ZephyrID[currentZephyr][1], ZephyrID[currentZephyr][7]+20)
                        ZephyrID[currentZephyr][7] = 0
                    else: ##Otherwise clear the last set of data
                        ZephyrID[currentZephyr][0] = clearBuffer(ZephyrID[currentZephyr][0], 63)
                        ZephyrID[currentZephyr][1] = clearBuffer(ZephyrID[currentZephyr][1], 63)

                    ZephyrID[currentZephyr][3] = maxVal ##Set the new threshold based on the current

                if ZephyrID[currentZephyr][3] != None: ##Adjust the threshold down by the decay rate
                    ZephyrID[currentZephyr][3] = ZephyrID[currentZephyr][3]*ZephyrID[currentZephyr][5] ##Find a way to adjust this variably? Base on numReliable found in HR calculation?

ZephyrID = {}
r = redisConnect()

ECG()
