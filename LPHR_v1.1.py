#Author Tristan VanFossen
#v1.0
#5-17-2017

##Changelog
##Added ability to store data from multiple Zephyrs

#Needed:
###Additional buffer to output instant PtP/IBI for HRV
###Further refinements in peak detection, seems to have error in picking the correct position
###Statistical adjustment of the data - Hold an average heart rate or average HRV to adjust the decay rate based on typical changes
    #Sean's writeup makes mention of adjusting the decay based on the age of user...?
###More informative commenting....

import redis
from scipy.signal import butter, lfilter
import math
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time
import json
from collections import OrderedDict
from threading import Timer



def redisConnect():
    return redis.StrictRedis(host='localhost',port=6379,db=0)

def findRelMax(buffer, threshold): ##Finds the first relative max
    for i in range (0,len(buffer)):
        if buffer[i]>threshold and i+40<len(buffer): ##Checks the relative max against the next 60 samples (63 samples per 250ms, there should be no peak higher than the selected value for atleast 60 samples
            for j in range (i, i+40):
                if buffer[i]<buffer[j]:
                    break
            else:
                return buffer[i],i
    else:
        return min(buffer), len(buffer)-1 #If no relative max is found, then return the minimum (For new low threshold) and the length of the buffer (To clear all data)


def clearBuffer(buffer, count): ##Clears the given buffer to the given point
    temp = np.array([])
    for i in range(0,count):
            temp = np.append(temp,i)
    return np.delete(buffer, temp)

def calcHR(peaks):
    avg = 0
    numReliable = 0
    for i in range (1, len(peaks)):
        if peaks[i] - peaks[i-1]>290 and peaks[i] - peaks[i-1]<1500: ##Min/max IBI
            #print(peaks[i], peaks[i-1])
            avg += peaks[i]-peaks[i-1]
            numReliable+=1
    if numReliable>len(peaks)/2: ##If there is not enough reliable data, return -1 otherwise return HR and the number of reliable IBIs
        return 60/(avg/numReliable/1000), numReliable
    else:
        return -1, numReliable


def pushMsg(r, ZephyrID):
    for i in ZephyrID: ##Push both HRs to the DcsLogger feed every second
        t = time.mktime(datetime.now().timetuple())

        msg = OrderedDict(
            [("timestamp", int(t)), ("device", "LPHR"), ("deviceId", str(i)), ("messageType", "refHR"),
             ("data", str(int(ZephyrID[i][6])))])

        print("check")
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

    while True:
        temp = r.blpop(redisKey, timeout=5)

        if temp != None:
            temp = temp[1].decode('utf8')
            temp = json.loads(temp)

            for i in ZephyrID:
                if i == str(temp["deviceId"]):
                    break
            else:
                numZephyr+=1
                ZephyrID[str(temp["deviceId"])] = [np.array([]),np.array([]),np.array([]),500, 0.9,pow(10, math.log10(0.9) * 0.1), -1,0, numZephyr]
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

                if len(ZephyrID[currentZephyr][0])>10*63: ##When 10 sets of 63 samples are in the buffer, begin calculations

                    maxVal, ZephyrID[currentZephyr][7] = findRelMax(ZephyrID[currentZephyr][0], ZephyrID[currentZephyr][3])

                    ZephyrID[currentZephyr][2] = np.append(ZephyrID[currentZephyr][2], ZephyrID[currentZephyr][1][ZephyrID[currentZephyr][7]])

                    if len(ZephyrID[currentZephyr][2])>=11: ##Adjust the number of peaks to detect from 11 to 2 for rough PtP calculation - Further refinement of PtP for HRV
                        ZephyrID[currentZephyr][6], numReliable = calcHR(ZephyrID[currentZephyr][2])

                        print(currentZephyr + ": " + str(ZephyrID[currentZephyr][6]))
                        
                        ZephyrID[currentZephyr][2] = clearBuffer(ZephyrID[currentZephyr][2], 1) ##Clear the last value from the peaks

                        if numReliable == 10 or numReliable == 9: # number of reliable peaks (IBI between 290 and 1500ms) is 9 or 10, further decrease the decay rate
                            if ZephyrID[currentZephyr][4]<0.95:
                                ZephyrID[currentZephyr][4] +=0.001
                            ZephyrID[currentZephyr][5] = pow(10, math.log10(ZephyrID[currentZephyr][4]) * 0.1)
                        elif numReliable <9  and numReliable>5: # Not enough reliable peaks, increase decay rate to catch more variation
                            if ZephyrID[currentZephyr][4] > 0.80:
                                ZephyrID[currentZephyr][4] -=0.003
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
