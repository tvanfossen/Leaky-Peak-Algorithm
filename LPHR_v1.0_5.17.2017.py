#Author Tristan VanFossen
#v1.0
#5-17-2017

##Changelog
##Fixed error in clearing peaks from buffer
##Added rough self monitoring (Adjusting threshold decay variable up and down according to the number of reliable PtP values
##Removed maxValTime function and replaced with a second return from findRelMax function
##Added push to redis client and changed the push to unix time

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

def redisConnect():
    return redis.StrictRedis(host='localhost',port=6379,db=0)

def findRelMax(buffer, threshold): ##Finds the relative
    for i in range (0,len(buffer)):
        if buffer[i]>threshold and i+60<len(buffer): ##Checks the relative max against the next 60 samples (63 samples per 250ms, there should be no peak higher than the selected value for atleast 60 samples
            for j in range (i, i+60):
                if buffer[i]<buffer[j]:
                    break
            else:
                return buffer[i],i
    else:
        return min(buffer), len(buffer)-1


def clearBuffer(buffer, count): ##Clears the given buffer to the given point
    temp = np.array([])
    for i in range(0,count):
            temp = np.append(temp,i)
    return np.delete(buffer, temp)

def calcHR(peaks,times):
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

def ECG():

    r = redisConnect()
    #redisKey = 'DcsLogger' ##Key to pull data from
    redisKey = 'DcsECG' ##Key to pull data from
    r.flushdb()

    ZephyrID = None


    buffer = np.append([],[])
    times = np.append([],[]) ##Stored timestamps for each of the samples in buffer
    peaks = np.append([],[])

    maxVal = 0
    maxValTime = 0
    HR = 0
    numReliable = 0

    flag = False

    threshold = 500
    decayPerSec = 0.90 ##Adjust this value to for the decay rate ---- 0.85 works great for Tristan & Jim, 0.95 worked well for John... Work out self monitoring method
    thresholdDecay = pow(10, math.log10(decayPerSec)*0.1)

    #plt.ion()

    while True:
        #temp = r.blpop(redisKey)
        temp = r.lpop(redisKey)

        if temp != None:
            temp2 = temp.decode('utf8')
            temp = json.loads(temp2)

            print(temp2)

            if temp["device"] == 'BH3' and temp["messageType"] == 'EDP':
                ZephyrID = temp["deviceId"] # Saves off the current Zephyr ID

                temp2 = temp['data']

                for i in range(0,len(temp2)-1): ##takes the current set of samples and appends them to the current buffer
                    if (int(temp2[i])!=None):
                        buffer = np.append(buffer, int(temp2[i]))
                        times = np.append(times, (datetime.now().hour*60*60 + datetime.now().minute*60+datetime.now().second)*1000+datetime.now().microsecond/1000)  # Take the time in milliseconds for the current sample

                #plt.clf() ## Horizontal line on plot indicates the threshold adjusting over time, vertical line the position of the peak
                #plt.plot(buffer)
                #plt.axhline(threshold)
                #plt.axvline(maxValTime)
                #plt.pause(0.05)
                print(HR)

                if len(buffer)>10*63 or flag: ##When 15 sets of 63 samples are in the buffer, begin calculations
                    t = time.mktime(datetime.now().timetuple())*1000

                    msg = OrderedDict([("timestamp",int(t)),("device","LPHR"),("deviceId", str(ZephyrID)),("messageType", "refHR"),("data", str(int(HR)))])
                    json_msg = json.dumps(msg)

                    r.lpush('DcsLogger', json_msg)  ##Push the data set to redis 'Dcslogger' with timestamp, ID, and HR

                    maxVal, maxValTime = findRelMax(buffer, threshold)

                    if (len(buffer<3*63)): ##If the buffer ever drops below 7 sets of data stop the calculation
                        flag = False
                    else:
                        flag = True

                    if len(peaks) == 0:
                        peaks = np.append(peaks, times[maxValTime])
                    elif maxValTime != None:
                        peaks = np.append(peaks, times[maxValTime])


                    if len(peaks)>=11: ##Adjust the number of peaks to detect from 11 to 2 for rough PtP calculation - Further refinement of PtP for HRV
                        HR, numReliable = calcHR(peaks,times)
                        peaks = clearBuffer(peaks, 1) ##Clear the last value from the peaks

                        if numReliable == 10 or numReliable == 9 or numReliable<=5:
                            if decayPerSec<0.95:
                                decayPerSec +=0.001 ##Adjust the decay rate based on the numReliable
                            thresholdDecay = pow(10, math.log10(decayPerSec) * 0.1)
                        elif numReliable <9  and numReliable>5:
                            if decayPerSec > 0.80:
                                decayPerSec -=0.003
                            thresholdDecay = pow(10, math.log10(decayPerSec) * 0.1)

                    if maxValTime != 0: ##If there is a clear peak found, clear buffer to that peak
                        buffer = clearBuffer(buffer, maxValTime+20)
                        times = clearBuffer(times, maxValTime+20)
                    else: ##Otherwise clear the last set of data
                        buffer = clearBuffer(buffer, 63)
                        times = clearBuffer(times, 63)
                    threshold = maxVal ##Set the new threshold based on the current

                if threshold != None: ##Adjust the threshold down by the decay rate
                    threshold = threshold*thresholdDecay ##Find a way to adjust this variably? Base on numReliable found in HR calculation?



ECG()
