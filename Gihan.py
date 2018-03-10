# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 21:59:02 2018

@author: gihan
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def Period(state,t):
    # unpack the initial state vector
    M = state[0] #concentration of mRNA
    Fc = state[1] #freuqency protein synthesised in cytoplasm
    Fn = state[2] #successively built up protein concentration
    # these are our constants
    #Vs = 2.0
    Vs = 1.5
    Vm = 1.5
    Vd = 1.0
    ks = 0.5
    K = 0.2
    Km = 0.15
    Kd = 0.15
    kin = 0.02
    kout = 0.1
    n = 4 #Hill number
    # compute state derivatives w.r.t time
    dM = (Vs*(K**n))/(K**n+Fn**n)-(Vm*M/(Km+M))
    dFc = ks*M - (Vd*(Fc/(Kd+Fc))) - (kin*Fc) + (kout*Fn)
    dFn = kin*Fc - kout*Fn
    # return the state derivatives
    return [dM, dFc, dFn]

def viewCircadianClock(result_odeint):
    fig, axes = plt.subplots(1,1, figsize=(5, 5))
    plt.plot(result_odeint[:,1], result_odeint[:,0])
    plt.xlabel("concentration of protein in cytoplasm")
    plt.ylabel("concentration of mRNA")
    plt.show()
    
def viewProductConcentrationEvolution(result_odeint, time):
    fig, axes = plt.subplots(1,1, figsize=(10, 10))
    labels = ["concentration of mRNA", "protein level in cytoskeleton", "protein level in nucleus"]
    plt.plot (time, result_odeint[:,0], 'b', time, result_odeint[:,1], 'g', time, result_odeint[:,2], 'r')
    #plt.plot(time, result_odeint[:,0], 'b')
    #plt.plot(time, result_odeint[:,1], 'g')
    #plt.plot(time, result_odeint[:,2], 'r')
    plt.legend(labels, loc='best') 
    plt.show()

"""
state0 = [0, 0, 0]
time = np.arange(0.0, 100, 1)
result_odeint = odeint(Period, state0, time)
"""

#viewCircadianClock(result_odeint)
#viewProductConcentrationEvolution(result_odeint, time)


def periodWithVarying_kin(state,t, kin):
    
    # Unpack initial state vector
    M = state[0] # mRNA concentration
    Fc = state[1] # frequency protein synthesised in cytoplasm
    Fn = state[2] # Successively built up protein concentration
    
    # Constants/Parmeters
    #Vs = 2.0
    Vs = 1.5
    Vm = 1.5
    Vd = 1.0
    ks = 0.5
    K = 0.2
    Km = 0.15
    Kd = 0.15
    kin = kin
    kout = 0.1
    n = 4 #Hill number
    
    # Compute state derivatives w.r.t time
    dM = (Vs*(K**n))/(K**n+Fn**n)-(Vm*M/(Km+M))
    dFc = ks*M - (Vd*(Fc/(Kd+Fc))) - (kin*Fc) + (kout*Fn)
    dFn = kin*Fc - kout*Fn
    
    # Return state derivatives
    return [dM, dFc, dFn]

# Calculate the mean period and return coordinates of peaks
# xlist is the corresponding time points of peaks
# maxindices is indices in the array representing peaks
def getMeanPeriodandPeakInfo(array,t):
    
    # Get the indices of the maxima
    # Index represents the no. of steps from start time
    maxindices = argrelextrema(array, np.greater)
    
    # List of all timepoints where M(t) is maximum
    xlist = []

    if len(maxindices[0]) > 1:
        for item in maxindices[0]:
            xlist.append(t[item])
            
    # Get the differences between values
    # Represents local period
    differences = np.diff(xlist)
    
    # Peak differences
    peakdiff = np.diff(maxindices[0])

    # Quarter values
    print(len(t))
    print("peakdiff/4",peakdiff/4)
    print("MAX",(len(t)-maxindices[0][-1])//(peakdiff/4))
    m = len(t)-maxindices[0][-1]
    q0value = np.array([maxindices[0][0]]*3 * np.array([0.25,0.5,0.75]))
    q1value = maxindices[0][:-1] + peakdiff/4
    q2value = maxindices[0][:-1] + peakdiff/4*2
    q3value = maxindices[0][:-1] + peakdiff/4*3
    q1last = array[(maxindices[0][-1]*0.25).astype(int)] if (m//(maxindices[0][-1]*0.25)!=0) else np.empty((0))
    q2last = array[(maxindices[0][-1]*0.5).astype(int)] if (m//(maxindices[0][-1]*0.5)!=0) else np.empty((0))
    q3last = array[(maxindices[0][-1]*0.75).astype(int)] if (m//(maxindices[0][-1]*0.75)!=0) else np.empty((0))
    print("Q1LAST",q1last)
    allqpositions = np.concatenate((q0value.astype(int),q1value.astype(int),q2value.astype(int),q3value.astype(int)))
    q0value = array[q0value.astype(int)]
    q1value = array[q1value.astype(int)]
    q2value = array[q2value.astype(int)]
    q3value = array[q3value.astype(int)]
    
    # Get the mean period
    meanP = np.average(differences)

    # Quarter time points
    q0 = np.array([xlist[0]]*3 * np.array([0.25,0.5,0.75])) # quarters before first peak
    q1 = xlist[:-1] + differences/4
    q2 = xlist[:-1] + differences*2/4
    q3 = xlist[:-1] + differences*3/4
    #qlast = 
    
    # Array values for quarters
    quartertimes = np.concatenate((q0,q1,q2,q3))
    quartervalues = np.concatenate((q0value,q1value,q2value,q3value))
    
    return meanP, xlist, maxindices, quartertimes, quartervalues, allqpositions
    

def findPeriod(minkin, maxkin,stepsize,state0, timeGap):
    
    # Start and end time in hours
    starttime = 0.0
    endtime = 60
    t = np.arange(starttime, endtime, timeGap)
    
    results = []
    
    # Evaluate on all kin values with certain stepsize
    kinvalues = np.arange(minkin,maxkin,stepsize)
    
    # Iterate through all the kin values
    for x in kinvalues:
        
        #print("kin: ",x)
        result_odeint = odeint(periodWithVarying_kin, state0, t, args=(x,))
        
        # Get results for M, Fc, and Fn
        currentkinresults = []
        for i in range(3):
            
            array = result_odeint[:,i]
            meanP,xlist,maxindices,quartertimes,quartervalues,allqpositions = getMeanPeriodandPeakInfo(array,t)

            #if 21.48 < meanP < 21.53:
            # Append a tuple of (kin, mean period, whole result, time array, max F values)
            currentkinresults.append((x, meanP, result_odeint, t,  xlist, array[maxindices], maxindices, quartertimes, quartervalues, allqpositions))

        results.append(currentkinresults)
    
    # Results is array of array
    # results[0] is results for first kin evaluated
    # results[0][0] is results for M for first kin
    # results[0][1] is results for Fc for first kin
    # results[0][2] is results for Fn for first kin
    #print(len(results))
    return results

def viewResultForPeriodsForProteinWithVarying_kin(results):
    
    kin = results[1][0]
    # Based on Fc
    period = results[1][1]
    result_odeint = results[1][2]
    time = results[1][3]
    peaksF = results[1][5]
    peaksM = results[0][2][:,0][results[1][6]]

    # Plot M(t) vs Fc(t)
    # With quarters (evening, midnight, morning, and noon)
    fig, axes = plt.subplots(1,1, figsize=(5, 5))
    plt.plot(result_odeint[:,1], result_odeint[:,0], c = 'blue')
    # Plot the quarter points
    plt.plot(results[1][8],results[0][2][:,0][results[1][9]], '.', c = 'orange')
    # Plot the eves based on Fc
    plt.plot(peaksF,peaksM, '.', c = 'black')
    plt.xlabel("concentration of protein in cytoplasm")
    plt.ylabel("concentration of mRNA")
    plt.show()
    
    # Plot all M, Fc, and Fn over time
    fig, axes = plt.subplots(1,1, figsize=(10, 10))
    labels = ["concentration of mRNA", "protein level in cytoskeleton", "protein level in nucleus"]
    plt.plot (time, result_odeint[:,0], 'blue', time, result_odeint[:,1], 'green', time, result_odeint[:,2], 'red')
    # Plot M quarter points - index is 0
    plt.plot(results[0][4],results[0][5], '.', c= 'black')
    plt.plot(results[0][7],results[0][8], '.', c= 'orange')
    # Plot Fc quarter points - index is 1
    plt.plot(results[1][4],results[1][5], '.', c= 'black')
    plt.plot(results[1][7],results[1][8], '.', c= 'orange')
    # Plot Fn quarter points - index is 2
    plt.plot(results[2][4],results[2][5], '.', c= 'black')
    plt.plot(results[2][7],results[2][8], '.', c= 'orange')
    plt.show()
    
minkin = 0.01
maxkin = 0.1
stepsize = .01
state0 = [0, 0, 0]
timeGap = 0.001


results = findPeriod(minkin, maxkin, stepsize, state0, timeGap)
targetkinresults= results[8]
#print(len(resultsM))
viewResultForPeriodsForProteinWithVarying_kin(targetkinresults)
#viewResultForPeriodsForProteinWithVarying_kin(results)

"""
def findPeriodBasedOnM(minkin, maxkin,stepsize,state0, timeGap):
    starttime = 0.0
    endtime = 1000
    t = np.arange(starttime, endtime, timeGap)
    Ps = []
    kinvalues = np.arange(minkin,maxkin,stepsize)
    # Iterate through all the kin values
    for x in kinvalues:
        result_odeint = odeint(periodWithVarying_kin, state0, t, args=(x,))
        # Array containing the values of M
        Marray = result_odeint[:,0]
        # Get the indices of the maxima
        # Index represents the no. of steps from start time
        maxindices = argrelextrema(Marray, np.greater)
         # List of all timepoints where M(t) is maximum
        xlist = []
        if len(maxindices[0]) > 1:
            for item in maxindices[0]:
                xlist.append(t[item])
        # Get the differences between values
        differences = np.diff(xlist)
        #print(differences)
        # Get the mean period
        meanP = np.average(differences)
        #print(meanP)
        if 21.48 < meanP < 21.53:
            # Append a tuple of (kin, mean period)
            Ps.append((x,meanP, result_odeint, t))
    return Ps           

def viewResultForPeriodsWithVarying_kin(results):
    for result in results:
        kin = result[0]
        period = result[1]
        result_odeint = result[2]
        time = result[3]
        print(kin, period)
        viewCircadianClock(result_odeint)
        viewProductConcentrationEvolution(result_odeint, time)
        




results = findPeriodBasedOnM(minkin, maxkin, stepsize, state0, timeGap)
print(len(results))
viewResultForPeriodsWithVarying_kin(results)
"""