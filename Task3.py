# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 19:49:09 2018

@author: gihan
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


def Period12h12h(state,t,Vs):
  # unpack the initial state vector
  M = state[0]
  Fc = state[1]
  Fn = state[2]
  
  # Constants/Parameters
  
  Vm = 1.5
  Vd = 1.0
  ks = 0.5
  K = .2
  Km = 0.15
  Kd = 0.15
  kin = 0.02
  kout = 0.1
  n = 4 #Hill number
  
  # Transcription Rate
  Vs = Vs
  # Compute the derivatives
  dM = (Vs*(K**n))/(K**n + Fn**n) - (Vm*M/(Km+M))
  dFc = ks*M - (Vd*(Fc/(Kd+Fc))) - (kin*Fc) + (kout*Fc)
  dFn = kin*Fc - kout*Fn

  # Return the values
  return [dM, dFc, dFn]


#state0 = [0.6, 0.4, 0.4]
#tryp = Period12h12h(state0)

def Period(state,t):
  # unpack the initial state vector
  
  M = state[0]
  Fc = state[1]
  Fn = state[2]

  # these are our constants
  
  Vs = 2.0
  Vm = 1.5
  Vd = 1.0
  ks = 0.5
  K = .2
  Km = 0.15
  Kd = 0.15
  kin = 0.02
  kout = 0.1
  n = 4 #Hill number
  
  # compute state derivatives
  
  dM = (Vs*(K**n))/(K**n + Fn**n) - (Vm*M/(Km+M))
  dFc = ks*M - (Vd*(Fc/(Kd+Fc))) - (kin*Fc) + (kout*Fc)
  dFn = kin*Fc - kout*Fn

  # return the state derivatives
  return [dM, dFc, dFn]


# Function for finding the kin value that produces a period of P = 21.5
def findPeriod(state0):
    
    h = 0.01
    starttime = 0.0
    endtime = 100
    t = np.arange(starttime, endtime, h)
    #numsteps = (endtime - starttime)/h

    state0 = state0
    #stepsize = stepsize

    # Solve the equations
    state = odeint(Period, state0, t)

    # Plot the M(t) vs F(t)
    plt.figure(1)
    plt.plot(state[:,0],state[:,1])
    #plt.plot(t, y(t), label='exact' )
    plt.title("Goldbeter Model")
    plt.xlabel('F(t)') 
    plt.ylabel('M(t)')
    #plt.legend(loc=4)
    plt.grid()
    #plt.show()
    #plt.savefig('Mt_vs_Ft.png', fmt='PNG', dpi=100)
    """
    # Plot the M(t) vs t
    plt.figure(2)
    plt.plot(t,state[:,0])
    #plt.plot(t, y(t), label='exact' )
    plt.title("Goldbeter Model")
    plt.xlabel('t') 
    plt.ylabel('M(t)')
    #plt.legend(loc=4)
    plt.grid()
    #plt.savefig('Mt_vs_t.png', fmt='PNG', dpi=100)
    #plt.show()
    """
    # Array containing the values of M
    Farray = state[:,1]
    
    # Get the indices of the maxima
    # Index represents the no. of steps from start time
    maxindices = argrelextrema(Farray, np.greater)
    
    # Convert to hours
    xlist = []
    for item in maxindices[0]:
        xlist.append(item*h + starttime)
        
    # Get the differences between values
    differences = np.diff(xlist)
    
    # Get the mean period
    meanP = np.average(differences)
    
    return meanP

"""  
state0 = [0.6, 0.4, 0.4]
results = findPeriod(state0)
"""

# Function for finding the minimum Vs step increase 
# that produces a period of P = 21.5
def findVsStep(state0):
    
    # Time parameters
    h = 0.01 # Stepsize for time
    starttime = 0.0
    endtime = 100
    # Each t point is a step through time in the plot
    # the value of t is an index, not an absolute value in hours
    # Here, len(t) = (end-start)/stepsize = 100/0.01 = 10,000
    t = np.arange(starttime, endtime, h)
    
    # Split the timepoints into 12-hour arrays
    # Each value is a number in hours
    wheretosplit = [i for i in range(0,endtime,12)][1:]
    
    """
    wheretosplit
    Out[429]: [12, 24, 36, 48]
    """
    
    # Convert to specific indices in the t array
    indices = [int(x/h) for x in wheretosplit]
    print("Indices",indices)
    # Array of splitted timepoints based on 12 h cycle
    timearrays= np.split(t,indices)
    
    # Vs is the transcription rate
    step = 0.5
    vs0 = 1.5
    vs1 = vs0 + step
    Vs = vs0
    
    # Create the Vs vector (alternating vs values through light/dark)
    # For plotting later
    # numsteps = (endtime - starttime)/stepsize
    vsarray = []
    s = starttime
    stepsizeh = h
    for i in range(len(wheretosplit)):
        #Convert to num
        print("WHERE TO SPLIT",wheretosplit[i])
        numsteps = int((wheretosplit[i]-s)/stepsizeh)
        #print(numsteps)
        if (i%2) == 0:
            curvs = vs0
        else:
            curvs = vs1
        
        s = wheretosplit[i]
        currentvs = [curvs]*numsteps
        #print("Currentvsarray",len(currentvs))
        #print(currentvs)
        vsarray = vsarray + currentvs
        #currentvs = [i for i in range(int(numsteps))]
        #print("Currentvs",currentvs)
        # Current end time is next start time
    lastvalues = endtime - wheretosplit[-1]
    print("Lastvalues",lastvalues)
    numsteps = int((endtime - wheretosplit[-1])/stepsizeh)
    last = len(wheretosplit)+1
    if (last%2==0):
        lastarray = [vs0]*numsteps
    else:
        lastarray = [vs1]*numsteps
    vsarray = vsarray+lastarray

    # Numpy array where all results are appended
    # Initialize empty numpy array
    allstates = np.empty(shape=(0,3))
    #print(len(allstates))
    
    currstate = state0
    
    # Calculate stepwise 
    for x in enumerate(timearrays):
  
        # 12 h:12 h light:dark conditions
        if (x[0] % 2 == 0):
            # Even --> dark conditions
            Vs = vs0
        else:
            # Odd --> light conditions
            Vs = vs1
        
        # Solve the equations
        # x[1] is the current timearray
        state = odeint(Period12h12h, currstate, x[1], args=(Vs,))
        #print(state,len(state))
        print(allstates.shape)
        #print(state.shape)
        # Plot the M(t) vs F(t)
        plt.figure(1)
        plt.plot(state[:,0],state[:,1])
        #plt.plot(t, y(t), label='exact' )
        plt.title("Goldbeter Model")
        plt.xlabel('F(t)') 
        plt.ylabel('M(t)')
        #plt.legend(loc=4)
        plt.grid()
        #plt.show()
        #plt.savefig('Mt_vs_Ft.png', fmt='PNG', dpi=100)
        
        # Update the current state
        currstate = state[-1]
        
        allstates = np.concatenate((allstates,state))
        
        """
        # Plot the M(t) vs t
        plt.figure(2)
        plt.plot(t,state[:,0])
        #plt.plot(t, y(t), label='exact' )
        plt.title("Goldbeter Model")
        plt.xlabel('t') 
        plt.ylabel('M(t)')
        #plt.legend(loc=4)
        plt.grid()
        #plt.savefig('Mt_vs_t.png', fmt='PNG', dpi=100)
        #plt.show()
        """
    # Final added data
    print(allstates.shape)
    plt.figure(2)
    plt.plot(allstates[:,0],allstates[:,1])
    plt.figure(3)
    #plt.plot(allstates[:,0],allstates[:,1])
    plt.plot(t,vsarray)
    #plt.plot(allstates[:,0],allstates[:,1])

    print(vsarray)           
     
    """
    # Get initial state
    state0 = state0
    
    # Step size = increment of Vs during light conditions
    stepsize = 0.05
    maxinc = 2.5
    incvalues = np.arange(0.0,maxinc,stepsize)
    
    Ps = []

    # Try different values of the increment/step size
    for inc in incvalues:
        
        
        
        # Solve the equations
        state = odeint(Period12h12h, state0, t, args=(inc,))

        # Plot the M(t) vs F(t)
        plt.figure(1)
        plt.plot(state[:,0],state[:,1])
        #plt.plot(t, y(t), label='exact' )
        plt.title("Goldbeter Model")
        plt.xlabel('F(t)') 
        plt.ylabel('M(t)')
        #plt.legend(loc=4)
        plt.grid()
        #plt.show()
        #plt.savefig('Mt_vs_Ft.png', fmt='PNG', dpi=100)
        
        # Plot the M(t) vs t
        plt.figure(2)
        plt.plot(t,state[:,0])
        #plt.plot(t, y(t), label='exact' )
        plt.title("Goldbeter Model")
        plt.xlabel('t') 
        plt.ylabel('M(t)')
        #plt.legend(loc=4)
        plt.grid()
        #plt.savefig('Mt_vs_t.png', fmt='PNG', dpi=100)
        #plt.show()
        
        # Array containing the values of F
        Farray = state[:,1]
        
        # Get the indices of the maxima
        # Index represents the no. of steps from start time
        maxindices = argrelextrema(Farray, np.greater)
        
        # Convert to hours
        xlist = []
        for item in maxindices[0]:
            xlist.append(item*h + starttime)
        print(xlist)
            
        # Get the differences between values
        differences = np.diff(xlist)
        
        # Get the mean period
        meanP = np.average(differences)
        
        # Append a tuple of (kin, mean period)
        Ps.append((inc,meanP))
    
    return Ps
    """

state0 = [0.6, 0.4, 0.4]
results = findVsStep(state0)
