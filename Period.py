# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 22:04:51 2018

@author: gihan
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# https://www.gribblelab.org/compneuro/2_Modelling_Dynamical_Systems.html

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
  kin = 0.03
  kout = 0.1
  n = 4 #Hill number
  
  # compute state derivatives
  
  dM = (Vs*(K**n))/(K**n + Fn**n) - (Vm*M/(Km+M))
  dFc = ks*M - (Vd*(Fc/(Kd+Fc))) - (kin*Fc) + (kout*Fc)
  dFn = kin*Fc - kout*Fn

  # return the state derivatives
  return [dM, dFc, dFn]

state0 = [0.6, 0.4, 0.4]
t = np.arange(0.0, 100, 0.01)

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
plt.savefig( 'Mt_vs_Ft.png', fmt='PNG', dpi=100 )

# Plot the M(t) vs t
plt.figure(2)
plt.plot(t,state[:,0])
#plt.plot(t, y(t), label='exact' )
plt.title("Goldbeter Model")
plt.xlabel('t') 
plt.ylabel('M(t)')
#plt.legend(loc=4)
plt.grid()
plt.savefig( 'Mt_vs_t.png', fmt='PNG', dpi=100 )

# Task 2b

# Find the local maxima for M(t) vs t
from scipy.signal import argrelextrema

# Array containing the values of M
Marray = state[:,1]

# Get the indices of the maxima#
# Index represents the no. of steps from start time
maxindices = argrelextrema(Marray, np.greater)

# Convert to hours
xlist = []
for item in maxindices[0]:
    xlist.append(item*stepsize + starttime)
print("xlist: ",xlist)
print(xlist[1]-xlist[0])

"""
https://stackoverflow.com/questions/25891972/solve-a-ode-for-arrays-of-parameters-python
# do some fancy 3D plotting
from mpl_toolkits.mplot3d import Axes3D
fig = figure()
ax = fig.gca(projection='3d')
ax.plot(state[:,0],state[:,1],state[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
show()

#FINDING LOCAL EXTREMA
import numpy as np
from scipy.signal import argrelextrema

x = np.random.random(12)

# for local maxima
argrelextrema(x, np.greater)

# for local minima
argrelextrema(x, np.less)

plt.figure(2)
# Plot the M(t) vs t
plt.plot(t,state[:,0],label='approximation')
#plt.plot(t, y(t), label='exact' )
#plt.title( "Euler's Method Example, N="+str(N) )
plt.xlabel('t') 
plt.ylabel('M(t)')
plt.legend(loc=4)
plt.grid()
plt.savefig( 'Mt_vs_t.png', fmt='PNG', dpi=100 )
"""