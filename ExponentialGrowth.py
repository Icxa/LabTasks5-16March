# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:03:53 2018

@author: gihan
"""

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

"""
r = 0
K = 0

# Righthand side of differential equation
def f(N):
    return r*N*(1-(N/K))

# Define initial condition
x0 = 0

# Time step
dt = 0.01

# Solve the differential equation from time 0 to time T
T = 5

# Define discretized time ; assumes dt divides nicely into T
t = np.arange(0, T,int(T/dt)+1)

# An array to store the solution
x = np.zeros(len(t))

# Integrate the differential equation using Euler ’s method
x[0] = x0
for i in xrange (1 , len(t)):
x[i] = x[i -1] + f(x[i -1])* dt
"""

# limits: 0.0 <= t <= 150
a = 0
b = 150
 
# Number of steps
nsteps = 1000

# Carrying capacity
K = 189.4

# Growth rate
r = 0.029443
 
# step-size
h = (b-a)/nsteps
 
# initial value: y(0.0) = 0.5
initial = (0.0,5.3)
 
# Logistic growth rate differential equation
def f(N,K,r):
    return r*N*(1-(N/K))
 
# Create arrays to hold values of t and N
t = np.arange(a,b+h,h)
narray = np.zeros((nsteps+1,))
 
# set the initial values
t[0], narray[0] = initial
 
# Apply Euler's method
for i in range(1,nsteps+1):
    narray[i] = narray[i-1] + h*f(narray[i-1], K, r)
    
# Exact solution
def y(t):
    return K / (1 + (((K-5.3)/5.3)*np.exp(-r*t)))

#plt.plot(t, narray, label='approximation' )
plt.plot(t, y(t), label='exact' )
#plt.title( "Euler's Method Example, N="+str(N) )
plt.xlabel('t') 
plt.ylabel('y(t)')
plt.legend(loc=4)
plt.grid()
plt.savefig( 'euler_example.png', fmt='PNG', dpi=100 )

"""

# time intervals
tt = np.arange(0, 10, 0.5)

# initial condition
xx = [0.1]

def (x):
    return x * (1.-x)

# loop over time
for t in tt[1:]:
    xx.append(xx[-1] + 0.5 * f(xx[-1]))

# plotting
plot(tt, xx, '.-')
ta = arange(0, 10, 0.01)
plot(ta, 0.1 * exp(ta)/(1+0.1*(exp(ta)-1.)))
xlabel('t')
ylabel('x')
legend(['approximation', 'analytical solution'], loc='best',)

"""
"""
r = .25 # growth rate / yr
K = 100 # carrying capacity
t = 40 # number of years
num = np.zeros(t+1)
num[0] = 1
for i in range(t):
    num[i+1] = num[i]+r*num[i]*(1-num[i]/K)
plt.plot(range(t+1),num, 'b')
plt.xlabel('Year')
plt.ylabel('Number')
plt.title('Growth rate: 0.25, Carrying Capacity = 100')
plt.axvline(np.argmax(np.diff(num)),  color = 'k'  )
plt.show()
"""