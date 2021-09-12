# Code made for Sergio Andrés Díaz Ariza
# 05 Abril 2021
# License MIT
# Introduction to Control: Python Program Assignment 1

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control as co
import sympy as sp

import seaborn as sns
sns.set()


# Define Transafer Function

G1 = co.tf([1],[1,2])
G2 = co.tf([1],[1,10])

# Convert Transfer function to Steady State
G1_ss = co.tf2ss(G1)
G2_ss = co.tf2ss(G2)

# Looking time Responses

# Impulse Response
t = np.linspace(0,4,1000)
t1,y1 = co.impulse_response(G1,t)

plt.figure(1)
plt.plot(t1,y1)
plt.xlabel("time $[s]$")
plt.ylabel("Amplitude")


# Step Response

t = np.linspace(0,4,1000)
t2 = np.linspace(0,4,1000)
t1,y1 = co.step_response(G1,t)
t2,y2 = co.step_response(G2,t2)

plt.figure(2)
plt.plot(t1,y1,label='G1')
plt.plot(t2,y2,label='G2')
plt.xlabel("time $[s]$")
plt.ylabel("Amplitude")
plt.legend()


####################################################
# Using Scipy and compare

lti = signal.lti([1],[1,2])
lti2 = signal.lti([1],[1,10])

t3, y3 = signal.step2(lti)
t4, y4 = signal.step2(lti2)

plt.figure(3)
plt.plot(t3, y3)
plt.plot(t4, y4)
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Step response for 1. Order ')

##########################################################
# Step1.2 Poitn 1 of Assignment
lti5 = signal.lti([12],[6,1])
lti6 = signal.lti([12],[12,1])
t5, y5 = signal.step2(lti5)
t6, y6 = signal.step2(lti6)

plt.figure(4)
plt.plot(t5, y5,label=r'$\tau=5$')
plt.plot(t6, y6,label=r'$\tau=10$')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Step response for 1. Order ')
plt.legend()


##########################################################
# Step1.2 Point 1 of Assignment
lti7 = signal.lti([4],[7,1])
t7, y7 = signal.step2(lti7)

plt.figure(5)
plt.plot(t7, y7,label=r'$\tau=10$')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Step response for 1. Order ')

# Using sympy library to find the inverse laplace
s = sp.symbols('s')
t = sp.symbols('t')

F = 8/(6.4*s+1)

F_1 = sp.inverse_laplace_transform(F,s,t)

print(F_1)

plt.show()
