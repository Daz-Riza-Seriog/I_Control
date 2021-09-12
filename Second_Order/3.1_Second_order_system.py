# Code made for Sergio Andrés Díaz Ariza
# 08 Abril 2021
# License MIT
# Introduction to Control: Python Program Assignment 1

import numpy as np
import matplotlib.pyplot as plt
import control as co
import seaborn as sns
from sympy import Piecewise

sns.set()

# Create a Transfer Functions without Feedback

G1 = co.tf([1], [1, 0])
G2 = co.tf([1], [1, 0.01])
G3 = co.tf([1], [1, -0.01])

# Step Response

t_1 = np.linspace(0, 10, 1000)  # Times to evaluate
t_2 = np.linspace(0, 10, 1000)
t_3 = np.linspace(0, 10, 1000)

t1, y1 = co.step_response(G1, t_1)  # Steps Response for each T.F
t2, y2 = co.step_response(G2, t_2)
t3, y3 = co.step_response(G3, t_3)

plt.figure(1)  # Plotting Signal
plt.plot(t1, y1, label='G1')
plt.plot(t2, y2, label='G2')
plt.plot(t3, y3, label='G3')
plt.xlabel("time $[s]$")
plt.ylabel("Amplitude")
plt.legend()

###################################################################
# Create a Transfer Functions with Feedback
# Recall is G1/1+G1 --> for systems whitout T.F in Controller
##################################################################

G4 = co.feedback(G1, 1)  # Create a New T.F with FeedBack
G5 = co.feedback(G2, 1)
G6 = co.feedback(G3, 1)
# Step Response

t4, y4 = co.step_response(G4, t_1)  # Steps Response for each T.F
t5, y5 = co.step_response(G5, t_2)
t6, y6 = co.step_response(G6, t_3)

plt.figure(2)  # Plotting Signal
plt.plot(t4, y4, label='G4')
plt.plot(t5, y5, label='G5')
plt.plot(t6, y6, label='G6')
plt.xlabel("time $[s]$")
plt.ylabel("Amplitude")
plt.legend()

###################################################################
# Create a Transfer Functions with Feedback and Gain K (Controller)
###################################################################

G7 = co.tf([5], [1])  # Gain for the Controller in T.F

G8 = G1 * G7  # Gain times F.T initial
G9 = G2 * G7
G10 = G3 * G7

G_4 = co.feedback(G8, 1)  # Create a New T.F with FeedBack with Gain
G_5 = co.feedback(G9, 1)
G_6 = co.feedback(G10, 1)

t8, y8 = co.step_response(G_4, t_1)  # Steps Response for each T.F with K gain of controller
t9, y9 = co.step_response(G_5, t_2)
t10, y10 = co.step_response(G_6, t_3)

plt.figure(3)  # Plotting Signal
plt.plot(t8, y8, label='G4')
plt.plot(t9, y9, label='G5')
plt.plot(t10, y10, label='G6')
plt.xlabel("time $[s]$")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
