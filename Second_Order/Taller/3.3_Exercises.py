# Code made for Sergio Andrés Díaz Ariza
# 25 August March 2021
# License MIT
# Introduction to Control: Python Program Linealization interchange


import numpy as np
import matplotlib.pyplot as plt
import sympy
import sympy as sp
from sympy import *
import control as co
from scipy.integrate import solve_ivp
import seaborn as sns

sns.set()

#    1   POINT
# Define Sentiment Time

Kp = co.tf([3], [1])
Ki = co.tf([3], [1, 0])
K = co.tf([3], [3, 1])

K_ = Kp + Ki
Sys = co.series(K_, K)
G_fbck = co.feedback(Sys)

pole = co.pole(G_fbck)

plt.figure(1)
roots = co.rlocus(G_fbck)

print(Sys)
print(G_fbck)
print(pole)

#   3 POINT
G1 = co.tf([90], [1, 3, 90])
print("Pole G1", co.pole(G1))
G2 = co.tf([10], [1, 10, 30])
print("Pole G2", co.pole(G2))
G3 = co.tf([10], [1, 5, 34])
print("Pole G3", co.pole(G3))
G4 = co.tf([1], [1, 1, 1])
print("Pole G4", co.pole(G4))
plt.show()
