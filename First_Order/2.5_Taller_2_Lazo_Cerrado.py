# Code made for Sergio Andrés Díaz Ariza
# 28 April 2021
# License MIT
# Introduction to Control: Python Program Assignment 1

from numpy.polynomial import Polynomial as P
import numpy as np
import matplotlib.pyplot as plt
import control as co
import seaborn as sns

sns.set()

# Create a Transfer

G_p = co.tf([5], [4, 1])
G_c = co.tf([2.16, 1.44], [1.5, 0])
L = co.series(G_c, G_p)

Y_s = co.feedback(L, 1)  # Create a New T.F with FeedBack
Y_s_2 = co.tf([5], [4, 8.2, 4.8])
roots = co.rlocus(Y_s_2)

# Plotting Roots, you don't need specify the plot
plt.figure(1)

# Check the roots with numpy
pol = P([4.8, 8.2, 4])
print(pol.roots())
print(Y_s)

plt.show()
