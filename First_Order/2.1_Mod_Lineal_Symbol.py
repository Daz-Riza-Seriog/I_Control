# Code made for Sergio Andrés Díaz Ariza
# 18 March 2021
# License MIT
# Introduction to Control: Python Program Lineal and Jacobian

# Validez de un modelo linealizado
# En este código se comprueba la validez de un modelo linealizado en
# un punto de operación dado a partir de un modelo no lineal para un
# sistema de nivel

# Universidad Nacional de Colombia
# Lina Gómez 02/2010
# Modificado: Mario Giraldo 01/2015


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import seaborn as sns
sns.set()

# Params of the model

Atanque = 0.1963  #[m2]
g = 9.8           #[m/s]

# Esteady State of the System Estado estacionario del sistema

uss = 0.00049     #[m2]
dss = 0.001534 #[m/s]
yss = 1.1         #[m]

# Disturbance
duss = 0


## Comparisson between Lineal Model and Not-Lineal Model

# Variables Model
ynl = np.arange(0,2,0.1)      #variable del modelo no lineal
yl = np.arange(-1.1,0.9,0.1)  #variable de desviación   y1=deltay=y-yss

# Lineal Model
fnl=1/Atanque*(dss-np.sqrt(g)*uss*ynl**(1/2))

# Linealized Model(Deviation Variables)
fl=1/Atanque*(dss-1/2*np.sqrt(g)*uss*yss**(-1/2)*yl-np.sqrt(g)*yss**(1/2)*uss)

# Resoults Plots

plt.figure(1)
plt.plot(ynl,fnl,'b',label='Modelo no lineal')
plt.plot(ynl,fl,'r',label='Modelo linealizado')
plt.legend()

## Linealización de modelos usando comando jacobian

# Definición variables simbólicas
A = sp.Symbol('A')
G = sp.Symbol('G')
U = sp.Symbol('U')
D = sp.Symbol('D')
Y = sp.Symbol('Y')

# Creación de la función (modelo)
Fnl = 1/A*(D-sp.sqrt(G)*U*Y**(1/2))

# Linealización del modelo
F = sp.diff(Fnl,Y)

print(F)
plt.show()