# Code made for Sergio Andrés Díaz Ariza
# 18 March 2021
# License MIT
# Introduction to Control: Python Program Lineal and Jacobian

# Solucion de una EDO por
# Limage   octubre-2013
# Tanque de nivel dL/dt=Fo/Atanque-Cv*sqrt(g*L)/Atanque


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import seaborn as sns
sns.set()

#===========================================
#   VALORES DE ESTADO ESTACIONARIO
#===========================================
Lee=1          #   [m]
dee=0.001534   #   [m3/s]
uee=0.00049    #   [m2]

# =========================================================================
#     RESUELVE LA ECUACION DIFERENCIAL CON EL COMANDO ODE45
# =========================================================================
# Primero debo ingresar la ecuacion diferencial atraves de una funcion, se
# puede hacer en el mismo programa o en otro punto.m Aca se hace en el mismo
# programa

u = uee
d = dee
Tfinal = 3500

#===========================================
#   CONDICIONES INICIALES
#============================================
Lo = 0

#===========================================
#   DEFINIR LA EDO
#============================================

# Derivadas ordinarias de los estados T1 y T2
def Level_Tank(t,X,u,d):
    # =======se ingresan los parametros
    g = 9.8  # gravedad [m/s]
    Dtanque = 0.5  # [m]
    Atanque = (Dtanque / 2) ** 2 * np.pi  # [m2]

    # ======se ingresa las derivadas
    y = d / Atanque - u * np.sqrt(g * X) / Atanque
    return y

sol = solve_ivp(Level_Tank,[0,Tfinal],[Lo],args=(u,d),t_eval=np.arange(0,Tfinal,0.1))
print(sol.y[0])
# Ahora si se llama al comado ODE45
#[t X]=ode45(@EDO,[0 Tfinal],[Lo],[],u,d)
#L=X


plt.figure(1)  # efecto de la condición inicial
plt.plot(sol.t,sol.y[0])
plt.xlabel('Tiempo [s]')
plt.ylabel('Nivel [m]')
plt.title ('Evolución temporal del nivel en el tanque')
plt.show()
