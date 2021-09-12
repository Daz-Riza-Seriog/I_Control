# Code made for Sergio Andrés Díaz Ariza
# 24 August 2021
# License MIT
# Introduction to Control: Python Program Linealization interchange

import numpy as np
import matplotlib.pyplot as plt
import sympy
from scipy import signal
import sympy as sp
from sympy import *
import control as co
from scipy.integrate import solve_ivp
import seaborn as sns

sns.set()
# ======================================================================
#     CSTR-Exothermic from Chemistry Engineering
# ======================================================================

#             ==================
#  C_air,Tir,F=     CSTR       =   C_ao,To
#             ==================
#                   =
#                 = = =
#                =  =  =
#                   =     Q
#                   =
#             ==================
#   Fcr,Tci,  =     tubos       =  Tcs
#             ==================

#  Mass Balance reactant A
#   Vi*dX1/dt=F(C_air-X1)-V*K_o*np.exp((-E/R)*(1/X2))*X1**2

#  Energy Balance of Reactant A
#  roo_a*Vi*Cp_a*dX2/dt=mo_a*Cp_a*(Ti-X2)-Delta_Hr*V*(-r_a)-Q
#  r_a = K_o*np.exp((-E/R)*(1/X2))*X1**2
#  Q = Uor*Ao*(X2 - ((Ao*Uor*X2-Fcr*roo_c*Cp_c*Tci)/(Fcr*roo_c*Cp_c+Ao*Uor)))

# ======================================================================
#    CSTR-Exothermic from Control Engineering
# ======================================================================

#  Declaración de variables
#   U = Fcr
#   d1 = C_air
#   d2 = Ti
#   d3 = Tci
#   y = X2 = T

#   Diagrama de bloques

#                         |
#                         | d=C_air,Ti,Tci
#                         V
#                 =================
#      u=Fcr ---> =               = ---> y = X2 = T
#                 =================
#
#

#   Balance in Control language
#   Vi*dX1/dt=F(d1-X1)-V*K_o*np.exp((-E/R)*(1/X2))*X1**2

#   roo_a*Vi*Cp_a*dX2/dt=mo_a*Cp_a*(d2-X2)-Delta_Hr*V*(-r_a)-Q
#   r_a = K_o*np.exp((-E/R)*(1/X2))*X1**2
#   Q = Uor*Ao*(X2 - ((Ao*Uor*X2-U*roo_c*Cp_c*d3)/(U*roo_c*Cp_c+Ao*Uor)))

#   Physical Properties
rho = 850  # [Kg/m^3] Density of the fluid
rho_c = 997  # [Kg/m^3] Density of the coolant
Cp = 3800  # [J/Kg*K] Heat Capacity of the fluid
Cp_c = 4184  # [J/Kg*K] Heat Capacity of the fluid
Delta_Hr = -450000  # [J/mol] Reaction Enthalpy
E_R = 8200  # [K] Activation Energy/Universal Constant Gas
Ko = 1100  # [m^3/mol*s] Rate Constant

#   Design and Operating Conditions
Uor = 438  # [W/m^2*K] Overall Heat Transfer Coefficient
Ao = 18  # [m^2] Outer Surface Area
V = 0.785  # [m^3] Volume of Reactant
F = 0.015  # [m^3/s] Caudal or Flux of Reactant
Tci = 300  # [K] Coolant Temperature inlet

# Reference Conditions
Tir = 330  # [K] Temperature inlet reference
C_air = 2300  # [mol/m^3] Concentration reference
Fcr = 0.008  # [m^3/s] Caudal or flux of Reference of the coolant (Water)


# ==========================================================================
#        NO-LINEAR SYSTEM
# ==========================================================================

# Solve the Non-Linear System
def Reactor_Exo(t, state):
    #   Physical Properties
    rho = 850  # [Kg/m^3] Density of the fluid
    rho_c = 997  # [Kg/m^3] Density of the coolant
    Cp = 3800  # [J/Kg*K] Heat Capacity of the fluid
    Cp_c = 4184  # [J/Kg*K] Heat Capacity of the fluid cooler
    Delta_Hr = -450000  # [J/mol] Reaction Enthalpy
    E_R = 8200  # [K] Activation Energy/Universal Constant Gas
    Ko = 1100  # [m^3/mol*s] Rate Constant

    #   Design and Operating Conditions
    Uor = 438  # [W/m^2*K] Overall Heat Transfer Coefficient
    Ao = 18  # [m^2] Outer Surface Area
    V = 0.785  # [m^3] Volume of Reactant
    F = 0.015  # [m^3/s] Caudal or Flux of Reactant
    Tci = 300  # [K] Coolant Temperature inlet

    # Reference Conditions
    Tir = 330  # [K] Temperature inlet reference
    C_air = 2300  # [mol/m^3] Concentration reference
    Fcr = 0.008  # [m^3/s] Caudal or flux of Reference of the coolant (Water)

    # Assign each ODE to a Vector
    C, T = state

    # Define each ODE
    dCdt = (F / V) * (C_air - C) - Ko * sp.exp(-E_R * (1 / T)) * C ** 2
    dTdt = (F / V) * (Tir - T) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / T) * C ** 2) - \
           (Uor * Ao / rho * V * Cp) * (
                   T - ((Ao * Uor * T - Fcr * rho_c * Cp_c * Tci) / (Fcr * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

    return [dCdt, dTdt]


# initial conditions
X0 = [2296, 326.7]
t_span = (0.0,1000.0)
t = np.linspace(0, 100, 100000)

# Solve the ODEs
sol_non_lin_sys = solve_ivp(Reactor_Exo,t_span , X0, method='LSODA', t_eval=t)

C_sol = sol_non_lin_sys.y[0, :]
T_sol = sol_non_lin_sys.y[1, :]
print(T_sol)

plt.figure(4)
plt.axhline(y=sol_non_lin_sys.y[0, :], color='b', linestyle='-')
plt.axhline(y=sol_non_lin_sys.y[1, :], color='r', linestyle='-')
plt.title('INTERCHANGER\nTemperature $[^{\circ}C]$ & Concentration $[s]$', fontsize=16)
plt.xlabel("Time \t$[s]$ ", fontsize=14)
plt.ylabel("Temperature $[^{\circ}C]$\t&\tConcentration $C_A$ ", fontsize=14)
plt.legend()
plt.show()

# ==========================================================================
#         LINEALIZATION
# ==========================================================================
# Definición variables simbólicas

# se declaran simbolicas
U = Fcr  # uEE
d1 = C_air  # d1EE
d2 = Tir  # d2EE
d3 = Tci  # d3EE

X1, X2 = sp.symbols('X1,X2', real=True)
f1 = (F / V) * (d1 - X1) - Ko * sp.exp(-E_R * (1 / X2)) * X1 ** 2  # f1=dx1/dt
f2 = (F / V) * (d2 - X2) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / X2) * X1 ** 2) - \
     (Uor * Ao / rho * V * Cp) * (
             X2 - ((Ao * Uor * X2 - U * rho_c * Cp_c * d3) / (U * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

A_a = sympy.Matrix([f1, f2]).jacobian([X1, X2])
s = (X1, X2)
A_a = sp.lambdify(s, A_a, "numpy")
X1 = 2296  # x1EE
X2 = 326.7  # x2EE
Aa = A_a(X1, X2)

U = sp.symbols('U', real=True)
f1 = (F / V) * (d1 - X1) - Ko * sp.exp(-E_R * (1 / X2)) * X1 ** 2  # f1=dx1/dt
f2 = (F / V) * (d2 - X2) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / X2) * X1 ** 2) - \
     (Uor * Ao / rho * V * Cp) * (
             X2 - ((Ao * Uor * X2 - U * rho_c * Cp_c * d3) / (U * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

B = sympy.Matrix([f1, f2]).jacobian([U])
s1 = U
B_b = sp.lambdify(s1, B, "numpy")
U = Fcr  # uEE
Bb = B_b(U)

d1 = sp.symbols('d1', real=True)
f1 = (F / V) * (d1 - X1) - Ko * sp.exp(-E_R * (1 / X2)) * X1 ** 2  # f1=dx1/dt
f2 = (F / V) * (d2 - X2) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / X2) * X1 ** 2) - \
     (Uor * Ao / rho * V * Cp) * (
             X2 - ((Ao * Uor * X2 - U * rho_c * Cp_c * d3) / (U * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

F1 = sympy.Matrix([f1, f2]).jacobian([d1])
s2 = d1
F_f1 = sp.lambdify(s2, F1, "numpy")
d1 = C_air  # dEE
Ff1 = F_f1(d1)

d2 = sp.symbols('d2', real=True)
f1 = (F / V) * (d1 - X1) - Ko * sp.exp(-E_R * (1 / X2)) * X1 ** 2  # f1=dx1/dt
f2 = (F / V) * (d2 - X2) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / X2) * X1 ** 2) - \
     (Uor * Ao / rho * V * Cp) * (
             X2 - ((Ao * Uor * X2 - U * rho_c * Cp_c * d3) / (U * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

F2 = sympy.Matrix([f1, f2]).jacobian([d2])
s3 = d2
F_f2 = sp.lambdify(s3, F2, "numpy")
d2 = Tir  # dEE
Ff2 = F_f2(d2)

d3 = sp.symbols('d3', real=True)
f1 = (F / V) * (d1 - X1) - Ko * sp.exp(-E_R * (1 / X2)) * X1 ** 2  # f1=dx1/dt
f2 = (F / V) * (d2 - X2) - (Delta_Hr / rho * Cp) * (Ko * sp.exp(-E_R / X2) * X1 ** 2) - \
     (Uor * Ao / rho * V * Cp) * (
             X2 - ((Ao * Uor * X2 - U * rho_c * Cp_c * d3) / (U * rho_c * Cp_c + Ao * Uor)))  # f2=dx2/dt

F3 = sympy.Matrix([f1, f2]).jacobian([d3])
s4 = d3
F_f3 = sp.lambdify(s4, F3, "numpy")
d3 = Tci  # dEE
Ff3 = F_f3(d3)

C = [0, 1]
D = 0

#   Función de transferencia Y(s)/U(s)
Ys_Us_ss = co.ss(Aa, Bb, C, D)
Ys_Us = co.ss2tf(Ys_Us_ss)
print(Ys_Us)

#   Función de transferencia Y(s)/D1(s) --> C_air
Ys_D1s_ss = co.ss(Aa, Ff1, C, D)
Ys_D1s = co.ss2tf(Ys_D1s_ss)
print(Ys_D1s)

#   Función de transferencia Y(s)/D2(s) --> Tir
Ys_D2s_ss = co.ss(Aa, Ff2, C, D)
Ys_D2s = co.ss2tf(Ys_D2s_ss)
print(Ys_D2s)

#   Función de transferencia Y(s)/D(s) --> Tci
Ys_D3s_ss = co.ss(Aa, Ff3, C, D)
Ys_D3s = co.ss2tf(Ys_D3s_ss)
print(Ys_D3s)

# Here define the time of Analysis
t_i = 0
t_ch_U = 100
t_ch_d1 = 200
t_ch_d2 = 400
t_ch_d3 = 600

t = np.linspace(t_i, t_ch_U, 1000)
t2 = np.linspace(t_ch_U, t_ch_d1, 1000)
t3 = np.linspace(t_ch_d1, t_ch_d2, 1000)
t4 = np.linspace(t_ch_d2, t_ch_d3, 1000)
t_1, y1 = co.forced_response(Ys_Us, t, 1.01 * U)  # --> Change in Service Flux
t_2, y2 = co.forced_response(Ys_D1s, t2, 1.033 * d1)  # --> Change in C_air
t_3, y3 = co.forced_response(Ys_D2s, t3, 1.02 * d2)  # --> Change in Tir
t_4, y4 = co.forced_response(Ys_D3s, t4, 1.00001 * d3)  # --> Change in Tci

plt.figure(1)
plt.title('INTERCHANGER\nTemperature $[^{\circ}C]$ & Time $[s]$', fontsize=16)
plt.plot(t_1, Tir + y1, label='Step Response from U*0.01 $[T_i]$ work fluid open loop', alpha=0.75)
plt.plot(t_2, Tir + y1[999] + y2, label='Step Response from d1*0.01 $[\.m]$ armor open loop', alpha=0.75)
plt.plot(t_3, Tir + y1[999] + y2[999] + y3, label='Step Response from d2*0.01 $[\.m]$ armor open loop', alpha=0.75)
plt.plot(t_4, Tir + y1[999] + y2[999] + y3[999] + y4, label='Step Response from d3*0.01 $[\.m]$ armor open loop',
         alpha=0.75)
plt.xlabel("Time \t$[s]$ ", fontsize=14)
plt.ylabel("Temperature \t$[^{\circ}C]$ ", fontsize=14)
plt.legend()

# Roots
plt.figure(2)
roots = co.rlocus(Ys_Us)
print("Pole G1", co.damp(Ys_Us))
plt.figure(3)
G_u_Sys_ol = co.series(Ys_Us, Ys_D1s, Ys_D2s, Ys_D3s)
roots2 = co.rlocus(G_u_Sys_ol)
print("Pole G1", co.damp(G_u_Sys_ol))
plt.show()
