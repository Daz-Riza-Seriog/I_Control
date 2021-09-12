# Code made for Sergio Andrés Díaz Ariza
# 29 July March 2021
# License MIT
# Introduction to Control: Python Program Linealization interchange

#   Simulador de intercambiador de calor
#   Intercambiador de calor de tubos y coraza
#   Tubos circula el fluído de trabajo (fluido caliente por los tubos)
#   Coraza circula el refrigerante (fluido frio por la coraza).
#   coraza(externo, outter), i: tubos(interno, inner)
#   Objetivo de control consiste en mantener la temperatura del fluido de
#   trabajo en Tcs_sp=326.36 grados kelvin.
#   Octubre 2009

import numpy as np
import matplotlib.pyplot as plt
import sympy
import sympy as sp
from sympy import *
import control as co
import seaborn as sns

sns.set()

# ======================================================================
#    EL INTERCAMBIADOR DESDE LA INGENIERÍA QUÍMICA
# ======================================================================

# Diagrama de flujo de procesos

#           ==================
#  mo,Tfe   =     coraza     =   Tfs
#           ==================
#                   =
#                 = = =
#                =  =  =
#                   =     Q
#                   =
#           ==================
#   mi,Tce  =    tubos       =  Tcs
#           ==================

#   Balances de energía
#   Balance de energía del Benceno
#   roi*Vi*cpi*dx1/dt=mi*cpi*(To-Tf)-UA*(T1-Tf),

#   Balance de energìa del tolueno
#   roo*Vo*cpo*dx2/dt=mo*cpo*(T2-T1)+UA*(T1-Tf),

# ======================================================================
#    EL INTERCAMBIADOR DESDE LA INGENIERÍA DE CONTROL
# ======================================================================

#  Declaración de variables
#   u=mo
#   d=To
#   y=x1=Tf
#   z2=T1

#   Diagrama de bloques

#                         |
#                         | d=To
#                         V
#                 =================
#       u=mo ---> =               = ---> y=x1=Tf
#                 =================
#
#

#   Los balances en el lenguaje de control
#   roi*Vi*cpi*dx1/dt=mi*cpi*(d-x1)-UA*(x2-x1),
#   roo*Vo*cpo*dx2/dt=u*cpo*(T2-x2)+UA*(x2-x1),

#   OBSERVE QUE EL SISTEMA ES NO LINEAL

# ==================================================================
#       PUNTO DE OPERACION
# ==================================================================
mi = 13.88  # flujo masico en los tubos                    [kg/s]
mo = 8.33  # flujo masico en la coraza                    [kg/s]
Tfe = 290.15  # temperatura de entrada del refrigerante      [K]
Tce = 340.15  # temperatura de entrada del fluido de trabajo [K]
Tfs = 313.15  # temperatura de salida del refrigerante       [K]
Tcs = 326.36  # temperatura de salida del fluido de trabajo  [K]
MLTD = ((Tce - Tfs) - (Tcs - Tfe)) / (log((Tce - Tfs) / (Tcs - Tfe)))  # MLTD para contracorriente

# El punto de operación determina los valores de estado estacionario
# Son los valores a los cuales se desea que opere el proceso
uEE = mo
x1EE = Tcs
x2EE = Tfs
dEE = Tce

# %==================================================================
# %   PARAMETROS DEL INTERCAMBIADOR DE TUBOS Y CORAZA
# %==================================================================
# %Parametros geometricos
Nt = 353  # número de tubos
Pt = 0.197  #
dti = 0.01188  # diametro interno de los tubos [m]
Dc = 0.489  # diametro interno de la coraza [m]
espesor = 0.002  # [m]
dto = dti + 2 * espesor  # diametro externo del tubo     [m]

# Parametros relacionados con las propiedades de los fluidos, en ambos casos agua

roi = 979  # densidad del agua del lado de los tubos     [kg/m3]
roo = 999  # densidad del agua del lado de la coraza     [kg/m3]
miui = 420e-6  # viscosidad del fluido de los tubos
miuo = 1080e-6  # viscosidad del fluido de la coraza
cpi = 4188  # capacidad calorifica del fluido de los tubos
cpo = 4184  # capacidad calorifica del fluido de la coraza
ki = 660e-3  # conductividad termica del fluido de los tubos
ko = 598e-3  # conductividad termica del fluido de la coraza
kw = 60  # conductividad termica de la pared del tubo
Toe = 290.15  # es la misma Tfe
Tie = 340.15  # es la misma Tce
Tos = 313.15  # es la misma Tfs
Tis = 326.36  # es la misma Tcs

# ==================================================================
#     CALCULO DEL COEFICIENTE GLOBAL DE TRANSFERENCIA DE CALOR
# ==================================================================
# Calculo del calor transferido
Q = abs(mi * cpi * (Tis - Tie))

# Calculo del coeficiente convectivo para los tubos
Rei = (4 * mi) / (pi * dti * miui)
Pri = (cpi * miui / ki)
f = ((0.790 * log(Rei) - 1.64) ** -2)
hi = (ki / dti) * ((f / 8) * (Rei - 1000) * Pri) / (1 + 12.7 * (f / 8) ** 0.5 * (Pri ** (2 / 3) - 1))

# Calculo para el coeficiete convectivo para la coraza
De = 4 * (0.5 * Pt * 0.86 * Pt - 0.5 * pi * dto ** 2 / 4) / (0.5 * pi * dto)  # diametro equivalente para paso
# triangular
B = 0.4 * Dc  # Recomendacion de diseño (entre 0.4 y 0.6 del diametro interno de la coraza, para corte de baffle del
# 25%)
asf = (Dc * B * (Pt - dto)) / Pt  # area de flujo para la coraza [m^2]
Reo = De * mo / (asf * miuo)
Pro = (cpo * miuo / ko)
JH = 10 ** (0.36897 - 0.8541 * log(Reo) / log(10) + 5.6609e-4 * (log(Reo) / log(10)) ** 2)
ho = (ko / dto) * JH * Reo * Pro ** (1 / 3)

# Calculo para elcoeficiente global de transferencia, tubos cilindricos
U = 1 / (dto / (dti * hi) + dto * log(dto / dti) / (2 * kw) + 1 / ho)

# Calculo del area de transferencia,Vo y Vi
At = Q / (U * MLTD)  # Area de transferencia de calor
Lt = At / (Nt * pi * dto)  # Longitud del IC, en funcion del numero de tubos y el diametro
Vi = Nt * (pi / 4) * dti ** 2 * Lt  # volumen ocupado por el fluido en los tubos
Vo = (pi / 4) * Dc ** 2 * Lt - Nt * (pi / 4) * dto ** 2 * Lt  # volumen ocupado por el fluido en la coraza
# ==========================================================================
#         LINEALIZACION
# ==========================================================================
# Definición variables simbólicas

# se declaran simbolicas
uu = 8.33  # uEE
d = 340.15  # dEE

x1, x2 = sp.symbols('x1,x2', real=True)

f1 = (1 / (roi * Vi * cpi)) * (
            mi * cpi * (d - x1) - U * At * ((d - x2) - (x1 - Tfe)) / sp.log((d - x2) / (x1 - Tfe)))  # f1=dx1/dt
f2 = (1 / (roo * Vo * cpo)) * (
            uu * cpo * (Tfe - x2) + U * At * ((d - x2) - (x1 - Tfe)) / (sp.log((d - x2) / (x1 - Tfe))))  # f2=dx2/dt

A_a = sympy.Matrix([f1, f2]).jacobian([x1, x2])
s = (x1, x2)
A_a = sp.lambdify(s, A_a, "numpy")
x1 = 326.36  # x1EE
x2 = 313.15  # x2EE
Aa = A_a(x1, x2)

uu = sp.symbols('uu', real=True)

f1 = (1 / (roi * Vi * cpi)) * (
            mi * cpi * (d - x1) - U * At * ((d - x2) - (x1 - Tfe)) / (sp.log((d - x2) / (x1 - Tfe))))  # f1=dx1/dt
f2 = (1 / (roo * Vo * cpo)) * (
            uu * cpo * (Tfe - x2) + U * At * ((d - x2) - (x1 - Tfe)) / (sp.log((d - x2) / (x1 - Tfe))))  # f2=dx2/dt

B = sympy.Matrix([f1, f2]).jacobian([uu])
s1 = uu
B_b = sp.lambdify(s1, B, "numpy")
uu = 8.33  # uEE
Bb = B_b(uu)
print(Bb)

d = sp.symbols('d', real=True)

f1 = (1 / (roi * Vi * cpi)) * (
            mi * cpi * (d - x1) - U * At * ((d - x2) - (x1 - Tfe)) / (sp.log((d - x2) / (x1 - Tfe))))  # f1=dx1/dt
f2 = (1 / (roo * Vo * cpo)) * (
            uu * cpo * (Tfe - x2) + U * At * ((d - x2) - (x1 - Tfe)) / (sp.log((d - x2) / (x1 - Tfe))))  # f2=dx2/dt

F = sympy.Matrix([f1, f2]).jacobian([d])
s2 = d
F_f = sp.lambdify(s2, F, "numpy")
d = 340.15  # dEE
Ff = F_f(d)

C = [1, 0]
D = 0

#   Función de transferencia Y(s)/U(s)
Ys_Us_ss = co.ss(Aa, Bb, C, D)
Ys_Us = co.ss2tf(Ys_Us_ss)
print(Ys_Us)

#   Función de transferencia Y(s)/D(s)
Ys_Ds_ss = co.ss(Aa, Ff, C, D)
Ys_Ds = co.ss2tf(Ys_Ds_ss)
print(Ys_Ds)

# Here define the time of Analysis
t_i = 50
t_ch = 150
t_f = 180
step_op = 1

t = np.linspace(t_i, t_ch, 1000)
t_1, y1 = co.forced_response(Ys_Ds, t, 0.01 * d)
t_2, y2 = co.forced_response(Ys_Us, t, 0.1 * uu)


plt.figure(1)
plt.title('INTERCHANGER\nTemperature $[^{\circ}C]$ & Time $[s]$', fontsize=16)
plt.plot(t_1, y1, label='Step Response from d*0.01 $[T_i]$ work fluid open loop', alpha=0.75)
plt.plot(100+t_2, y1[999]+y2, label='Step Response from u*0.1 $[\.m]$ armor open loop', alpha=0.75)
plt.xlabel("Time \t$[s]$ ", fontsize=14)
plt.ylabel("Temperature \t$[^{\circ}C]$ ", fontsize=14)
plt.legend()
plt.show()
