# Validez de un modelo lineal
# Sistema de almacenamiento
# Por Sergio Andres Diz Ariza  Marzo 04/2021
# del codigo de Prof Lina Gómez febrero de 2010 en Matlab

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

###############################
#           Parametros        #
###############################

Atanque = 0.1963 # [m2]
g = 9.8 #[m/s]
r = np.sqrt(Atanque/np.pi) #[m]
Diametro = 2*r  #[m]

###############################
#       Estado estacionario   #
###############################

uss=0.00049 #[m2] Cvee
dss=0.001534 #[m/s] Foee
yss=1.1 #[m]

#################################
#   Compara Modelo no lineal    #
#   con modelo lineal           #
#################################

y = np.arange(0,2,0.1)
y1 = np.arange(-1.1,0.9,0.1)  #variable de desviación y1=deltay=y-yss
fnl = 1/Atanque*(dss-np.sqrt(g)*uss*y**(1/2))
fl = 1/Atanque*(dss - 1/2*np.sqrt(g)*uss*yss**(-1/2)*y1 - np.sqrt(g)*yss**(1/2)*uss) #está en variable de desviación
plt.figure(1)
plt.plot(y,fnl)
plt.plot(y,fl)
uss=0.00049 #[m2]
dss=0.001534 #[m/s]
yss=1.1 #[m]

#################################
#   Compara Modelo no lineal    #
#   con modelo lineal           #
#################################

y = np.arange(0,2,0.1)
y1 = np.arange(-1.1,0.9,0.1) #variable de desviación y1=deltay=y-yss
fnl = 1/Atanque*(dss-np.sqrt(g)*uss*y**(1/2))
fl = 1/Atanque*(dss - 1/2*np.sqrt(g)*uss*yss**(-1/2)*y1 - np.sqrt(g)*yss**(1/2)*uss) #está en variable de desviación
plt.figure (1)
plt.plot(y,fnl)
plt.plot(y,fl)
plt.title("Nivel del Tanque")
plt.xlabel("Label X")
plt.ylabel("Label Y")

plt.show()