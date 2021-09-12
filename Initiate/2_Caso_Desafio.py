# Ejemplo Desafio Fourier
# Por Sergio Andres Diz Ariza  Marzo 04/2021

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

x = np.arange(0,1,0.001)
n = np.arange(1,1000,1)
# Funcion dominio tiempo
ft = x**2 - x + 1/6

# Function
# in dominio Fourier

ff = [(lambda n: np.sum(np.cos(2*np.pi*a*n)/(np.pi**2 * n**2)))(n) for a in x ]

plt.figure()
plt.plot(x,ft, label="Fn domain time")
plt.plot(x,ff, label="Fn domain fourier")
plt.legend()
plt.show()