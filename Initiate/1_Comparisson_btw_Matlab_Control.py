import numpy as np

A = np.array([[1, 2, 3],[4, 5 ,6],[7, 8 ,9]])
B = np.array([[4,5,6],[7,8,9],[1,2,3]])

print("\nSuma de las Matrices componenete a componente:\n",A+B)
print("\nSuma de las Matrices:\t",np.sum(A+B))

# Hagamos un vector lineal de 100 digitos
X = np.arange(0,100)
print("\nVectopr X escrito uniformemente:\n",X)
print("\nSuma del Vector:\t",np.sum(X))



