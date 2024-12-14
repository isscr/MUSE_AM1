import numpy as np

V = np.array([1,2,3])

pV=V
pV[0]=4

print(V)

U = V.copy
#U[0]=18

print(U)
