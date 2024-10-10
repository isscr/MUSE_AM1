from numpy import array, concatenate
from numpy.linalg import norm

#para comentar varias lineas - ctrl+k
#para descomentar varias lineas - ctrl-u

#lista de valores con array

U = array([1,2,3,4])

print("U=", U[0:2]) #U= [1 2] --> dos valores

#indicar que las dos primeras componentes del vector U pertenecen al vector r y lss otars dos al vecor rdot

r = U[0:2] #se pone un rango hasta 2 porque N=2 (el ultimo) no se incluye en el intervalo de valores
rdot = U[2:4]

F = concatenate((rdot,-r/norm(r)**3), axis=0) #concatenate = para unir dos vectores
                                                #axis = dimension, 0= one-dimensional

print("F=", F)

#U_n+1 = Un + dt*F(Un,tn)