 
from numpy import array, zeros, linspace
import matplotlib.pyplot as plt

from miscellaneous import decorators 

@decorators.profiling
def first_version(): 

    U = array( [ 1, 0, 0, 1 ]) #vector de Kepler de 4 componentes
    
    N = 200 #componentes
    x = array( zeros(N) ) #vector de N huecos (ceros)
    y = array( zeros(N) )
    x[0] = U[0] #el primer valor de x (en instante temporal 0) coincide con el primer componente del vector U (1)
    y[0] = U[1]
    
    for i in range(1, N): #range va desde i = 1 hasta i = N - 1 (intervalo cerrado al inicio y abierto al final)
    
      F = array( [ U[2], U[3], -U[0]/(U[0]**2+U[1]**2)**1.5, -U[1]/(U[0]**2+U[1]**2)**1.5 ] ) 
      dt = 0.1 
      U = U + dt * F 
      x[i] = U[0] 
      y[i] = U[1]
    
    plt.plot(x, y)
    plt.show()


@decorators.profiling
def abstraction_for_F(): 

    U = array( [ 1, 0, 0, 1 ])
    
    N = 200 
    x = array( zeros(N) )
    y = array( zeros(N) )
    t = array( zeros(N) )
    x[0] = U[0] 
    y[0] = U[1]
    t[0] = 0 
    
    for i in range(1, N): 

      dt = 0.1 
      t[i] = dt*i
      F = Kepler( U, t[i-1])
    
      U = U + dt * F 
      x[i] = U[0] 
      y[i] = U[1]
    
    plt.plot(x, y)
    plt.show()


@decorators.profiling
def abstraction_for_F_and_Euler(): 

    U = array( [ 1, 0, 0, 1 ])
    
    N = 200 
    x = array( zeros(N) )
    y = array( zeros(N) )
    t = array( zeros(N) )
    x[0] = U[0] 
    y[0] = U[1]
    t[0] = 0 
    
    for i in range(1, N): 

      dt = 0.1 
      t[i] = dt*i
      U = Euler(U, dt, t, Kepler)
      x[i] = U[0] 
      y[i] = U[1]
    
    plt.plot(x, y)
    plt.show()

def Kepler(U, t): 

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = ( x**2  +y**2 )**1.5

    return  array( [ dxdt, dydt, -x/d, -y/d ] ) 

def Euler(U, dt, t, F): 

    return U + dt * F(U, t)



first_version() #hay que llamar a la funciones
abstraction_for_F()
abstraction_for_F_and_Euler()




