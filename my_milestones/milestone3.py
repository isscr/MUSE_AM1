from numpy import array, zeros, linspace, log10, polyfit
import matplotlib.pyplot as plt
from scipy.optimize import newton
from numpy.linalg import norm

from milestone2 import CI, Kepler_acceleration, Kepler_derivatives, Cauchy_problem, Euler, Euler_Inverso, RK4, CN

## FUNCIONES

# REFINACIÓN DE MALLA
def refine_mesh(t1): 
    '''''''''''
    Refinación de malla: dada la partición t1 (con N+1 puntos), obtiene la partición t2 (que tiene 2N+1 puntos)
    Los nodos pares de t2 seguirán siendo los mismos que los de t1, y los nodos impares serán los puntos medios de los intervalos de t1
    
    INPUTS:
     - t1: partición temporal con N+1 puntos
    '''''''''''
    N = len(t1) - 1  
    t2 = zeros(2*N +1) 

    for i in range(0,N): 
           t2[2*i] =  t1[i] #Nodos pares
           t2[2*i+1] = ( t1[i]  + t1[i+1] )/ 2 #Nodos impares
    t2[2*N] = t1[N]  #Se añade el último nodo, hay que añadirlo a mano.

    return t2 

# ERROR DE CONVERGENCIA DE LOS ESQUEMAS NUMÉRICOS: Extrapolación de Richardson
def schemes_error(U0, t, F, problem, scheme):
    '''''''''''
    INPUTS:
        - U0: vector de condiciones iniciales
        - F(U,t): función a resolver
        - scheme(U, F, t): función que representa el esquema numérico a utilizar
        - problem(Esquema, F, U0, t): Función que representa el problema a resolver (Cauchy hasta el momento)
        - t: partición temporal 
    '''''''''''

    t1 = t
    t2 = refine_mesh(t1)

    U_1 = problem(F, t1, U0, scheme) # Solución del problema (en este caso el de Cauchy) con la malla original
    U_2 = problem(F, t2, U0, scheme) # Solución del problema (en este caso el de Cauchy) con la malla modificada (más fina), es decir, con una malla refinada
    
    N = len(t) - 1
    Nv = len(U0)
    Error = zeros((N+1, len(U0)))
    
    # Para calcular el error se hace la resta, pero un vector no se puede restar de otro si uno mide N+1 y el otro N, por eso se hace la resta en los nodos pares
    for i in range (0, N+1): 
        Error[i, :] = U_2[2*i, :] - U_1[i, :] # en este caso los dos puntitos significan: para todas las variables
        
    return U_1, Error

# CONVERGENCIA DE LOS ESQUEMAS NUMÉRICOS
def convergence_rate(t, F, U0, error, problem, scheme):
    '''''''''''
    INPUTS:
        - U0: Vector del estado inicial
        - F: Función a resolver
        - error(U0, F, Problema, Esquema, t): Función que devuelve un vector con el error de un esquema en cada paso temporal
        - scheme: Esquema temporal a resolver
        - t: partición temporal
    OUTPUTS:
        - logN:  vector for the different number of time partitions 
        - logE:  error for each time partition      
        - order:  order of Error of the temporal scheme 
    '''''''''''
    
    N = len(t) - 1
    np = 6  #Número de puntos de la regresión, si se sube más tarda MUCHO en converger
    logE = zeros(np)
    logN = zeros(np)

    t1 = t
    for i in range(0,np):
        
        N = len(t1) - 1
        U, Error = error(U0, t1, F, problem, scheme)  #Asumiendo que Error devuelve U_1 y Error
        logE[i] = log10(norm(Error[N, :])) 
        logN[i] = log10(float(N))
        
        t1 = refine_mesh(t1) 
        
    y = logE[ logE > -12 ]
    x = logN[ 0:len(y) ]
    order, b = polyfit(x, y, 1) # Regresión lineal para encontrar la pendiente de la recta que mejor se ajusta a los datos
    
    # print("Order =", Order, "b =", b)
    
    logE = logE - log10( 1 - 1./2**abs(order) )

    return logN, logE, order

## MAIN

# Parámetros del problema
tf = 10 # Tiempo final
N = 1000 # Numero de intervalos de tiempo
t = linspace(0, tf, N)  # Intervalo de tiempo

# Condiciones iniciales
U0 = CI()

# Error
U_E, Error_E = schemes_error(U0, t, Kepler_derivatives, Cauchy_problem, Euler)
U_RK4, Error_RK4 = schemes_error(U0, t, Kepler_derivatives, Cauchy_problem, RK4)
U_CN, Error_CN = schemes_error(U0, t, Kepler_derivatives, Cauchy_problem, CN)
U_EI, Error_EI = schemes_error(U0, t, Kepler_derivatives, Cauchy_problem, Euler_Inverso)

# Convergencia
logN_E, logE_E, Order_E = convergence_rate(t, Kepler_derivatives, U0, schemes_error, Cauchy_problem, Euler)
logN_RK4, logE_RK4, Order_RK4 = convergence_rate(t, Kepler_derivatives, U0, schemes_error, Cauchy_problem, RK4)
logN_CN, logE_CN, Order_CN = convergence_rate(t, Kepler_derivatives, U0, schemes_error, Cauchy_problem, CN)
logN_EI, logE_EI, Order_EI = convergence_rate(t, Kepler_derivatives, U0, schemes_error, Cauchy_problem, Euler_Inverso)
print ("Order Euler =", Order_E)
print ("Order RK4 =", Order_RK4)
print ("Order CN =", Order_CN)
print ("Order EI =", Order_EI)

# GRAFICAS
# Gráfica de las todas soluciones
plt.plot(t, U_E[:, 0], label="Euler")
plt.plot(t, U_RK4[:, 0], label="RK4")
plt.plot(t, U_CN[:, 0], label="Crank-Nickolson")
plt.plot(t, U_EI[:, 0], label="Euler Inverso")
plt.plot(t, Error_E[:, 0],  label="Error Euler")
plt.plot(t, Error_RK4[:, 0],  label="Error RK4")
plt.plot(t, Error_CN[:, 0],  label="Error CN")
plt.plot(t, Error_EI[:, 0],  label="Error EI")
plt.title("Soluciones de todos los esquemas numéricos y sus errores")
plt.legend()
plt.xlabel("t")
plt.show()

# Gráfica de Euler y su error
plt.plot(t, U_E[:, 0], label="Euler")
plt.plot(t, Error_E[:, 0],  label="Error Euler")
plt.title("Solución del esquema de Euler y su error")
plt.legend()
plt.xlabel("t")
plt.show()

# Gráfica de RK4 y su error
plt.plot(t, U_RK4[:, 0], label="RK4")
plt.plot(t, Error_RK4[:, 0],  label="Error RK4")
plt.title("Solución del esquema de RK4 y su error")
plt.legend()
plt.xlabel("t")
plt.show()

# Gráfica de Crank-Nickolson y su error
plt.plot(t, U_CN[:, 0], label="Crank-Nickolson")
plt.plot(t, Error_CN[:, 0],  label="Error CN")
plt.title("Solución del esquema de Crank-Nickolson y su error")
plt.legend()
plt.xlabel("t")
plt.show()

# Gráfica de Euler Inverso y su error
plt.plot(t, U_EI[:, 0], label="Euler Inverso")
plt.plot(t, Error_EI[:, 0],  label="Error EI")
plt.title("Solución del esquema de Euler Inverso y su error")
plt.legend()
plt.xlabel("t")
plt.show()

# Gráfica de convergencia de todos los esquemas
plt.axis('equal') # Cada unidad en el eje 𝑥 x tiene la misma escala visual que cada unidad en el eje y
plt.xlabel('logN')
plt.ylabel('logE')
plt.plot(logN_E, logE_E, '-b')
plt.plot(logN_RK4, logE_RK4, '-g')
plt.plot(logN_CN, logE_CN, '-r')
plt.plot(logN_EI, logE_EI, '-m')
plt.title("Convergencia de todos los esquemas numéricos")
plt.show()

# Gráfica de convergencia de Euler
plt.axis('equal') 
plt.xlabel('logN_E')
plt.ylabel('logE_E')
plt.plot(logN_E, logE_E, '-b')
plt.title("Convergencia del esquema de Euler")
plt.show()

# Gráfica de convergencia de RK4
plt.axis('equal') 
plt.xlabel('logN_RK4')
plt.ylabel('logE_RK4')  
plt.plot(logN_RK4, logE_RK4, '-g')
plt.title("Convergencia del esquema de RK4")
plt.show()

# Gráfica de convergencia de Cranck-Nickolson
plt.axis('equal') 
plt.xlabel('logN_CN')
plt.ylabel('logE_CN')
plt.plot(logN_CN, logE_CN, '-r')
plt.title("Convergencia del esquema de Crank-Nickolson")
plt.show()

# Gráfica de convergencia de Euler Inverso
plt.axis('equal') 
plt.xlabel('logN_EI')
plt.ylabel('logE_EI')
plt.plot(logN_E, logE_E, '-m')
plt.title("Convergencia del esquema de Euler Inverso")
plt.show()