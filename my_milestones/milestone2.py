from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
from scipy.optimize import newton

## FUNCIONES

# CONDICIONES INICIALES
def CI():
    x0 = 1
    y0 = 0
    vx0 = 0
    vy0 = 1

    U0 = array([x0, y0, vx0, vy0]) 
    return U0

# KEPLER
def Kepler_acceleration(x, y):
    "Calcula la aceleración según las leyes de Kepler."
    r = (x**2 + y**2)**(3/2)
    ax = -x / r
    ay = -y / r
    
    return ax, ay

def Kepler_derivatives(U,t):
    "Calcula las derivadas para el problema de Kepler."
    x, y, vx, vy = U

    ax, ay = Kepler_acceleration(x, y)

    return  array( [vx, vy, ax, ay] ) # Devuelve F(U,t)

## PROBLEMA DE CAUCHY
def Cauchy_problem(F, t, U0, scheme):  
    """Resuelve un problema de Cauchy (de condiciones iniciales) genérico dado un esquema numérico temporal.
    
    Parametros:
    F: Función que define el sistema de ecuaciones diferenciales dU/dt = F(U, t).
    t: Array de puntos en el tiempo donde se calculará la solución.
    U0: Condiciones iniciales del sistema.
    scheme: Esquema numérico que actualiza la solución en cada paso.

    Devuelve:
    U: Matriz donde cada fila representa el estado del sistema en un instante de tiempo.
    """
    Np = len(t) - 1  # Número de pasos 
    Nv = len(U0)    # Número de variables del sistema
    
    "len() devuelve la longitud/número de elementos de un objeto"

    # Matriz para almacenar las soluciones
    U = zeros((Np + 1, Nv))
    U[0, :] = U0  # Condiciones iniciales
    
    # Iteración para avanzar en el tiempo
    for i in range(Np):
        dt = t[i+1] - t[i]  # Tamaño del paso temporal
        U[i+1, :] = scheme(U[i, :], dt, t[i], F)

    return U

## ESQUEMAS NUMÉRICOS

# ESQUEMA EULER
def Euler(U, dt, t, F):
    """Esquema de Euler explícito."""
    return U + dt * F(U, t)

# ESQUEMA RUNGE-KUTTA ORDEN 4
def RK4(U, dt, t, F):
    """Esquema de Runge-Kutta de orden 4."""
    K1 = F(U, t)
    K2 = F(U + (1/2) * dt * K1, t + (1/2) * dt)
    K3 = F(U + (1/2) * dt * K2, t + (1/2) * dt)
    K4 = F(U + (1/2) * dt * K3, t + (1/2) * dt)
    return U + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4)

# ESQUEMA CRANK-NICKOLSON
def CN(U, dt, t, F): 
    """Esquema de Crank-Nickolson implícito."""
    def G(X): # Siendo X = U(n-1)
        return X - U - dt/2 * (F(X, t) + F(U, t))
    return newton(G, U) #U no se modifica al "pasar" por la función, la U que sale es la que devuelve Newton

# Esquema implícito EULER IMPLÍCITO
def Euler_Inverso(U, dt, t, F): 
    """Esquema de Euler implícito."""
    def G(X): # Siendo X = U(n-1)
        return X - U - dt * F(X, t)
    return newton(G, U) 

## MAIN

# Parámetros del problema
tf = 10 # Tiempo final
N = 1000 # Numero de intervalos de tiempo
t = linspace(0, tf, N)  # Intervalo de tiempo

# Condiciones iniciales
U0 = CI()

# Resolver problema de Cauchy con diferentes esquemas
U_Euler = Cauchy_problem(Kepler_derivatives, t, U0, Euler)

U_RK4 = Cauchy_problem(Kepler_derivatives, t, U0, RK4)

U_CN = Cauchy_problem(Kepler_derivatives, t, U0, CN)

U_EulerInv = Cauchy_problem(Kepler_derivatives, t, U0, Euler_Inverso)

# Extraer posiciones
x_E = U_Euler[:, 0]
y_E = U_Euler[:, 1]

x_RK4 = U_RK4[:, 0]
y_RK4 = U_RK4[:, 1]

x_CN = U_CN[:, 0]
y_CN = U_CN[:, 1]

x_EInv = U_EulerInv[:, 0]
y_EInv = U_EulerInv[:, 1]

# REPRESENTACION RESULTADOS EN GRAFICAS
plt.figure(figsize=(10, 8))
plt.plot(x_E, y_E, label="Euler Explícito")
plt.plot(x_RK4, y_RK4, label="Runge-Kutta Orden 4")
plt.plot(x_CN, y_CN, label="Crank-Nickolson")
plt.plot(x_EInv, y_EInv, label="Euler Inverso")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la EDO para el problema de Kepler empleando diferentes métodos numéricos")
plt.legend()
plt.grid()
plt.show()