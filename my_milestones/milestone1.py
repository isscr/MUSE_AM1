from numpy import array, zeros, linspace
import matplotlib.pyplot as plt

## DATOS

#Condiciones iniciales
x0 = 1
y0 = 0
vx0 = 0
vy0 = 1

#Parametros temporales
delta_t = 15 # Intervalo de tiempo de la simulacion, s
dt = 0.01 # Paso del tiempo, s

#Numero de pasos
n_pasos = int(delta_t/dt)

## EULER para el Problema de Kepler

#Lista/Array de condiciones iniciales
U0 = [x0,y0,vx0,vy0]

#Listas U_Euler) y F_Euler
U_Euler = zeros((n_pasos,4)) #inicializar un array de ceros para almacenar posiciones y velocidades de columnas [x, y, vx, vy]
U_Euler[0] = U0 #Inicializar el array
F_Euler = zeros(4)

#Calculos del Euler Explicito
for i in range(n_pasos-1):
    #Posicion y velocidad en cada instante de tiempo
    x, y, vx, vy = U_Euler[i]

    #Calculo de las componentes de velocidad de F_Euler
    F_Euler[0] = vx
    F_Euler[1] = vy

    #Calculo de las componentes de aceleracion de F_Euler
    F_Euler[2] = -x/(x**2 + y**2)**(3/2)
    F_Euler[3] = -y/(x**2 + y**2)**(3/2)

    #Actualizar U_Euler
    U_Euler[i + 1, 0] = U_Euler[i][0] + F_Euler[0] * dt
    U_Euler[i + 1, 1] = U_Euler[i][1] + F_Euler[1] * dt 
    U_Euler[i + 1, 2] = U_Euler[i][2] + F_Euler[2] * dt 
    U_Euler[i + 1, 3] = U_Euler[i][3] + F_Euler[3] * dt

#Extraer posiciones
x = U_Euler[:, 0]
y = U_Euler[:, 1]

## RUNGE-KUTTA DE ORDEN 4 para el Problema de Kepler
U_RK4 = zeros((n_pasos,4))
U_RK4[0] = U0
K1 = zeros(4)
K2 = zeros(4)
K3 = zeros(4)
K4 = zeros(4)

for i in range(n_pasos-1):
    x, y, vx, vy = U_RK4[i]

    #Calcular K1
    K1[0] = vx
    K1[1] = vy
    K1[2] = -x/(x**2 + y**2)**(3/2)
    K1[3] = -y/(x**2 + y**2)**(3/2)

    #Calcular K2
    K2[0] = K1[0] + (1/2) * dt * K1[2]
    K2[1] = K1[1] + (1/2) * dt * K1[3]
    K2[2] = K1[2] + (1/2) * dt * K1[2]
    K2[3] = K1[3] + (1/2) * dt * K1[3]

    #Calcular K3
    K3[0] = K1[0] + (1/2) * dt * K2[2]
    K3[1] = K1[1] + (1/2) * dt * K2[3]
    K3[2] = K1[2] + (1/2) * dt * K2[2]
    K3[3] = K1[3] + (1/2) * dt * K2[3]

    #Calcular K4
    K4[0] = K1[0] + (1/2) * dt * K3[2]
    K4[1] = K1[1] + (1/2) * dt * K3[3]
    K4[2] = K1[2] + (1/2) * dt * K3[2]
    K4[3] = K1[3] + (1/2) * dt * K3[3]

    #Actualizar U_RK4
    U_RK4[i + 1, 0] = U_RK4[i, 0] + (dt/6) * (K1[0] + 2*K2[0] + 2*K3[0] + K4[0])
    U_RK4[i + 1, 1] = U_RK4[i, 1] + (dt/6) * (K1[1] + 2*K2[1] + 2*K3[1] + K4[1]) 
    U_RK4[i + 1, 2] = U_RK4[i, 2] + (dt/6) * (K1[2] + 2*K2[2] + 2*K3[2] + K4[2]) 
    U_RK4[i + 1, 3] = U_RK4[i, 3] + (dt/6) * (K1[3] + 2*K2[3] + 2*K3[3] + K4[3])

#Extraer posiciones
x_RK4 = U_RK4[:, 0]
y_RK4 = U_RK4[:, 1]

## CRANK-NICKOLSON para el Problema de Kepler
U_CN = zeros((n_pasos,4))
U_CN[0] = U0
dU_CN = zeros(4)
dU_CN_nueva = zeros(4)

for i in range(n_pasos-1):
    x, y, vx, vy = U_CN[i]

    #Calcular aceleracion
    dU_CN[0] = vx
    dU_CN[1] = vy
    dU_CN[2] = -x/(x**2 + y**2)**(3/2)
    dU_CN[3] = -y/(x**2 + y**2)**(3/2)

    #Actualizar posiciones U_CN
    U_CN[i + 1, 0] = U_CN[i, 0] + dt * dU_CN[0] + (1/2) * dt**2 * dU_CN[2]
    U_CN[i + 1, 1] = U_CN[i, 1] + dt * dU_CN[1] + (1/2) * dt**2 * dU_CN[3] 

    #Calcular nueva aceleracion
    dU_CN_nueva[0] = vx
    dU_CN_nueva[1] = vy
    dU_CN_nueva[2] = -U_CN[i + 1, 0]/(U_CN[i + 1, 0]**2 + U_CN[i + 1, 1]**2)**(3/2)
    dU_CN_nueva[3] = -U_CN[i + 1, 1]/(U_CN[i + 1, 0]**2 + U_CN[i + 1, 1]**2)**(3/2)

    #Actualizar velocidades U_CN
    U_CN[i + 1, 2] = U_CN[i, 2] + (1/2) * dt * (dU_CN[2] + dU_CN_nueva[2]) 
    U_CN[i + 1, 3] = U_CN[i, 3] + (1/2) * dt * (dU_CN[3] + dU_CN_nueva[3])

#Extraer posiciones
x_CN = U_CN[:, 0]
y_CN = U_CN[:, 1]

## REPRESENTACIÓN RESULTADOS EN GRÁFICAS
plt.figure(figsize=(8, 8))
plt.plot(x, y, label="Euler Explícito")
plt.plot(x_RK4, y_RK4, label="Runge-Kutta Orden 4")
plt.plot(x_CN, y_CN, label="Crank-Nickolson")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la EDO para el problema de Kepler empleando diferentes métodos numéricos")
plt.legend()
plt.grid()
plt.show()