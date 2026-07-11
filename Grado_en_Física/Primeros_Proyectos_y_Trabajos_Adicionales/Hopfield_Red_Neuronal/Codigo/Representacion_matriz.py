import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Configuración de parámetros
N = 10  # Tamaño de la matriz
mu = 3  # Número de patrones
T = 0.03  # Temperatura del sistema
iteraciones = 100 * N * N
paso = N * N
semilla = 1934050348
np.random.seed(semilla)

# Inicializar matrices y parámetros
S = np.random.choice([0, 1], size=(N, N))  # Matriz inicial aleatoria
epsilon = np.random.choice([0, 1], size=(mu, N, N))  # Patrones aleatorios
a = np.zeros(mu)
theta = np.zeros((N, N))
solapamiento = np.zeros(mu)

# Funciones auxiliares
def calcular_a(epsilon):
    return np.mean(epsilon, axis=(1, 2))

def calcular_omega(epsilon, a):
    omega = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    omega[i, j, k, l] = np.sum((epsilon[:, i, j] - a) * (epsilon[:, k, l] - a)) / mu
    return omega

def calcular_theta(omega):
    return np.sum(omega, axis=(2, 3))

def calc_DeltaE(S, omega, theta, n, m):
    return -2 * S[n, m] * (theta[n, m] + np.sum(S * omega[n, m]))

def calcularp(S, T, DeltaE, n, m):
    return min(1, np.exp(-DeltaE / T))

def calcular_solapamiento(S, epsilon, a):
    return np.mean((epsilon - a[:, None, None]) * (S - a[:, None, None]), axis=(1, 2))

# Inicializar parámetros
a = calcular_a(epsilon)
omega = calcular_omega(epsilon, a)
theta = calcular_theta(omega)

# Crear listas para almacenar las matrices en el tiempo
matrices_tiempo = []

# Iterar sobre el sistema
for i in range(paso, iteraciones + 1):
    n, m = np.random.randint(0, N, size=2)
    DeltaE = calc_DeltaE(S, omega, theta, n, m)
    p = calcularp(S, T, DeltaE, n, m)

    if np.random.rand() < p:
        S[n, m] = 1 - S[n, m]

    if i % paso == 0:
        matrices_tiempo.append(S.copy())

# Crear animación
def actualizar_imagen(frame):
    matriz.set_array(matrices_tiempo[frame])
    return [matriz]

fig, ax = plt.subplots()
matriz = ax.imshow(matrices_tiempo[0], cmap="binary", interpolation="none")
plt.colorbar(matriz)
ani = animation.FuncAnimation(fig, actualizar_imagen, frames=len(matrices_tiempo), blit=True)

# Guardar el video
output_video = "hopfield_simulation.mp4"
ani.save(output_video, fps=10, writer="ffmpeg")

print(f"Simulación completada. Video guardado como {output_video}.")
