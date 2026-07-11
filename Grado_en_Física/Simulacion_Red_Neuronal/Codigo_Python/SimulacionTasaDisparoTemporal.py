import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Generar la matriz W

def generar_matriz_W(filas, columnas, w, g):
    matriz_J = np.zeros((filas, columnas))  # Inicializar matriz con cero
    matriz_M = np.zeros((filas, columnas)) 
    matriz_W = np.zeros((filas, columnas)) 
    
    # Llenar el 80% de las columnas con valores aleatorios entre 0 y 1
    for i in range(columnas):
        if i < 0.8 * columnas:
            matriz_J[:, i] = np.random.uniform(-g*1./2., g*1./2., filas)
            matriz_M[:, i] = g*1./2.
            
        else:
            # Llenar el 20% restante con valores aleatorios entre 0 y w
            matriz_J[:, i] = np.random.uniform(-g*w/2., g*w/2., filas)
            matriz_M[:, i] = -g*w/2.
    
    for i in range(len(matriz_J)):
        I  = np.sum(matriz_J[i, :])/N
        for j in range(len(matriz_J)):
            matriz_W[i][j] = matriz_J[i][j] + matriz_M[i][j] - I
    
    return matriz_W


# Condicion solo de sparsity
def generar_matriz_W2(filas, columnas, w, g):
    matriz_W = np.zeros((filas, columnas))
    
    # Determinar el punto de corte para las columnas
    cutoff = int(0.8 * columnas)
    
    # Llenar las columnas
    matriz_W[:, :cutoff] = np.random.uniform(0., g, (filas, cutoff))
    matriz_W[:, cutoff:] = np.random.uniform(-g*w, 0., (filas, columnas - cutoff))
    
    return condicion1(filas, columnas, matriz_W)


def generar_matriz_W1y2(filas, columnas, w, g):
    matriz_J = np.zeros((filas, columnas))  # Inicializar matriz con cero
    matriz_M = np.zeros((filas, columnas)) 
    matriz_W = np.zeros((filas, columnas)) 
    
    # Llenar el 80% de las columnas con valores aleatorios entre 0 y 1
    for i in range(columnas):
        if i < 0.8 * columnas:
            matriz_J[:, i] = np.random.uniform(-g*1./2., g*1./2., filas)
            matriz_M[:, i] = g*1./2.
            
        else:
            # Llenar el 20% restante con valores aleatorios entre 0 y w
            matriz_J[:, i] = np.random.uniform(-g*w/2., g*w/2., filas)
            matriz_M[:, i] = -g*w/2.

    # Aplicar la condición de sparsity sincronizada a ambas matrices
    condicion1(filas, columnas, matriz_J)
    
    # Sumar las matrices J y M para obtener W
    for i in range(len(matriz_J)):
        I  = np.sum(matriz_J[i, :])/N
        for j in range(len(matriz_J[0])):
            matriz_W[i][j] = matriz_J[i][j] + matriz_M[i][j] - I  

    return matriz_W

def condicion1(filas, columnas, matriz):
    # Determinar el número de elementos a convertir en cero
    num_zeros = int(0.8 * filas * columnas)
    
    # Generar posiciones aleatorias para establecer en cero
    indices = np.random.choice(filas * columnas, num_zeros, replace=False)
    
    # Convertir las posiciones a índices de fila y columna
    filas_indices, columnas_indices = np.unravel_index(indices, (filas, columnas))
    
    # Establecer en cero los elementos seleccionados
    matriz[filas_indices, columnas_indices] = 0
    
    return matriz

def condicion2(matriz):
    # Calcular la media de cada fila
    media_filas = np.mean(matriz, axis=1)
    
    # Restar la media de cada fila
    matriz_resultante = matriz - media_filas[:, np.newaxis]
    
    return matriz_resultante



def condicion1y2(filas, columnas, matriz1, matriz2): 
    # Determinar el número de elementos a convertir en cero
    num_zeros = int(0.8 * filas * columnas)
    
    # Generar posiciones aleatorias para establecer en cero
    indices = np.random.choice(filas * columnas, num_zeros, replace=False)
    
    # Convertir las posiciones a índices de fila y columna
    filas_indices, columnas_indices = np.unravel_index(indices, (filas, columnas))
    
    # Establecer en cero los elementos seleccionados
    matriz1[filas_indices, columnas_indices] = 0
    matriz2[filas_indices, columnas_indices] = 0
  
   # Convertir la lista a un array de numpy si no lo es
    matriz = np.array(matriz1)
    
    # Crear una copia de la matriz para trabajar en ella
    matriz_resultante = np.copy(matriz)
    
    # Iterar sobre cada fila de la matriz
    for i in range(matriz.shape[0]):
        # Extraer los elementos no cero de la fila
        fila = matriz[i]
        elementos_no_cero = fila[fila != 0]
        
        # Calcular la media de los elementos no cero si existen elementos no cero
        if len(elementos_no_cero) > 0:
            media_no_cero = np.mean(elementos_no_cero)
        else:
            media_no_cero = 0
        
        # Restar la media de los elementos no cero a los elementos no cero en la fila
        matriz_resultante[i] = np.where(fila != 0, fila - media_no_cero, 0)
    
    return matriz_resultante

# Definir la función epsilon (siempre devuelve 0 en este caso)
def epsilon(t):
    return 0

# Definir la derivada de r con respecto al tiempo
def drdt(r, t, matriz, g, I):
    sum_term = (matriz) @ r
    return -r + np.tanh(sum_term + I) + epsilon(t)

# Realizar un paso del método Runge-Kutta de cuarto orden (RK4)
def rk4_step(r, t, dt, matriz, g, I):
    k1 = drdt(r, t, matriz, g, I)
    k2 = drdt(r + 0.5 * k1 * dt, t + 0.5 * dt, matriz, g, I)
    k3 = drdt(r + 0.5 * k2 * dt, t + 0.5 * dt, matriz, g, I)
    k4 = drdt(r + k3 * dt, t + dt, matriz, g, I)
    return r + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

# Calcular la covarianza 
def calcular_covarianza(r):
    N = len(r)
    media_r = np.mean(r)
    suma = np.sum(((r - media_r)**2))
    covarianza = np.sqrt(suma /N)
    return covarianza

N = 1000
filas = N
columnas = N
I = 0
timesteps = 2000
dt = 0.01
w=6
g=0.15

matriz = generar_matriz_W1y2(filas, columnas, w, g)
r = np.random.rand(N)  # Condición inicial de r aleatoria entre 0 y 1 para cada r_i  
t = 0.0
r_results = []
time_points = []
spatial_means = []  # Lista para almacenar las medias espaciales
Cov = []  # Lista para almacenar las covarianzas
r_medio=0
# Integrar en el tiempo
for _ in range(timesteps):
    if _ % 1000 == 0:
        print(_)
    r_medio+=np.mean(r)/timesteps
    r_results.append(r.copy())
    time_points.append(t)
    spatial_means.append(np.mean(r))  # Calcular y almacenar la media espacial en cada paso
    Cov.append(calcular_covarianza(r))  # Calcular y almacenar la covarianza en cada paso
    r = rk4_step(r, t, dt, matriz, g, I)
    t += dt
print('r_media:', r_medio)
# Convertir los resultados en arrays de numpy para un análisis más fácil
r_results = np.array(r_results)
time_points = np.array(time_points)
spatial_means = np.array(spatial_means)
Cov = np.array(Cov)

# Seleccionar aleatoriamente 20 trayectorias
selected_indices = np.random.choice(N, size=200, replace=False)
selected_r_values = r_results[:, selected_indices]

# Guardar los resultados en un archivo
with open('resultados.txt', 'w') as file:
    file.write(f"{'Tiempo':<10}\t{'Media Espacial':<15}\t{'Covarianza':<15}\t")
    file.write("\t".join([f"r{i}" for i in selected_indices]))
    file.write("\n")
    for i in range(timesteps):
        file.write(f"{time_points[i]:<10.3f}\t{spatial_means[i]:<15.3f}\t{Cov[i]:<15.3f}\t")
        file.write("\t".join([f"{x:.3f}" for x in r_results[i, selected_indices]]))
        file.write("\n")

# Visualizar los resultados usando matplotlib
plt.figure(figsize=(10, 6))
for i in range(selected_r_values.shape[1]):
    plt.plot(time_points, selected_r_values[:, i], linestyle='--')  # Plot de las trayectorias seleccionadas
plt.plot(time_points, spatial_means, label='Media Espacial', color='red', linewidth=4)  # Plot de la media espacial
plt.plot(time_points, Cov, label='Covarianza', color='black', linewidth=4)  # Plot de la covarianza
r_medio = [r_medio] * len(time_points)
plt.plot(time_points, r_medio, label='r_{medio}', color='orange', linewidth=4)
plt.xlabel('t(s)')
plt.ylabel('$r(t)$')
plt.legend(loc='upper right')
plt.title('Evolución de $r$ en el tiempo')
plt.show()