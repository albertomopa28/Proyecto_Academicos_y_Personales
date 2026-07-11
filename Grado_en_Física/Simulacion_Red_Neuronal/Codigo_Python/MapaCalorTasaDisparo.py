import numpy as np
import matplotlib.pyplot as plt

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

# Realizar un paso del método de Euler
def rk4_step(r, t, dt, matriz, g, I):
    k1 = drdt(r, t, matriz, g, I)
    k2 = drdt(r + 0.5 * k1 * dt, t + 0.5 * dt, matriz, g, I)
    k3 = drdt(r + 0.5 * k2 * dt, t + 0.5 * dt, matriz, g, I)
    k4 = drdt(r + k3 * dt, t + dt, matriz, g, I)
    return r + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

# Calcular la covarianza 
def calcular_covarianza(r, N):
    media_r = np.mean(r)
    suma = np.sum(((r - media_r)**2))
    covarianza = np.sqrt(suma /N)
    return covarianza

# Función para calcular el módulo del autovalor máximo de una matriz dada
def calcular_modulo_autovalor_maximo(W):
    autovalores = np.linalg.eigvals(W)
    modulo_maximo = np.max(np.real(autovalores))
    return modulo_maximo

# Configuración de parámetros
N = 1000
p=0.2
filas = N
columnas = N
I = 0
timesteps = 300
dt = 0.2
M = 3  # Número de simulaciones a realizar para cada combinación de g y w

g_values = np.arange(0, 0.15, 0.00375)
w_values = np.arange(0, 6, 0.15)
mean_results = np.zeros((len(g_values), len(w_values)))
covariance_results = np.zeros((len(g_values), len(w_values)))
lambda_max = np.zeros((len(g_values), len(w_values)))
# Realizar simulaciones para cada combinación de g y w
for g_index, g in enumerate(g_values):
    print('g=', g)
    for w_index, w in enumerate(w_values):
        print('w=',w)
        # Inicializar acumuladores para medias y varianzas
        S = 0
        covariance_accumulator = 0
        
        for _ in range(M):
            matriz = generar_matriz_W1y2(filas, columnas, w, g)
            #Calculamos y almacenamos el autovalor max
            lambda_max[g_index, w_index]+=calcular_modulo_autovalor_maximo(matriz)/M
            
            r = np.random.rand(N)  # Condición inicial de r random
            t = 0.0
            spatial_mean_r_squared=0
            spatial_covariance_r=0
            for _ in range(timesteps):
                r = rk4_step(r, t, dt, matriz, g, I)
                t += dt
                S +=np.dot(r,r) #Calculamos S a cada paso, acumulamos sum_i(r_i^2) que es el prod vect de r*r
                covariance_accumulator += calcular_covarianza(r, N) #igual con la cov

        # Guardar los resultados promedio
        mean_results[g_index, w_index] = S/(M*timesteps*N)
        covariance_results[g_index, w_index] = covariance_accumulator/(M*timesteps)
        
# Definimos las funciones f(x) y g(x)
def f(w):
    return 2/(N*p*(0.8-0.2*w))

def g(w):
    return np.sqrt(1/(N*(p/3.-p**2/4)*(0.8+0.2*w**2)))

#RA + Sparsity
def g2(w):
    return 1./np.sqrt(N*p*((1-0.2)/12+0.2*(w)**2/12))

def f2(w):
    return 2/(N*(0.8-0.2*w))

# Generamos los valores de x 
#Para RA
#w1 = np.linspace(0, 3.82325, 400)
#w2 = np.linspace(3.82325, 6, 400)

#Para RA + Sparsity
w1 = np.linspace(0, 3.6222, 400)
w2 = np.linspace(3.6222, 6, 400)

#Con Sparsity
#w1 = np.linspace(0, 3.344, 400)
#w2= np.linspace(3.344 , 6, 400)

# Calculamos los valores de f(x) y g(x)
y_f = f2(w1)
y_g = g2(w2)

# Ajustar los pasos de los ejes para una mejor visualización
#El tercer argumento es el paso
w_ticks = np.arange(0, 6.001, 1)
g_ticks = np.arange(0, 0.15, 0.02)

# Plot de la media temporal de las medias espaciales
plt.figure(figsize=(10, 8))
plt.imshow(mean_results, cmap='inferno', aspect='auto', extent=[0, 6, 0, 0.15], origin='lower')
plt.colorbar(label='S(g,w)')
plt.xticks(w_ticks)
plt.yticks(g_ticks)
plt.xlabel('w')
plt.ylabel('g')
plt.title('S(g,w)')

# Superponer las funciones f(w) y g(w)
plt.plot(w1, y_f, color='blue', label=r'\lambda(g,w)=1', linewidth=4)
plt.plot(w2, y_g, color='red', label=r'R(g,w)=1', linewidth=4)
plt.legend()

plt.show()

# Plot de la media temporal de las varianzas espaciales
plt.figure(figsize=(10, 8))
plt.imshow(covariance_results, cmap='inferno', aspect='auto', extent=[0, 6, 0, 0.15], origin='lower')
plt.colorbar(label='Covarianzas')
plt.xticks(w_ticks)
plt.yticks(g_ticks)
plt.xlabel('w')
plt.ylabel('g')
plt.title('Covarianzas')

# Superponer las funciones f(w) y g(w)
plt.plot(w1, y_f, color='blue', label=r'$\lambda(g,w)=1$', linewidth=4)
plt.plot(w2, y_g, color='red', label=r'R(g,w)=1', linewidth=4)
plt.legend()
plt.show()

# Crear una figura
plt.figure(figsize=(10, 8))

# Mostrar el mapa de calor
plt.imshow(lambda_max, cmap='inferno', aspect='auto', extent=[0, 6, 0, 0.15], origin='lower')
plt.colorbar(label=r'$\lambda_{max}$')

# Superponer los puntos en verde fosforito para lambda_max entre 0.7 y 1.3
lambda_mask = (lambda_max >= 0.8) & (lambda_max <= 1.2)
X, Y = np.meshgrid(np.linspace(0, 6, lambda_max.shape[1]), np.linspace(0, 0.15, lambda_max.shape[0]))

plt.scatter(X[lambda_mask], Y[lambda_mask], color='lime', marker='o', s=10, label=r'$0.8 \leq \lambda_{max} \leq 1.2$')

# Superponer las funciones f(w) y g(w)
plt.plot(w1, y_f, color='blue', label=r'$\lambda(g,w)=1$', linewidth=4)
plt.plot(w2, y_g, color='red', label=r'R(g,w)=1$', linewidth=4)

# Configuración de la gráfica
plt.xticks(w_ticks)
plt.yticks(g_ticks)
plt.xlabel('w')
plt.ylabel('g')
plt.title('Módulo máximo de los autovalores')
plt.legend()
plt.show()