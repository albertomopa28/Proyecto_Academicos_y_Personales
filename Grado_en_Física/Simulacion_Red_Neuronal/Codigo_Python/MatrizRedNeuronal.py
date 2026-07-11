import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#La idea es generar una matriz W=J+M. Tenemos que generar primero una matriz J, que debe cumplir la condición de que la suma de sus 
#filas sea 0. Esto se hará de la siguiente manera: se generan los elementos de la matriz J, para cada fila, se suma para obtener el 
#valor total. Tras esto, se hace la suma de las filas y se resta a cada elemento suma/N.


#Primero, tenemos que hacer una función que cree la matriz aleatoria J, esto es, tiene que crear una matriz con solo el 20% de elementos
#distintos de 0, y que entre ese 20%, cada columna tiene una distribución uniforme diferente
def generar_matriz_W(filas, columnas, w):
    matriz_J = np.zeros((filas, columnas))  # Inicializar matriz con cero
    matriz_M = np.zeros((filas, columnas)) 
    matriz_W = np.zeros((filas, columnas)) 
    
    # Llenar el 80% de las columnas con valores aleatorios entre 0 y 1
    for i in range(columnas):
        if i < 0.8 * columnas:
            matriz_J[:, i] = np.random.uniform(-1./2., 1./2., filas)
            matriz_M[:, i] = 1./2.
            
        else:
            # Llenar el 20% restante con valores aleatorios entre 0 y w
            matriz_J[:, i] = np.random.uniform(-1.*w/2., 1.*w/2., filas)
            matriz_M[:, i] = -1.*w/2.
    
    #condicion2(matriz_J)
    
    #autovalores=calcular_autovalores(matriz_J)
    #R=np.sqrt(N*((1-a)*(1/12)*1.**2+a*(1/12)*(1.*w)**2))
    #plot_autovalores(autovalores, R)
    
    for i in range(len(matriz_J)):
        I  = np.sum(matriz_J[i, :])/N
        for j in range(len(matriz_J)):
            matriz_W[i][j] = matriz_J[i][j] + matriz_M[i][j] - I
    
    return matriz_W



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

def generar_matriz_W2(filas, columnas, w):
    
    matriz_W = np.zeros((filas, columnas)) 
    
    # Llenar el 80% de las columnas con valores aleatorios entre 0 y 1
    for i in range(columnas):
        if i < 0.8 * columnas:
            matriz_W[:, i] = np.random.uniform(0., 1., filas)
        else:
            # Llenar el 20% restante con valores aleatorios entre 0 y w
            matriz_W[:, i] = np.random.uniform(-1.*w, 0., filas)
    
    condicion1(filas, columnas, matriz_W)
    
    return matriz_W


def generar_matriz_W1y2(filas, columnas, w):
    matriz_J = np.zeros((filas, columnas))  # Inicializar matriz con cero
    matriz_M = np.zeros((filas, columnas)) 
    matriz_W = np.zeros((filas, columnas)) 
    
    # Llenar el 80% de las columnas con valores aleatorios entre 0 y 1
    for i in range(columnas):
        if i < 0.8 * columnas:
            matriz_J[:, i] = np.random.uniform(-1./2., 1./2., filas)
            matriz_M[:, i] = 1./2.
        else:
            # Llenar el 20% restante con valores aleatorios entre 0 y w
            matriz_J[:, i] = np.random.uniform(-1.*w/2., 1.*w/2., filas)
            matriz_M[:, i] = -1.*w/2.

    # Aplicar la condición de sparsity sincronizada a ambas matrices
    #matriz_J = condicion1y2(filas, columnas, matriz_J, matriz_M)
    condicion1(filas, columnas, matriz_J)
    
    # Sumar las matrices J y M para obtener W
    for i in range(len(matriz_J)):
        I  = np.sum(matriz_J[i, :])/N
        for j in range(len(matriz_J[0])):
            matriz_W[i][j] = matriz_J[i][j] + matriz_M[i][j] - I  

    return matriz_W

#La siguiente funcion hace 0 el 80% de elementos de la matriz    
def condicion1(filas, columnas, matriz):            
    # Determinar el número de elementos a convertir en cero
    num_zeros = int(0.8 * filas * columnas)
    
    # Generar posiciones aleatorias para establecer en cero
    indices = np.random.choice(filas * columnas, num_zeros, replace=False)
    filas_indices, columnas_indices = np.unravel_index(indices, (filas, columnas))
    
    # Establecer en cero los elementos seleccionados
    for fila, columna in zip(filas_indices, columnas_indices):
        matriz[fila, columna] = 0
       
    return matriz



#Ahora, creamos una función que haga la suma de todas las filas, y modifique los elementos de cada fila para que la suma sea 0.
def condicion2(matriz):
    # Convertir la lista a un array de numpy si no lo es
    matriz = np.array(matriz)
    
    # Calcular la suma de cada fila
    suma_filas = np.sum(matriz, axis=1)
    
    # Calcular la media de cada fila
    media_filas = suma_filas / matriz.shape[1]
    
    # Restar la media de cada fila a sus elementos correspondientes
    matriz_resultante = matriz - media_filas[:, np.newaxis]
    
    return matriz_resultante



#Creamos una función para hacer un histograma de los valores de la matriz

def histograma(matriz,w, metodo):
    data = matriz.flatten()

    # Crear el KDE Plot
    sns.kdeplot(data, shade=True)

    # Mostrar el gráfico
    plt.title(f'Distribución de pesos de la matriz W, w= {w}, {metodo}')
    plt.xlabel('Valor')
    plt.ylabel('Densidad')
    plt.show()

#Calcula los autovalores de una matriz
def calcular_autovalores(matriz):
    autovalores = np.linalg.eigvals(matriz)
    return autovalores

#Representa los autovalores de una matriz y el circulo de radio R
def plot_autovalores(complex_values, R, w, metodo):
    # Extraer partes reales e imaginarias
    real_parts = [z.real for z in complex_values]
    imaginary_parts = [z.imag for z in complex_values]

    # Crear una figura y un eje
    fig, ax = plt.subplots()

    # Graficar los valores complejos en el plano complejo
    ax.scatter(real_parts, imaginary_parts, color='blue', label='Valores complejos')

    # Crear un círculo de radio R
    theta = np.linspace(0, 2 * np.pi, 100)
    x_circle = R * np.cos(theta)
    y_circle = R * np.sin(theta)
    ax.plot(x_circle, y_circle, color='red')

    # Configuración de los ejes
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    ax.set_aspect('equal')

    # Añadir leyenda
    ax.legend()

    # Añadir títulos
    plt.title(f'Distribución de autovalores de la función W, w= {w}, {metodo}')
    plt.xlabel('Parte Real')
    plt.ylabel('Parte Imaginaria')

    # Mostrar gráfico
    plt.show()

#La siguiente función sirve para escribir la la matriz en un txt
def guardar_matriz_en_txt(matriz, nombre_archivo):
    with open(nombre_archivo, 'w') as archivo:
        for fila in matriz:
            fila_str = ' '.join(map(str, fila))
            archivo.write(fila_str + '\n')



#Ponemos las condiciones iniciales
N=1000
filas = N
columnas = N
w=3
g=0.06
a=0.2
p=0.2
nombre_archivo1 = 'matriz.txt'
metodo= "spars y filas"


#Generamos una matriz W a la que le imponemos la condicion de Rajann y Abott a J siendo W=J+W
#matriz = generar_matriz_W(filas, columnas, w)
#matriz=generar_matriz_W2(filas, columnas, w)
matriz=generar_matriz_W2(filas, columnas, w)
#condicion1(filas, columnas, matriz)
histograma(matriz, w, metodo)

#Vamos a representar los autovalores de W, para ello, calculamos tambien el radio del circulo con la expresion
#R=np.sqrt(N*((1-a)*(1/3.-1/4.)*g**2+a*(1/3.-1/4.)*(g*w)**2))  este es para gaussiana
#l = g/2.*N*p*(1-a)-g*w/2.*N*p*a
#print(l)
#Sparsity
R=g*np.sqrt(N*((1-a)*(p)+a*(w**2)*(p/3.-(p**2)/4.)))

#Para RA mas Sparsity
#R=np.sqrt(N*p*((1-a)*g**2/12+a*(g*w)**2/12))

#Calculamos los autovalores de la matriz W
autovalores=calcular_autovalores(g*matriz)

#Representamos los autovalores y el circulo de radio R
plot_autovalores(autovalores, R, w, metodo)


#guardar_matriz_en_txt(matriz, nombre_archivo_normal)



#print(f"Matriz guardada en '{nombre_archivo1}'")