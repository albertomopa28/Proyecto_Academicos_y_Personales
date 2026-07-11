import numpy as np
import matplotlib.pyplot as plt

# Nombre del archivo
import numpy as np

# Nombre del archivo
archivo = "Energia.dat"

# Leer datos del archivo
data = np.loadtxt(archivo)

# Extraer columnas
tiempo = data[:, 0]  # Primera columna (Tiempo)
E = data[:, 1]       # Segunda columna (Energía total)
T = data[:, 2]       # Tercera columna (Energía cinética)
V = data[:, 3]       # Cuarta columna (Energía potencial)

# Crear gráfica

# Leer datos del archivo
data = np.loadtxt(archivo)

# Extraer columnas
tiempo = data[:, 0]  # Primera columna (Tiempo)
E = data[:, 1]       # Segunda columna (Energía total)
T = data[:, 2]       # Tercera columna (Energía cinética)
V = data[:, 3]       # Cuarta columna (Energía potencial)

# Crear gráfica
plt.figure(figsize=(12, 6))

# Graficar las energías
plt.plot(tiempo, E, label="E (Energía total)", color="blue")
plt.plot(tiempo, T, label="T (Energía cinética)", color="orange")
plt.plot(tiempo, V, label="V (Energía potencial)", color="green")

# Etiquetas y leyenda
plt.xlabel("Tiempo")
plt.ylabel("Energía")
plt.legend()
plt.grid()

# Mostrar la gráfica
plt.tight_layout()
plt.show()
