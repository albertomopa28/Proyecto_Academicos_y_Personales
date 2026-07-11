import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Nombre del archivo con datos
archivo = "Cohete_respecto_a_Marte.dat"

# Leer datos del archivo
def leer_datos(archivo):
    datos = []
    with open(archivo, "r") as f:
        for linea in f:
            if linea.strip():  # Si la línea no está vacía
                try:
                    datos.append([float(x) for x in linea.split(",")])  # Convertir a float
                except ValueError:
                    print(f"Error al procesar la línea: {linea.strip()}")
    return np.array(datos)

# Procesar datos
datos = leer_datos(archivo)

# Depuración: mostrar los primeros puntos leídos
print("Ejemplo de datos leídos:")
print(datos[:5])  # Mostrar los primeros 5 puntos

# Crear figura para la animación
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(-200000, 200000)
ax.set_ylim(-200000, 200000)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Trayectoria del cohete respecto a Marte")

# Inicializar trayectoria y punto móvil
trayectoria_line, = ax.plot([], [], 'b-', label="Trayectoria")
cuerpo_point, = ax.plot([], [], 'ro', label="Cohete")

# Leyenda
ax.legend()

# Variables para acumular trayectorias
trayectoria_x, trayectoria_y = [], []

# Función de inicialización
def init():
    trayectoria_line.set_data([], [])
    cuerpo_point.set_data([], [])
    return trayectoria_line, cuerpo_point

# Función de actualización
def actualizar(frame):
    x, y = datos[frame]
    
    # Depuración: imprimir posiciones actuales
    print(f"Frame {frame}: Cohete ({x}, {y})")

    # Acumular trayectorias
    trayectoria_x.append(x)
    trayectoria_y.append(y)

    # Actualizar datos
    trayectoria_line.set_data(trayectoria_x, trayectoria_y)
    cuerpo_point.set_data([x], [y])

    return trayectoria_line, cuerpo_point

# Crear animación
fps = 120  # Ajustar FPS para la animación
interval = 1000 / fps
anim = FuncAnimation(fig, actualizar, frames=len(datos), init_func=init, interval=interval, blit=True)

# Guardar el video directamente
anim.save("Trayectoria_cohete_respecto_a_Marte.mp4", writer="ffmpeg", fps=fps)

# Cerrar la figura de la animación
plt.close(fig)

print("Video guardado como 'Trayectoria_cohete_respecto_a_Marte.mp4'")

# Crear figura para la trayectoria completa
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(-200000, 200000)
ax.set_ylim(-200000, 200000)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Trayectoria del cohete respecto a Marte")
ax.grid(True)

# Procesar todos los frames para acumular trayectorias completas
trayectoria_x, trayectoria_y = datos[:, 0], datos[:, 1]

# Dibujar la trayectoria completa
ax.plot(trayectoria_x, trayectoria_y, 'b-', label="Trayectoria")
ax.plot(trayectoria_x[-1], trayectoria_y[-1], 'ro', label="Posición final")  # Punto final

# Añadir leyenda
ax.legend()

# Mostrar la figura
plt.show()
