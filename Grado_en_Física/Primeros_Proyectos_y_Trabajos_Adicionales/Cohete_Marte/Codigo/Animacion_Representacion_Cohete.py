import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Nombre del archivo con datos
archivo = "cohete_data.dat"

# Leer datos del archivo
def leer_datos(archivo):
    datos = []
    with open(archivo, "r") as f:
        frame = []
        for linea in f:
            if linea.strip():  # Si la línea no está vacía
                try:
                    frame.append([float(x) for x in linea.split(",")])  # Convertir a float
                except ValueError:
                    print(f"Error al procesar la línea: {linea.strip()}")
            else:  # Línea vacía, indica un nuevo frame
                if frame:
                    datos.append(np.array(frame))
                frame = []  # Reiniciar el frame
        if frame:  # Agregar el último frame si no está vacío
            datos.append(np.array(frame))
    return datos

# Procesar datos
datos = leer_datos(archivo)

# Depuración: mostrar los primeros frames leídos
print("Ejemplo de datos leídos:")
for i, frame in enumerate(datos[:3]):  # Mostrar solo los primeros 3 frames
    print(f"Frame {i}: {frame}")

# Crear figura para la animación
fig_anim, ax_anim = plt.subplots(figsize=(10, 6))
ax_anim.set_xlim(-200000, 200000)
ax_anim.set_ylim(-200000, 200000)
ax_anim.set_xlabel("X")
ax_anim.set_ylabel("Y")
ax_anim.set_title("Movimiento del cohete, Tierra, Marte y Sol")

# Dibujar el Sol (posición fija en el origen)
sol, = ax_anim.plot(0, 0, 'yo', markersize=10, label="Sol")

# Inicializar trayectorias
cohete_line, = ax_anim.plot([], [], 'b-', label="Cohete")
tierra_line, = ax_anim.plot([], [], 'orange', label="Tierra")
marte_line, = ax_anim.plot([], [], 'g-', label="Marte")

# Inicializar puntos
cohete_point, = ax_anim.plot([], [], 'bo')
tierra_point, = ax_anim.plot([], [], 'o', color="orange")
marte_point, = ax_anim.plot([], [], 'go')

# Leyenda
ax_anim.legend()

# Variables para acumular trayectorias
cohete_x, cohete_y = [], []
tierra_x, tierra_y = [], []
marte_x, marte_y = [], []

# Función de inicialización
def init():
    cohete_line.set_data([], [])
    tierra_line.set_data([], [])
    marte_line.set_data([], [])
    cohete_point.set_data([], [])
    tierra_point.set_data([], [])
    marte_point.set_data([], [])
    return cohete_line, tierra_line, marte_line, cohete_point, tierra_point, marte_point

# Función de actualización
def actualizar(frame):
    cohete = datos[frame][0]  # Primera línea del frame
    tierra = datos[frame][1]  # Segunda línea del frame
    marte = datos[frame][2]   # Tercera línea del frame

    # Depuración: imprimir posiciones actuales
    print(f"Frame {frame}: Cohete {cohete}, Tierra {tierra}, Marte {marte}")

    # Asegurar que los datos sean flotantes (o listas de flotantes)
    cohete_x_val, cohete_y_val = float(cohete[0]), float(cohete[1])
    tierra_x_val, tierra_y_val = float(tierra[0]), float(tierra[1])
    marte_x_val, marte_y_val = float(marte[0]), float(marte[1])

    # Acumular trayectorias
    cohete_x.append(cohete_x_val)
    cohete_y.append(cohete_y_val)
    tierra_x.append(tierra_x_val)
    tierra_y.append(tierra_y_val)
    marte_x.append(marte_x_val)
    marte_y.append(marte_y_val)

    # Actualizar datos
    cohete_line.set_data(cohete_x, cohete_y)
    tierra_line.set_data(tierra_x, tierra_y)
    marte_line.set_data(marte_x, marte_y)

    cohete_point.set_data([cohete_x_val], [cohete_y_val])
    tierra_point.set_data([tierra_x_val], [tierra_y_val])
    marte_point.set_data([marte_x_val], [marte_y_val])

    return cohete_line, tierra_line, marte_line, cohete_point, tierra_point, marte_point

# Crear animación
fps = 120  # Ajustar FPS para la animación
interval = 1000 / fps
anim = FuncAnimation(fig_anim, actualizar, frames=len(datos), init_func=init, interval=interval, blit=True)

# Guardar el video directamente
anim.save("Trayectorias_planetas_cohete.mp4", writer="ffmpeg", fps=fps)
plt.close(fig_anim)  # Cerrar la figura de animación para que no se muestre

print("Video guardado como 'planetas_cohete.mp4'")

# Crear figura para la trayectoria completa
fig_trayectoria, ax_trayectoria = plt.subplots(figsize=(10, 6))
ax_trayectoria.set_xlim(-200000, 200000)
ax_trayectoria.set_ylim(-200000, 200000)
ax_trayectoria.set_xlabel("X")
ax_trayectoria.set_ylabel("Y")
ax_trayectoria.set_title("Trayectoria completa de los astros")
ax_trayectoria.grid(True)

# Dibujar el Sol (posición fija en el origen)
ax_trayectoria.plot(0, 0, 'yo', markersize=10, label="Sol")

# Variables para acumular trayectorias
tray_cohete_x, tray_cohete_y = [], []
tray_tierra_x, tray_tierra_y = [], []
tray_marte_x, tray_marte_y = [], []

# Procesar todos los frames para acumular trayectorias
for frame in datos:
    cohete = frame[0]  # Primera línea del frame
    tierra = frame[1]  # Segunda línea del frame
    marte = frame[2]   # Tercera línea del frame

    tray_cohete_x.append(cohete[0])
    tray_cohete_y.append(cohete[1])
    tray_tierra_x.append(tierra[0])
    tray_tierra_y.append(tierra[1])
    tray_marte_x.append(marte[0])
    tray_marte_y.append(marte[1])

# Dibujar las trayectorias completas
ax_trayectoria.plot(tray_cohete_x, tray_cohete_y, 'b-', label="Cohete")
ax_trayectoria.plot(tray_tierra_x, tray_tierra_y, 'orange', label="Tierra")
ax_trayectoria.plot(tray_marte_x, tray_marte_y, 'g-', label="Marte")

# Añadir leyenda
ax_trayectoria.legend()

# Mostrar la figura
plt.show()
