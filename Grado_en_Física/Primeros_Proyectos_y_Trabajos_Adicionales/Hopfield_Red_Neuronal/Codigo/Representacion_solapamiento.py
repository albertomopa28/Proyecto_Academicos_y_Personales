import numpy as np
import matplotlib.pyplot as plt

# Configuración de parámetros
T_inicial = 0.0
T_final = 0.2
DeltaT = 0.0005
mu = 3  # Número de patrones

# Simulación ficticia de datos para reemplazar la salida de f2 (puedes cargar datos reales desde un archivo CSV)
num_temperaturas = int((T_final - T_inicial) / DeltaT) + 1
temperaturas = np.linspace(T_inicial, T_final, num_temperaturas)
solapamientos = np.random.rand(num_temperaturas, mu)  # Generar solapamientos ficticios

# Visualización de los solapamientos
plt.figure(figsize=(10, 6))
for i in range(mu):
    plt.plot(temperaturas, solapamientos[:, i], label=f"Patrón {i+1}")

plt.xlabel("Temperatura (T)")
plt.ylabel("Solapamiento promedio")
plt.title("Evolución del solapamiento respecto a la temperatura")
plt.legend()
plt.grid()
plt.tight_layout()

# Guardar la figura
output_figure = "hopfield_solapamiento.png"
plt.savefig(output_figure)
plt.show()

print(f"Gráfica completada. Guardada como {output_figure}.")
