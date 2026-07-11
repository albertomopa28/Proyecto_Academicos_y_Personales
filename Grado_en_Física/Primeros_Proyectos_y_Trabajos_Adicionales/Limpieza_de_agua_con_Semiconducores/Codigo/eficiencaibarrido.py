import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from decimal import Decimal, getcontext  # para mayor precision

# Variables a alterar dependiendo del material
# W 3D es 100nm para 2D y 1D depende del grosor y es hasta donde llega la influencia del potencial
# E_B y E_gap
# Indice de refraccion n se busca en internet para cada material, no cambia mucho creo
# movilidad DE HUECOS tambien se busca en internet
# sum(indices_intervalo) Numero de puntos que se usan para el ajuste linea 108
# temperature = 5800  Temperatura del cuerpo negro en Kelvin (aproximadamente la del sol) Linea 245
# num_photons = 5   Número que quieres que se absorban  y luego reaacionen pensad que la eficiencia es muy pequeña
# Configurar la precisión deseada (por ejemplo, 28 dígitos decimales)

getcontext().prec = 28

# Constantes físicas
h = 6.62607015e-34  # Constante de Planck en J*s
Pi = 3.141592
c = 3.0e8  # Velocidad de la luz en m/s
k = 1.380649e-23  # Constante de Boltzmann en J/K
q = 1.6022e-19  # carga del electron

# Constantes propias del material
B = 0.2  # Probabilidad de recombinacion en 1 segundo de momento no la usamos pq no la hemos encontrado,
# seria cambiar la p_recombinacion por B*t (el tiempo de llegar a la superficie)
t_recombinacion = 100*10**(-9)  # tiempo de recombinacion antes de llegar a la superficie
t_recombinacionS = 100*10**(-9)  # tiempo de recombinacion antes de llegar a la superficie
t_reaccion = 2.5*10**(-16)  # tiempo de reaccion
E_B = 5.04*10**(-19)  # Energia de la barrera en J aunque sea negativa ponerla positiva
E_gap = 1.27*10**(-19)  # E gab en J
n = 4  # Indice de refraccion?
movlidad = 2*5*10**(-2)  # movilidad huecos m2 V-1s-
nombre_archivo = 'Ge30.txt'  # Archivo que se va a leer
W = 1*10**(-9)  # Es donde empieza el potencial de la barrera en m (dato sacado del paper 100 nm) lo he cambiado a 1000nm

def planck_distribution(wavelength, temperature):
    """
    Calcula la distribución de Planck para una longitud de onda dada y una temperatura.
    :param wavelength: Longitud de onda en metros.
    :param temperature: Temperatura en Kelvin.
    :return: Intensidad espectral (radiancia) en W/m^3.
    """
    exponent = (h * c) / (wavelength * k * temperature)
    # Usamos una manera segura de calcular la intensidad para evitar overflow en np.exp
    with np.errstate(over='ignore'):
        exp_term = np.exp(exponent) - 1.0
    exp_term = np.where(exp_term == 0, np.inf, exp_term)  # evitar división por cero
    intensity = (2.0 * h * c**2) / (wavelength**5) / exp_term
    
    return intensity

# Función para calcular el tiempo exacto
def tiempo_llegada_superficie(x0):
    # La formula sale de v=-movilidad_huecos*E, integras y sale eso
    return np.log(W/x0)/(2*E_B/q/(W**2)*movlidad)

def leer_datos_archivo(nombre_archivo):
    with open(nombre_archivo, 'r') as archivo:
        lineas = archivo.readlines()

    ev_data = []
    ec_data = []
    en_ev = False
    en_ec = False

    for linea in lineas:
        linea = linea.strip()
        if 'EV' in linea:
            en_ev = True
            en_ec = False
        elif 'EC' in linea:
            en_ev = False
            en_ec = True
        elif linea.startswith('k ((nm)-1), E (eV)'):
            continue
        elif '------------------------------------------------------------' in linea:
            continue
        elif linea:
            try:
                valores = linea.split(',')
                k = float(valores[0].strip())
                e = float(valores[1].strip())
                if en_ev:
                    ev_data.append((k, e))
                elif en_ec:
                    ec_data.append((k, e))
            except ValueError as ve:
                print(f"No se puede convertir la línea: {linea}, error: {ve}")
    
    return ev_data, ec_data

def parabola(x, a, b, c):
    return a * x**2 + b * x + c

def encontrar_intervalo(datos, tipo, margen=0.1):
    k, e = zip(*datos)
    k = np.array(k)
    e = np.array(e)
    
    if tipo == 'min':
        indice_extremo = np.argmin(e)
    else:
        indice_extremo = np.argmax(e)

    k_extremo = k[indice_extremo]
    indices_intervalo = (k >= k_extremo - margen) & (k <= k_extremo + margen)
    
    # Ajustar margen si no hay suficientes puntos
    while np.sum(indices_intervalo) < 6:  # Numero de pumtos que usamos para el ajuste
        margen += 0.05
        indices_intervalo = (k >= k_extremo - margen) & (k <= k_extremo + margen)
    
    return k[indices_intervalo], e[indices_intervalo]

def ajustar_parabola(x, y):
    if len(x) < 3:
        raise ValueError("No hay suficientes puntos de datos para realizar un ajuste parabólico.")
    popt, _ = curve_fit(parabola, x, y)
    return popt

def calcular_masas_efectivas(nombre_archivo):
    ev_data, ec_data = leer_datos_archivo(nombre_archivo)
    
    # Ajuste para los datos EV alrededor del máximo
    k_ev_intervalo, e_ev_intervalo = encontrar_intervalo(ev_data, 'max')
    if len(k_ev_intervalo) < 3:
        raise ValueError("No hay suficientes puntos de datos en el intervalo de EV para realizar un ajuste parabólico.")
    a_v, b_v, c_v = ajustar_parabola(k_ev_intervalo, e_ev_intervalo)
    
    # Ajuste para los datos EC alrededor del mínimo
    k_ec_intervalo, e_ec_intervalo = encontrar_intervalo(ec_data, 'min')
    if len(k_ec_intervalo) < 3:  # Numero puntos que se usan en el ajuste
        raise ValueError("No hay suficientes puntos de datos en el intervalo de EC para realizar un ajuste parabólico.")
    a_c, b_c, c_c = ajustar_parabola(k_ec_intervalo, e_ec_intervalo)
    
    # Conversión de unidades: de nm^(-2) * eV a m^(-2) * J
    eV_to_J = 1.60218e-19
    nm_to_m = 1e-9
    conversion_factor = eV_to_J**2*(nm_to_m**2)
    
    a_c_converted = a_c * conversion_factor
    a_v_converted = a_v * conversion_factor
    
    return (h/(2*Pi))**2/a_c_converted, -(h/(2*Pi))**2/a_v_converted

def probabilidad(t, tau):
    """Calcula la probabilidad P(t) = 1 - exp(-t/tau)."""
    if t == 0:
        return 0  # Para evitar división por cero
    return 1 - np.exp(-t / tau)

def bucle_temporal(tau_recom, tau_reaccion):
    t = tau_recom / 1000  # Tiempo inicial
    while True:
        # Generar un número aleatorio y comparar con la probabilidad de recombinación

        # Generar otro número aleatorio y comparar con la probabilidad de reacción
        random_number2 = np.random.rand()
        P_reaccion = probabilidad(t, tau_reaccion)
        if random_number2 < P_reaccion:
            return 1
        
        random_number1 = np.random.rand()
        P_recom = probabilidad(t, tau_recom)
        if random_number1 < P_recom:
            return 0
        
        # Aumentar el tiempo en tau_recom / 10
        t += tau_recom / 20

def generate_photons_with_probability_and_recombination(temperature, num_photons, wavelength_min=1e-9, wavelength_max=3e-6):
    """
    Genera fotones con una distribución de acuerdo con la ley de Planck y probabilidad adicional dependiente de la longitud de onda
    y una probabilidad de recombinación constante.
    :param temperature: Temperatura del cuerpo negro en Kelvin.
    :param num_photons: Número de fotones a generar.
    :param wavelength_min: Longitud de onda mínima en metros.
    :param wavelength_max: Longitud de onda máxima en metros.
    :param recombination_prob: Probabilidad de recombinación (constante).
    :return: Longitudes de onda de los fotones generados.
    """
    A = 0  # numero de no absorciones
    R = 0  # Numero de recombinciones
    NR = 0  # Numero de veces que no reacciones
    
    # Generar una muestra de longitudes de onda
    wavelengths = np.linspace(wavelength_min, wavelength_max, 10000)
    
    # Calcular la distribución de Planck para cada longitud de onda
    intensities = planck_distribution(wavelengths, temperature)
    
    # Normalizar las intensidades para obtener una distribución de probabilidad
    probabilities = intensities / np.sum(intensities)
    
    photon_wavelengths = []
    total_attempts = 0
    
    while len(photon_wavelengths) < num_photons:
        total_attempts += 1
        
        # Seleccionar un fotón de acuerdo a la distribución de probabilidad
        chosen_wavelength = np.random.choice(wavelengths, p=probabilities)
        
        # Calculamos el coeficiente de absorcion
        # print('\n Energia foton:',h*c/chosen_wavelength)
        if h*c/chosen_wavelength - E_gap < 0:
            # la raiz de la alpha es negativa
            generation_probability = 0
            A = A + 1
        else:            
            # Calculamos las masas efectivas y el parametro de absorcion A*. Las masas efectivas se calcuan realizado un ajuste de E(k) respecto a los maximos o minimos de EC o EV
            m_e, m_h = calcular_masas_efectivas(nombre_archivo)
            A_star = (q**2 * (2 * (m_h * m_e) / (m_h + m_e))**(3/2)) / (n * c * h**2 * m_e)
            alpha = A_star*np.sqrt(h*c/chosen_wavelength - E_gap)  # todo en julios 
            # Calcular la probabilidad de generación
            generation_probability = 1 - np.exp(-alpha*(W*10**2))   # La W la cambiamos a cm por la teoria
            
            # Generar un número aleatorio y comparar con la probabilidad de generación
            if np.random.rand() > generation_probability:
                A = A + 1
            else:
                # Generar un número x aleatorio entre 0 y 10
                x = np.random.uniform(0, W)  # Generar un número aleatorio x entre 0 y W
                # Realmente la probabilidad no es la misma para todo x, creo que en vd seria ir cogiendose x aleatorios segun la distribucion e^-alpha*x pero como x es del orden de 10^-7 o 10^-9 y alpha de 10^4 va a salir que para todo x la probabilidad es 1
    
                # Ahora calculamos el tiempo que tardaría en llegar hasta la superficie teniendo en cuenta su posicion inicial
                t_lleagada = tiempo_llegada_superficie(W-x)    
                recombination_prob = probabilidad(t_lleagada, t_recombinacion) 
                # print('t_recombinacion: ', t_recombinacion)
                # print('t_lleagada: ', t_lleagada)
                # print('t_recombinacion/t_lleagada: ', t_lleagada/t_recombinacion)
                # print('exp(-t_recombinacion/t_lleagada): ', np.exp(-t_lleagada/t_recombinacion))
                # print('recombination_prob: ',recombination_prob)                    
                if np.random.rand() < recombination_prob:
                    print('Se recombina antes de llegar')
                    R = R + 1
                else:
                    if bucle_temporal(t_recombinacionS, t_reaccion) == 0:
                        print('No reacciona')
                        NR = NR + 1
                    else:
                        photon_wavelengths.append(chosen_wavelength)
                        print('REACCIONA')
    
    # Convertir los valores a Decimal para calcular la eficiencia
    efficiency = Decimal(len(photon_wavelengths)) / Decimal(total_attempts)
    
    PA = A / total_attempts
    PR = Decimal(R) / Decimal(total_attempts)
    PNR = Decimal(NR) / Decimal(total_attempts)
    print(f"Eficiencia: {100*efficiency:.16f} %")
    print(f"Porcentaje no absorcion: {100*PA:.16f} %")
    print(f"Porcentaje recombinacion: {100*PR:.16f} %")
    print(f"Porcentaje de no reaccion: {100*PNR:.16f} %")
    return np.array(photon_wavelengths), efficiency, PA, PR, PNR

# Parámetros
temperature = 5800  # Temperatura del cuerpo negro en Kelvin (aproximadamente la del sol)
num_photons = 5  # Número que quieres que se absorban pensad que la eficiencia es muy pequeña
recombination_prob = B  # Probabilidad de recombinación

# Listas para almacenar los valores
t_recombinacion_values = []
efficiency_values = []
PA_values = []
PR_values = []
PNR_values = []

t_recombinacion = 10**(-19)  # tiempo de recombinacion antes de llegar a la superficie
t_reaccion = 2.5*10**(-19)  # tiempo de reaccion

while t_recombinacion < 10**-7:
    # Generar fotones con la probabilidad de generación adicional y recombinación
    photon_wavelengths, efficiency, PA, PR, PNR = generate_photons_with_probability_and_recombination(temperature, num_photons)
    
    # Almacenar los valores en las listas
    t_recombinacion_values.append(t_recombinacion)
    efficiency_values.append(efficiency)
    PA_values.append(PA)
    PR_values.append(PR)
    PNR_values.append(PNR)
    
    # Incrementar los tiempos de recombinación y reacción
    t_recombinacion *= 10  # tiempo de recombinacion antes de llegar a la superficie
    t_reaccion *= 10  # tiempo de reaccion
    num_photons += 10

# Convertir las listas a arrays para facilidad de procesamiento
t_recombinacion_values = np.array(t_recombinacion_values, dtype=float)
efficiency_values = np.array(efficiency_values, dtype=float) * 100
PA_values = np.array(PA_values, dtype=float) * 100
PR_values = np.array(PR_values, dtype=float) * 100
PNR_values = np.array(PNR_values, dtype=float) * 100

# Crear el gráfico
plt.figure(figsize=(10, 6))
plt.plot(t_recombinacion_values, efficiency_values, label='Eficiencia')
plt.plot(t_recombinacion_values, PA_values, label='Porcentaje no absorcion')
plt.plot(t_recombinacion_values, PR_values, label='Porcentaje recombinacion')
plt.plot(t_recombinacion_values, PNR_values, label='Porcentaje de no reaccion')

plt.xscale('log')
plt.xlabel('Tiempo de recombinacion (s)')
plt.ylabel('Porcentaje (%)')
plt.title('Porcentajes en función del tiempo de recombinacion')
plt.legend()
plt.grid(True)
plt.show()



 

 