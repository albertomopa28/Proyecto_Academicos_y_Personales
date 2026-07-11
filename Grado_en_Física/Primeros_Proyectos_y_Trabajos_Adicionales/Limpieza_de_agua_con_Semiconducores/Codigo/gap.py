import matplotlib.pyplot as plt


E_f_dis = -2.78  # Energia de fermi de la disolución en eV

# Constantes
eV_to_J = 1.60218e-19  # Conversión de eV a julios

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
            print("Entrando en sección EV")
        elif 'EC' in linea:
            en_ev = False
            en_ec = True
            print("Entrando en sección EC")
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

def obtener_extremos(ev_data, ec_data):
    e_ev_values = [e for k, e in ev_data]
    e_ec_values = [e for k, e in ec_data]
    e_ev_max = max(e_ev_values) if e_ev_values else None
    e_ec_min = min(e_ec_values) if e_ec_values else None
    return e_ev_max, e_ec_min

def calcular_longitud_onda(gap):
    h = 4.135667696e-15  # eV·s (constante de Planck)
    c = 299792458  # m/s (velocidad de la luz)
    longitud_onda = h * c / gap * 1e9  # Convertir de metros a nanómetros
    return longitud_onda

def verificar_espectro_visible(longitud_onda):
    if 400 <= longitud_onda <= 700:
        return True
    return False

def graficar_bandas(ev_data, ec_data):
    if ev_data:
        k_ev, e_ev = zip(*ev_data)
    else:
        k_ev, e_ev = [], []
    if ec_data:
        k_ec, e_ec = zip(*ec_data)
    else:
        k_ec, e_ec = [], []

    plt.figure(figsize=(10, 6))
    plt.plot(k_ev, e_ev, label='EV', color='blue')
    plt.plot(k_ec, e_ec, label='EC', color='red')

    # Obtener valores extremos y dibujar líneas horizontales
    e_ev_max, e_ec_min = obtener_extremos(ev_data, ec_data)
    if e_ev_max is not None:
        plt.axhline(y=e_ev_max, color='blue', linestyle='--', label=f'EV Max: {e_ev_max:.2f} eV')
    if e_ec_min is not None:
        plt.axhline(y=e_ec_min, color='red', linestyle='--', label=f'EC Min: {e_ec_min:.2f} eV')
    
    # Rellenar el hueco entre las líneas horizontales
    if e_ev_max is not None and e_ec_min is not None:
        plt.fill_between([min(k_ev + k_ec), max(k_ev + k_ec)], e_ev_max, e_ec_min, color='yellow', alpha=0.3, label='Gap')
        
        # Calcular la diferencia y añadir texto al gráfico
        gap = abs(e_ev_max - e_ec_min)  # Valor absoluto del gap
        plt.text(max(k_ev + k_ec) * 0.5, (e_ev_max + e_ec_min) * 0.5, f'Gap: {gap:.2f} eV', 
                 horizontalalignment='center', verticalalignment='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

        # Calcular la longitud de onda
        longitud_onda = calcular_longitud_onda(gap)
        esta_en_visible = verificar_espectro_visible(longitud_onda)

        # Calcular valores adicionales
        E_F = (e_ev_max + e_ec_min) / 2
        E_B = E_F - E_f_dis
        gap_J = gap * eV_to_J
        E_F_J = E_F * eV_to_J
        E_B_J = E_B * eV_to_J

        # Imprimir en la terminal
        print(f"Gap energético: {gap:.2f} eV")
        print(f"Longitud de onda del gap: {longitud_onda:.2f} nm")
        if esta_en_visible:
            print("La longitud de onda está dentro del espectro visible.")
        else:
            print("La longitud de onda NO está dentro del espectro visible.")
            
        print(f"{gap} ")
        print(f"{E_F} ")

    plt.xlabel('k ((nm)^(-1)')
    plt.ylabel('E (eV)')
    plt.title('GaAs 2D')
    plt.legend()
    plt.grid(False)
    plt.show()

# Nombre del archivo de texto
nombre_archivo = 'Ge30.txt'
# Leer los datos del archivo
ev_data, ec_data = leer_datos_archivo(nombre_archivo)

# Graficar las bandas
graficar_bandas(ev_data, ec_data)
