import numpy as np
import pandas as pd
import pvlib
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Configuración Geográfica (EPN, Quito)
LAT, LON = -0.2105367, -78.491614
TZ = 'America/Guayaquil'
ALTITUDE = 2800 

def obtener_vector_solar(fecha_hora):
    """
    Calcula el vector solar cartesiano en el marco ENU usando pvlib.
    """
    location = pvlib.location.Location(LAT, LON, tz=TZ, altitude=ALTITUDE)
    solpos = location.get_solarposition(fecha_hora)
    
    # Conversión a radianes
    az_rad = np.radians(solpos['azimuth'].values)
    zen_rad = np.radians(solpos['zenith'].values)
    
    # Conversión a Vector Cartesiano ENU (x=Este, y=Norte, z=Arriba)
    # Nota: Usamos .item() para asegurar que sean escalares y no arrays de 1 elemento
    Sx = np.sin(az_rad) * np.sin(zen_rad)
    Sy = np.cos(az_rad) * np.sin(zen_rad)
    Sz = np.cos(zen_rad)
    
    # CORRECCIÓN 1: Devolver el array con los valores calculados
    return np.array([Sx.item(), Sy.item(), Sz.item()])

def calculo_analitico_pitch_roll(vector_solar):
    """
    Calcula ángulos pitch y roll usando álgebra vectorial exacta.
    """
    Sx, Sy, Sz = vector_solar
    
    # Pitch (rotación eje X)
    theta_rad = np.arcsin(-Sy)
    
    # Roll (rotación eje Y)
    phi_rad = np.arctan2(Sx, Sz)
    
    return np.degrees(theta_rad), np.degrees(phi_rad)

def funcion_costo(angulos, vector_solar_objetivo):
    """Función a minimizar: 1 - coseno del ángulo entre vectores"""
    theta, phi = angulos # pitch, roll en radianes
    
    # CORRECCIÓN 2: Completar Matriz Pitch (X). Fila 1 es [1, 0, 0]
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(theta), -np.sin(theta)],
                   [0, np.sin(theta), np.cos(theta)]])
                   
    # CORRECCIÓN 3: Completar Matriz Roll (Y). Fila 2 es [0, 1, 0]
    Ry = np.array([[np.cos(phi), 0, np.sin(phi)],
                   [0, 1, 0],
                   [-np.sin(phi), 0, np.cos(phi)]])
    
    # CORRECCIÓN 4: Vector Normal Inicial (apunta arriba Z+)
    N0 = np.array([0, 0, 1])
    
    # Transformación: N = Ry * Rx * N0
    # Nota: El orden de multiplicación importa. Aquí asumimos rotación intrínseca o extrínseca
    # según tu diseño mecánico. Ry @ (Rx @ N0) rota primero en Pitch, luego en Roll.
    N_actual = Ry @ (Rx @ N0)
    
    # Producto punto (maximizar dot -> minimizar negativo)
    return -np.dot(N_actual, vector_solar_objetivo)

def calculo_numerico_pitch_roll(vector_solar):
    """Encuentra ángulos óptimos usando L-BFGS-B"""
    # CORRECCIÓN 5: Semilla inicial [pitch=0, roll=0]
    x0 = [0, 0]
    limites = [(-np.pi/2, np.pi/2), (-np.pi/2, np.pi/2)] 
    
    res = minimize(funcion_costo, x0, args=(vector_solar,), 
                   method='L-BFGS-B', bounds=limites, tol=1e-6)
    
    return np.degrees(res.x)

# Ejemplo de Uso
fechas = pd.date_range(start='2026-03-21 06:00', end='2026-03-21 18:00', freq='10min', tz=TZ)

# CORRECCIÓN 6: Inicializar lista vacía
resultados = []

for fecha in fechas:
    S = obtener_vector_solar(fecha)
    # Filtrar noche (Componente Z < 0 significa bajo el horizonte)
    if S[2] > 0: 
        pitch_ana, roll_anal = calculo_analitico_pitch_roll(S)
        pitch_num, roll_num = calculo_numerico_pitch_roll(S)
        resultados.append([fecha, pitch_ana, roll_anal, pitch_num, roll_num])

# CORRECCIÓN 7: Definir nombres de columnas
cols = ['Fecha', 'Pitch_Analitico', 'Roll_Analitico', 'Pitch_Numerico', 'Roll_Numerico']
df_res = pd.DataFrame(resultados, columns=cols)

print(df_res.head())

# Opcional: Graficar para verificar
plt.figure(figsize=(10, 5))
plt.plot(df_res['Fecha'], df_res['Pitch_Analitico'], label='Pitch (Analítico)', linestyle='--')
plt.plot(df_res['Fecha'], df_res['Pitch_Numerico'], label='Pitch (Numérico)', alpha=0.7)
plt.plot(df_res['Fecha'], df_res['Roll_Analitico'], label='Roll (Analítico)', linestyle='--')
plt.plot(df_res['Fecha'], df_res['Roll_Numerico'], label='Roll (Numérico)', alpha=0.7)
plt.legend()
plt.title('Comparación Cálculo Analítico vs Numérico - Quito')
plt.ylabel('Grados')
plt.xlabel('Hora')
plt.grid(True)
plt.show()