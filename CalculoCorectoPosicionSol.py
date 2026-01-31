import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pysolar.solar import get_azimuth, get_altitude
from datetime import datetime, timedelta
from pytz import timezone
import tkinter as tk
from tkinter import ttk


# Coordenadas geográficas de la EPN (Escuela Politécnica Nacional)
latitude = -0.2105367
longitude = -78.491614

# Definir el intervalo de tiempo en minutos como constante
INTERVALO_MINUTOS = 5  # Cambia este valor cuando desees modificar el intervalo


# Función para obtener la posición solar (azimut y elevación)
def get_solar_position(date, latitude, longitude):
    # Asegurarse de que la fecha esté consciente de la zona horaria
    date = date.astimezone(timezone('America/Guayaquil'))
    az = get_azimuth(latitude, longitude, date)
    el = get_altitude(latitude, longitude, date)
    return az, el


# Función para ejecutar la simulación y generar las gráficas
def run_simulation():
    # Obtener los valores de fecha, hora y duración de la simulación
    date_str = date_entry.get()
    start_time_str = time_entry.get()
    simulation_duration = float(duration_entry.get())

    # Convertir fecha y hora a datetime
    try:
        start_datetime = datetime.strptime(f"{date_str} {start_time_str}", "%Y-%m-%d %H:%M")
    except ValueError:
        print("Formato de fecha o hora incorrecto.")
        return

    start_datetime = timezone('America/Guayaquil').localize(start_datetime)

    # Intervalo de tiempo para la simulación (usando la constante definida)
    interval = timedelta(minutes=INTERVALO_MINUTOS)
    times = [start_datetime + interval * i for i in range(int(simulation_duration * 60 / INTERVALO_MINUTOS))]

    # Obtener los ángulos de elevación y azimut para cada intervalo de tiempo
    az_angles = []
    el_angles = []
    for t in times:
        az, el = get_solar_position(t, latitude, longitude)
        az_angles.append(az)
        el_angles.append(el)

    # Crear las gráficas
    fig, axs = plt.subplots(2, 1, figsize=(10, 7))

    # Primer gráfico: Ángulo de elevación con respecto al tiempo
    axs[0].plot(times, el_angles, label='Elevación (θ)', color='blue')
    axs[0].set_xlabel('Tiempo')
    axs[0].set_ylabel('Elevación (°)')
    axs[0].set_title('Ángulo de Elevación vs Tiempo')
    axs[0].grid(True)

    # Segundo gráfico: Ángulo de azimut con respecto al tiempo
    axs[1].plot(times, az_angles, label='Azimut (α)', color='red')
    axs[1].set_xlabel('Tiempo')
    axs[1].set_ylabel('Azimut (°)')
    axs[1].set_title('Ángulo de Azimut vs Tiempo')
    axs[1].grid(True)

    # Ajustar el diseño de las gráficas
    plt.tight_layout()

    # Mostrar las gráficas en Tkinter
    canvas = FigureCanvasTkAgg(fig, master=window)  # Crear un canvas para las gráficas
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=4)  # Mostrar en la interfaz
    canvas.draw()


# Configuración de la interfaz gráfica con Tkinter
window = tk.Tk()
window.title("Simulación del Movimiento del Sol")

# Etiquetas y entradas para la fecha, hora y duración de la simulación
ttk.Label(window, text="Fecha (YYYY-MM-DD):").grid(row=0, column=0)
date_entry = ttk.Entry(window)
date_entry.grid(row=0, column=1)

ttk.Label(window, text="Hora de inicio (HH:MM):").grid(row=1, column=0)
time_entry = ttk.Entry(window)
time_entry.grid(row=1, column=1)

ttk.Label(window, text="Duración simulación (horas):").grid(row=2, column=0)
duration_entry = ttk.Entry(window)
duration_entry.grid(row=2, column=1)

# Botón para ejecutar la simulación
simulate_button = ttk.Button(window, text="Simulación Automática", command=run_simulation)
simulate_button.grid(row=3, column=0, columnspan=2)

# Ejecutar la interfaz
window.mainloop()
