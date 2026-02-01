import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from pysolar.solar import get_azimuth, get_altitude
from datetime import datetime, timedelta
from pytz import timezone
import tkinter as tk
from tkinter import ttk
import matplotlib.animation as animation

def get_solar_position(date, latitude=-0.2105367, longitude=-78.491614):
    az = get_azimuth(latitude, longitude, date)
    el = get_altitude(latitude, longitude, date)
    return az, el


def simulate_sun_motion(date, duration, ax, sun_point, angle_text, sun_vector, time_text, times):
    interval = timedelta(minutes=10)
    times = [date + interval * i for i in range(int(duration * 60 / 10))]

    results = []

    def update_frame(i):
        az, el = get_solar_position(times[i])

        # Convertir los ángulos a radianes para la simulación 3D
        az_rad = np.deg2rad(az)
        el_rad = np.deg2rad(el)

        # Calcular la posición del sol en 3D
        sun_x = np.cos(el_rad) * np.sin(az_rad)
        sun_y = np.cos(el_rad) * np.cos(az_rad)
        sun_z = np.sin(el_rad)

        # Actualizar la posición del sol en la visualización
        sun_point.set_data([sun_x], [sun_y])
        sun_point.set_3d_properties([sun_z])

        # Actualizar el vector amarillo
        sun_vector.set_data([0, sun_x], [0, sun_y])
        sun_vector.set_3d_properties([0, sun_z])

        # Actualizar los ángulos de elevación y azimut en el gráfico
        angle_text.set_text(f"Azimut: {az:.2f}°\nElevación: {el:.2f}°")

        # Almacenar los resultados en la lista
        current_time = times[i].strftime("%H:%M:%S")
        results.append([current_time, az, el])

        # Actualizar la hora actual en el gráfico
        time_text.set_text(f"Hora: {current_time}")

        return sun_point, sun_vector, angle_text, time_text

    # Función para guardar los resultados al finalizar
    def save_results():
        # Imprimir los resultados en la shell solo una vez al final
        print(f"{'Hora':<12} {'Azimut (°)':<15} {'Elevación (°)':<15}")
        for result in results:
            print(f"{result[0]:<12} {result[1]:<15.2f} {result[2]:<15.2f}")

    # Guardar los resultados al terminar
    save_results()

    return update_frame


# Función para ejecutar la simulación al hacer clic en el botón
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

    # Crear la figura de la simulación 3D
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Limites y etiquetas de los ejes
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([0, 1])
    ax.set_xlabel("Este (E)")
    ax.set_ylabel("Norte (N)")
    ax.set_zlabel("Arriba (U)")
    ax.set_title("Simulación del Movimiento del Sol")

    # Crear el punto que representa al sol
    sun_point, = ax.plot([], [], [], 'o', color='yellow', markersize=10)

    # Crear un vector amarillo que seguirá la trayectoria del sol
    sun_vector, = ax.plot([], [], [], color='yellow', linewidth=2)

    # Crear texto para mostrar los ángulos
    angle_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

    # Crear texto para mostrar la hora
    time_text = ax.text2D(0.05, 0.90, "", transform=ax.transAxes)

    # Llamar a la función de simulación para actualizar el gráfico
    ani = FuncAnimation(fig, simulate_sun_motion(start_datetime, simulation_duration, ax, sun_point, angle_text, sun_vector, time_text, []),
                        frames=int(simulation_duration * 60 / 10),
                        interval=500, blit=False)

    # Guardar el GIF
    ani.save('solar_simulation.gif', writer='imagemagick')

    # Mostrar el gráfico en Tkinter
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=4)
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
simulate_button = ttk.Button(window, text="Simular Automáticamente", command=run_simulation)
simulate_button.grid(row=3, column=0, columnspan=2)

# Ejecutar la interfaz
window.mainloop()
