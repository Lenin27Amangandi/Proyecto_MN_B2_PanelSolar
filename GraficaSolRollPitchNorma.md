
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from pysolar.solar import get_azimuth, get_altitude
from datetime import datetime, timedelta
from pytz import timezone
import tkinter as tk
from tkinter import ttk


# Función para obtener la posición solar (azimut y elevación)
def get_solar_position(date, latitude=-0.2105367, longitude=-78.491614):
    az = get_azimuth(latitude, longitude, date)
    el = get_altitude(latitude, longitude, date)
    return az, el


# Función para calcular el error entre el vector solar y la normal del panel
def calculate_error(s, n):
    return np.linalg.norm(s - n)**2


# Función para actualizar los ángulos de pitch y roll usando el gradiente descendente
def gradient_descent(s, n, learning_rate=0.1, iterations=100):
    pitch = 0.0  # Inicializar el pitch
    roll = 0.0   # Inicializar el roll

    for _ in range(iterations):
        # Calcular el error
        error = calculate_error(s, n)

        # Calcular gradientes con respecto a pitch y roll
        grad_pitch = 2 * np.dot(s - n, [
            -np.cos(pitch) * np.cos(roll),
            0,
            -np.sin(pitch) * np.cos(roll)
        ])
        grad_roll = 2 * np.dot(s - n, [
            np.sin(pitch) * np.sin(roll),
            np.cos(roll),
            -np.cos(pitch) * np.sin(roll)
        ])

        # Actualizar los valores de pitch y roll
        pitch -= learning_rate * grad_pitch
        roll -= learning_rate * grad_roll

    return pitch, roll


# Función para simular el movimiento del sol y el panel solar
def simulate_sun_motion(date, duration, ax, sun_point, sun_vector, roll_vector, pitch_vector, normal_vector, angle_text, time_text, pitch_roll_text, times):
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

        # Actualizar el vector amarillo (representa la trayectoria del sol)
        sun_vector.set_data([0, sun_x], [0, sun_y])
        sun_vector.set_3d_properties([0, sun_z])

        # Calcular la normal del panel
        s = np.array([sun_x, sun_y, sun_z])  # Vector solar
        n = np.array([0, 1, 0])  # Normal del panel (inicialmente en dirección Z)

        # Calcular los ángulos de pitch y roll que minimizan el error
        pitch, roll = gradient_descent(s, n)

        # Guardar los resultados
        current_time = times[i].strftime("%H:%M:%S")
        results.append([current_time, az, el, pitch, roll])

        # Actualizar la hora actual en el gráfico
        time_text.set_text(f"Hora: {current_time}")

        # Actualizar el ángulo de azimut y elevación
        angle_text.set_text(f"Azimut: {az:.2f}°\nElevación: {el:.2f}°")

        # Vectores de roll y pitch
        # El vector de roll (ajuste de la orientación sobre el eje horizontal)
        roll_vector.set_data([0, np.cos(roll)], [0, np.sin(roll)])
        roll_vector.set_3d_properties([0, 0])

        # El vector de pitch (ajuste de la inclinación sobre el eje vertical)
        pitch_vector.set_data([0, np.cos(pitch)], [0, np.sin(pitch)])
        pitch_vector.set_3d_properties([0, 0])

        # El vector de la normal del panel (perpendicular al sol), más corto que el vector amarillo
        # Hacer el vector verde más corto
        normal_vector.set_data([0, sun_x * 0.8], [0, sun_y * 0.8])  # Reducción de longitud al 80% del vector amarillo
        normal_vector.set_3d_properties([0, sun_z * 0.8])  # Reducir también la componente Z

        # Mostrar los valores de roll y pitch en texto
        pitch_roll_text.set_text(f"Pitch: {np.degrees(pitch):.2f}°\nRoll: {np.degrees(roll):.2f}°")

        return sun_point, sun_vector, roll_vector, pitch_vector, normal_vector, angle_text, time_text, pitch_roll_text

    # Función para guardar los resultados al finalizar
    def save_results():
        print(f"{'Hora':<12} {'Azimut (°)':<15} {'Elevación (°)':<15} {'Pitch (°)':<15} {'Roll (°)':<15}")
        for result in results:
            print(f"{result[0]:<12} {result[1]:<15.2f} {result[2]:<15.2f} {result[3]:<15.2f} {result[4]:<15.2f}")

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
    ax.set_title("Simulación del Movimiento del Sol y Panel Solar")

    # Crear el punto que representa al sol
    sun_point, = ax.plot([], [], [], 'o', color='yellow', markersize=10, label='Trayectoria del Sol')

    # Crear el vector amarillo para representar la trayectoria del sol
    sun_vector, = ax.plot([], [], [], color='yellow', linewidth=2)

    # Crear el vector azul para representar el roll
    roll_vector, = ax.plot([], [], [], color='blue', linewidth=2, label='Vector de Roll')

    # Crear el vector rojo para representar el pitch
    pitch_vector, = ax.plot([], [], [], color='red', linewidth=2, label='Vector de Pitch')

    # Crear el vector verde para representar la normal del panel
    normal_vector, = ax.plot([], [], [], color='green', linewidth=2, label='Normal al Panel')

    # Crear texto para mostrar los ángulos
    angle_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

    # Crear texto para mostrar la hora
    time_text = ax.text2D(0.05, 0.90, "", transform=ax.transAxes)

    # Crear texto para mostrar los valores de roll y pitch
    pitch_roll_text = ax.text2D(0.05, 0.85, "", transform=ax.transAxes)

    # Llamar a la función de simulación para actualizar el gráfico
    ani = FuncAnimation(fig, simulate_sun_motion(start_datetime, simulation_duration, ax, sun_point, sun_vector, roll_vector, pitch_vector, normal_vector, angle_text, time_text, pitch_roll_text, []),
                        frames=int(simulation_duration * 60 / 10),
                        interval=500, blit=False)

    # Añadir leyenda para explicar qué representa cada vector
    ax.legend(loc="upper left")

    # Guardar el GIF
    ani.save('solar_simulation.gif', writer='imagemagick')

    # Mostrar el gráfico en Tkinter
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=4)
    canvas.draw()


# Configuración de la interfaz gráfica con Tkinter
window = tk.Tk()
window.title("Simulación del Movimiento del Sol y Panel Solar")

# Obtener la fecha y hora actuales
now = datetime.now()
current_date = now.strftime("%Y-%m-%d")
current_time = now.strftime("%H:%M")

# Etiquetas y entradas para la fecha, hora y duración de la simulación
ttk.Label(window, text="Fecha (YYYY-MM-DD):").grid(row=0, column=0)
date_entry = ttk.Entry(window)
date_entry.insert(0, current_date)  # Prellenar con la fecha actual
date_entry.grid(row=0, column=1)

ttk.Label(window, text="Hora de inicio (HH:MM):").grid(row=1, column=0)
time_entry = ttk.Entry(window)
time_entry.insert(0, current_time)  # Prellenar con la hora actual
time_entry.grid(row=1, column=1)

ttk.Label(window, text="Duración simulación (horas):").grid(row=2, column=0)
duration_entry = ttk.Entry(window)
duration_entry.grid(row=2, column=1)

# Botón para ejecutar la simulación
simulate_button = ttk.Button(window, text="Simular Automáticamente", command=run_simulation)
simulate_button.grid(row=3, column=0, columnspan=2)

# Ejecutar la interfaz
window.mainloop()
```
