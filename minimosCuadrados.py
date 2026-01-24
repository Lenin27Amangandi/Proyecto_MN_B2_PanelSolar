import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
from datetime import datetime, timedelta
import tkinter as tk
from tkinter import ttk, messagebox

# ===============================
# MÓDULO MATEMÁTICO (MÍNIMOS CUADRADOS)
# ===============================
LS_ITERS = 60
MU = 0.4

def least_squares_roll_pitch(az, el):
    s = np.array([
        np.cos(el) * np.sin(az),
        np.cos(el) * np.cos(az),
        np.sin(el)
    ])

    phi, psi = 0.0, 0.0

    for _ in range(LS_ITERS):
        n = np.array([
            -np.sin(phi) * np.cos(psi),
            np.sin(psi),
            np.cos(phi) * np.cos(psi)
        ])

        e = n - s

        dJ_dphi = 2 * np.dot(e, [
            -np.cos(phi) * np.cos(psi),
            0,
            -np.sin(phi) * np.cos(psi)
        ])

        dJ_dpsi = 2 * np.dot(e, [
            np.sin(phi) * np.sin(psi),
            np.cos(psi),
            -np.cos(phi) * np.sin(psi)
        ])

        phi -= MU * dJ_dphi
        psi -= MU * dJ_dpsi

    return phi, psi


# ===============================
# MODELO SOLAR SIMPLIFICADO
# ===============================
def solar_position(date):
    hour = date.hour + date.minute / 60
    az = np.deg2rad(15 * (hour - 12) + 180)
    el = np.deg2rad(max(0, 60 - abs(hour - 12) * 7))
    return az, el


# ===============================
# SIMULACIÓN COMPLETA
# ===============================
def run_simulation(date, duration_hours):
    t0 = datetime(date.year, date.month, date.day, 6, 0)
    times = [t0 + timedelta(minutes=5*i) for i in range(int(duration_hours*12))]

    az_hist, el_hist, phi_hist, psi_hist = [], [], [], []

    for t in times:
        az, el = solar_position(t)
        if el <= 0:
            continue
        phi, psi = least_squares_roll_pitch(az, el)

        az_hist.append(az)
        el_hist.append(el)
        phi_hist.append(phi)
        psi_hist.append(psi)

    return az_hist, el_hist, phi_hist, psi_hist


# ===============================
# GRÁFICAS TEMPORALES
# ===============================
def plot_time_graphs(az, el, phi, psi):
    t = np.arange(len(phi))

    plt.figure(figsize=(10,6))
    plt.plot(t, np.rad2deg(phi), label="Roll (°)")
    plt.plot(t, np.rad2deg(psi), label="Pitch (°)")
    plt.plot(t, np.rad2deg(az), label="Azimut (°)")
    plt.plot(t, np.rad2deg(el), label="Elevación (°)")
    plt.xlabel("Tiempo (pasos)")
    plt.ylabel("Ángulo (°)")
    plt.title("Ángulos del Sistema vs Tiempo")
    plt.legend()
    plt.grid()
    plt.show()


# ===============================
# SIMULACIÓN 3D INTERACTIVA
# ===============================
def interactive_simulation_manual():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(wspace=0.45)
    fig.subplots_adjust(top=0.82, bottom=0.25)

    # --- Gráfico 3D ---
    ax3d = fig.add_subplot(121, projection='3d')
    ax3d.set_xlim([-1,1])
    ax3d.set_ylim([-1,1])
    ax3d.set_zlim([0,1])
    ax3d.set_xlabel("Este (E)")
    ax3d.set_ylabel("Norte (N)")
    ax3d.set_zlabel("Arriba (U)")
    ax3d.set_title("Simulación Manual del Seguidor Solar")

    # --- Gráfico 2D Azimut vs Elevación ---
    ax2d = fig.add_subplot(122)
    ax2d.set_xlim([0,360])
    ax2d.set_ylim([0,90])
    ax2d.set_xlabel("Azimut (°)")
    ax2d.set_ylabel("Elevación (°)")
    ax2d.set_title("Coloca el Sol")
    ax2d.grid(True)

    sun_point_2d, = ax2d.plot([150], [25], 'o', color='gold', markersize=10)

    sun_vec = None
    panel_surface = None
    sun_point_3d = None

    # --- Sliders ---
    ax_az = plt.axes([0.25, 0.08, 0.5, 0.03])
    ax_el = plt.axes([0.25, 0.03, 0.5, 0.03])

    slider_az = Slider(ax_az, "Azimut (°)", 0, 360, valinit=150)
    slider_el = Slider(ax_el, "Elevación (°)", 0, 90, valinit=25)

    # --- Botón volver a posición inicial ---
    ax_reset = plt.axes([0.01, 0.05, 0.14, 0.04])
    btn_reset = Button(ax_reset, "Volver a posición inicial")

    def update(val):
        nonlocal sun_vec, panel_surface, sun_point_3d

        az = np.deg2rad(slider_az.val)
        el = np.deg2rad(slider_el.val)

        # --- Vector solar ---
        s = np.array([
            np.cos(el)*np.sin(az),
            np.cos(el)*np.cos(az),
            np.sin(el)
        ])

        # --- Mínimos cuadrados ---
        phi, psi = least_squares_roll_pitch(az, el)

        n = np.array([
            -np.sin(phi)*np.cos(psi),
            np.sin(psi),
            np.cos(phi)*np.cos(psi)
        ])

        # --- Limpiar ---
        if sun_vec:
            sun_vec.remove()
        if panel_surface:
            panel_surface.remove()
        if sun_point_3d:
            sun_point_3d.remove()

        # --- Flecha solar ---
        sun_vec = ax3d.quiver(0,0,0, s[0], s[1], s[2], color='green', linewidth=3)

        # --- Punto del sol ---
        sun_point_3d = ax3d.scatter(s[0], s[1], s[2], color='gold', s=120)

        # --- Plano del panel ---
        size = 0.6
        u = np.array([1,0,0])
        if abs(np.dot(u, n)) > 0.9:
            u = np.array([0,1,0])
        v = np.cross(n, u)
        u = np.cross(v, n)
        u /= np.linalg.norm(u)
        v /= np.linalg.norm(v)

        corners = [
            -size*u - size*v,
            -size*u + size*v,
            size*u + size*v,
            size*u - size*v
        ]

        panel_surface = ax3d.add_collection3d(
            Poly3DCollection([corners], alpha=0.5, color='dodgerblue')
        )

        ax3d.set_title(
            f"Manual — Az: {slider_az.val:.1f}°, El: {slider_el.val:.1f}°  |  "
            f"Roll: {np.rad2deg(phi):.1f}°, Pitch: {np.rad2deg(psi):.1f}°"
        )

        # --- Gráfico 2D ---
        sun_point_2d.set_data([slider_az.val], [slider_el.val])
        ax2d.set_title(f"Coloca el Sol — Az: {slider_az.val:.1f}°, El: {slider_el.val:.1f}°")

        fig.canvas.draw_idle()

    # --- Función del botón ---
    def reset_position(event):
        slider_az.set_val(150)
        slider_el.set_val(25)

    slider_az.on_changed(update)
    slider_el.on_changed(update)
    btn_reset.on_clicked(reset_position)

    update(None)
    plt.show()


# ===============================
# ANIMACIÓN Y GIF
# ===============================
def create_gif(az, el, phi, psi):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([0,1])
    ax.set_xlabel("Este")
    ax.set_ylabel("Norte")
    ax.set_zlabel("Arriba")

    sun_vec = None
    panel_vec = None

    def update(i):
        nonlocal sun_vec, panel_vec

        s = np.array([
            np.cos(el[i])*np.sin(az[i]),
            np.cos(el[i])*np.cos(az[i]),
            np.sin(el[i])
        ])

        n = np.array([
            -np.sin(phi[i])*np.cos(psi[i]),
            np.sin(psi[i]),
            np.cos(phi[i])*np.cos(psi[i])
        ])

        if sun_vec:
            sun_vec.remove()
        if panel_vec:
            panel_vec.remove()

        sun_vec = ax.quiver(0,0,0, s[0], s[1], s[2], color='gold', linewidth=3)
        panel_vec = ax.quiver(0,0,0, n[0], n[1], n[2], color='blue', linewidth=3)

        ax.set_title(f"Seguimiento Solar – Frame {i+1}")
        return sun_vec, panel_vec

    ani = FuncAnimation(fig, update, frames=len(az), interval=200)
    ani.save("simulacion_seguidor_solar.gif", writer="pillow", fps=10)
    plt.show()


# ===============================
# INTERFAZ GRÁFICA (GUI)
# ===============================
def launch_gui():
    root = tk.Tk()
    root.title("Simulador Seguidor Solar – Mínimos Cuadrados")

    ttk.Label(root, text="Fecha (YYYY-MM-DD):").grid(row=0, column=0, padx=5, pady=5)
    date_entry = ttk.Entry(root)
    date_entry.insert(0, datetime.now().strftime("%Y-%m-%d"))
    date_entry.grid(row=0, column=1, padx=5, pady=5)

    ttk.Label(root, text="Duración simulación (horas):").grid(row=1, column=0, padx=5, pady=5)
    duration_entry = ttk.Entry(root)
    duration_entry.insert(0, "10")
    duration_entry.grid(row=1, column=1, padx=5, pady=5)

    def run():
        try:
            date = datetime.strptime(date_entry.get(), "%Y-%m-%d")
            duration = float(duration_entry.get())
        except:
            messagebox.showerror("Error", "Ingrese valores válidos.")
            return

        az, el, phi, psi = run_simulation(date, duration)
        plot_time_graphs(az, el, phi, psi)
        create_gif(az, el, phi, psi)

    ttk.Button(root, text="Ejecutar Simulación", command=run)\
        .grid(row=2, column=0, columnspan=2, pady=10)

    ttk.Button(root, text="Simulación Manual (Interactiva)", command=interactive_simulation_manual)\
        .grid(row=3, column=0, columnspan=2, pady=10)

    root.mainloop()


# ===============================
# EJECUCIÓN
# ===============================
if __name__ == "__main__":
    launch_gui()
