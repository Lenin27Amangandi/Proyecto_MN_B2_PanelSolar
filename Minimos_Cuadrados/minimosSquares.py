import sys
# (Opcional) Mantener emojis en Windows si tu consola soporta UTF-8.
# Si no, puedes comentar esto y usar el print ASCII más abajo.
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.animation import FuncAnimation
from datetime import datetime, timedelta
import tkinter as tk
from tkinter import ttk, messagebox
import csv
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 🔹 Posición solar real
from pysolar.solar import get_altitude, get_azimuth
from pytz import timezone

# ===============================
# PARÁMETROS DEL LUGAR (EPN por defecto, SOLO INFORMATIVOS)
# ===============================
LATITUDE = -0.2105367
LONGITUDE = -78.491614
TZ = "America/Guayaquil"

# ===============================
# MÍNIMOS CUADRADOS (Roll y Pitch)
# ===============================
LS_ITERS = 236
MU = 0.4

def least_squares_roll_pitch(az, el):
    # Vector solar en coordenadas ENU
    s = np.array([
        np.cos(el) * np.sin(az),
        np.cos(el) * np.cos(az),
        np.sin(el)
    ])

    # Inicialización de ángulos
    phi, psi = 0.0, 0.0

    for _ in range(LS_ITERS):
        # Normal del panel
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
# POSICIÓN SOLAR REAL
# ===============================
def getSolarPosition(latitude, longitude, date):
    az = get_azimuth(latitude, longitude, date)   # grados
    el = get_altitude(latitude, longitude, date)  # grados
    return az, el

# ===============================
# SIMULACIÓN TEMPORAL
# ===============================
def run_simulation(date, start_hour, duration_hours):
    tz = timezone(TZ)
    t0 = tz.localize(datetime(date.year, date.month, date.day, start_hour, 0))
    times = [t0 + timedelta(minutes=5*i) for i in range(int(duration_hours*12))]

    az_hist, el_hist, phi_hist, psi_hist, time_hist = [], [], [], [], []

    print("\nHora\t\tAzimut(°)\tElevación(°)\tRoll(°)\t\tPitch(°)")
    print("-"*80)

    for t in times:
        az_deg, el_deg = getSolarPosition(LATITUDE, LONGITUDE, t)
        if el_deg <= 0:
            # Sol bajo el horizonte: omitir para evitar vectores degenerados
            continue

        az = np.deg2rad(az_deg)
        el = np.deg2rad(el_deg)

        phi, psi = least_squares_roll_pitch(az, el)

        az_hist.append(az)
        el_hist.append(el)
        phi_hist.append(phi)
        psi_hist.append(psi)
        time_hist.append(t)

        print(f"{t.strftime('%H:%M')}\t\t{az_deg:6.2f}\t\t{el_deg:6.2f}\t\t{np.rad2deg(phi):6.2f}\t\t{np.rad2deg(psi):6.2f}")

    # Guardar CSV (solo valores ASCII; no problemas de codificación)
    with open("resultados_simulacion.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Hora", "Azimut (°)", "Elevación (°)", "Roll (°)", "Pitch (°)"])
        for i in range(len(time_hist)):
            writer.writerow([
                time_hist[i].strftime("%H:%M"),
                np.rad2deg(az_hist[i]),
                np.rad2deg(el_hist[i]),
                np.rad2deg(phi_hist[i]),
                np.rad2deg(psi_hist[i])
            ])

    # Versión segura sin emoji (recomendada si tu consola no es UTF-8):
    print("\n[OK] Archivo 'resultados_simulacion.csv' generado con éxito.")
    # Si tu consola soporta UTF-8 y dejaste el reconfigure(), puedes usar:
    # print("\n📁 Archivo 'resultados_simulacion.csv' generado con éxito.")

    return az_hist, el_hist, phi_hist, psi_hist, time_hist

# ===============================
# GRÁFICAS SEPARADAS
# ===============================
def plot_separate_graphs(az, el, phi, psi, times):
    t = range(len(times))
    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axs[0].plot(t, np.rad2deg(az), color='blue')
    axs[0].set_ylabel("Azimut (°)")
    axs[0].grid(True)

    axs[1].plot(t, np.rad2deg(el), color='orange')
    axs[1].set_ylabel("Elevación (°)")
    axs[1].grid(True)

    axs[2].plot(t, np.rad2deg(phi), color='green')
    axs[2].set_ylabel("Roll (°)")
    axs[2].grid(True)

    axs[3].plot(t, np.rad2deg(psi), color='red')
    axs[3].set_ylabel("Pitch (°)")
    axs[3].set_xlabel("Tiempo (pasos)")
    axs[3].grid(True)

    fig.suptitle("Evolución de Ángulos del Seguidor Solar", fontsize=14)
    plt.tight_layout()
    plt.show()

# ===============================
# ANIMACIÓN 3D CON RELOJ Y COORDENADAS (ENCABEZADO)
# ===============================
def animate_simulation(az, el, phi, psi, times):
    import numpy as np
    from matplotlib.animation import FuncAnimation
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([0, 1])
    ax.set_xlabel("Este (E)")
    ax.set_ylabel("Norte (N)")
    ax.set_zlabel("Arriba (U)")

    # --- Plano suelo (opcional) y ejes de referencia ---
    ground = Poly3DCollection(
        [[[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]],
        alpha=0.12, color='gray'
    )
    ax.add_collection3d(ground)
    #ax.quiver(0, 0, 0, 0, 0.8, 0, color='k', lw=2)   # Norte
    #ax.text(0, 0.85, 0, "Norte 0°", color='k')
    #ax.quiver(0, 0, 0, 0.8, 0, 0, color='k', lw=2)   # Este
    #ax.text(0.85, 0, 0, "Este 90°", color='k')
    #ax.quiver(0, 0, 0, 0, 0, 0.8, color='k', lw=2)   # Zenith
    #ax.text(0, 0, 0.85, "Zenith 90°", color='k')

    # Encabezados
    clock_text = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, fontsize=12)
    ax.text2D(0.02, 0.90,
              f"Ubicación: Lat {LATITUDE:.4f}, Lon {LONGITUDE:.4f}",
              transform=ax.transAxes, fontsize=10)

    # Objetos dinámicos (se recrean cada frame)
    sun_vec = None       # vector verde
    proj_vec = None      # vector azul
    sun_point = None     # punto Sol (amarillo)
    vert_line = None     # línea vertical (proy -> Sol)
    panel_surface = None # panel

    def update(i):
        nonlocal sun_vec, proj_vec, sun_point, vert_line, panel_surface

        # --- Borrar elementos previos para que no dejen rastro ---
        for artist in (sun_vec, proj_vec, sun_point, vert_line, panel_surface):
            if artist is not None:
                try:
                    artist.remove()
                except Exception:
                    # Algunos devuelven listas; intentar remover el primero
                    try:
                        artist[0].remove()
                    except Exception:
                        pass
        sun_vec = proj_vec = sun_point = vert_line = panel_surface = None

        # --- Vector Sol y proyección horizontal ---
        s = np.array([
            np.cos(el[i]) * np.sin(az[i]),   # E
            np.cos(el[i]) * np.cos(az[i]),   # N
            np.sin(el[i])                    # U
        ])
        p = np.array([s[0], s[1], 0.0])

        # --- Dibujar vectores principales ---
        sun_vec = ax.quiver(0, 0, 0, s[0], s[1], s[2],
                            color='green', linewidth=3, arrow_length_ratio=0.1)
        proj_vec = ax.quiver(0, 0, 0, p[0], p[1], 0.0,
                             color='tab:blue', linewidth=2, arrow_length_ratio=0.1)

        # --- Punto del Sol y línea vertical entrecortada ---
        sun_point = ax.scatter(s[0], s[1], s[2], color='gold', s=120, edgecolor='k', linewidth=0.5)
        vert_line, = ax.plot([p[0], s[0]], [p[1], s[1]], [0, s[2]],
                             color='gray', linestyle='--', linewidth=1.5)

        # --- Panel (superficie) usando la normal n (de tus mínimos cuadrados) ---
        n = np.array([
            -np.sin(phi[i]) * np.cos(psi[i]),
            np.sin(psi[i]),
            np.cos(phi[i]) * np.cos(psi[i])
        ])
        size = 0.6
        u = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(u, n)) > 0.9:
            u = np.array([0.0, 1.0, 0.0])
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
        panel_surface = ax.add_collection3d(
            Poly3DCollection([corners], alpha=0.45, color='dodgerblue', edgecolor='none')
        )

        # --- Títulos ---
        ax.set_title(
            f"Az: {np.rad2deg(az[i]):.1f}°, El: {np.rad2deg(el[i]):.1f}°  |  "
            f"Roll: {np.rad2deg(phi[i]):.1f}°, Pitch: {np.rad2deg(psi[i]):.1f}°"
        )
        clock_text.set_text(f"Hora: {times[i].strftime('%H:%M')}")

        return sun_vec, proj_vec, sun_point, vert_line, panel_surface

    _ = FuncAnimation(fig, update, frames=len(az), interval=200, blit=False)
    plt.tight_layout()
    plt.show()

# ===============================
# SIMULACIÓN MANUAL INTERACTIVA
# ===============================
def interactive_simulation_manual():
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 6))

    # --- 3D ---
    ax3d = fig.add_subplot(121, projection='3d')
    ax3d.set_xlim([-1, 1])
    ax3d.set_ylim([-1, 1])
    ax3d.set_zlim([0, 1])
    ax3d.set_xlabel("Este (E)")
    ax3d.set_ylabel("Norte (N)")
    ax3d.set_zlabel("Arriba (U)")
    ax3d.set_title("Simulación Manual (vectores + Sol + panel)")

    # Plano suelo + ejes (opcionales)
    ground = Poly3DCollection(
        [[[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]],
        alpha=0.12, color='gray'
    )
    ax3d.add_collection3d(ground)
    #ax3d.quiver(0, 0, 0, 0, 0.8, 0, color='k', lw=2); ax3d.text(0, 0.85, 0, "Norte 0°", color='k')
    #ax3d.quiver(0, 0, 0, 0.8, 0, 0, color='k', lw=2); ax3d.text(0.85, 0, 0, "Este 90°", color='k')
    #ax3d.quiver(0, 0, 0, 0, 0, 0.8, color='k', lw=2); ax3d.text(0, 0, 0.85, "Zenith 90°", color='k')

    # --- 2D ---
    ax2d = fig.add_subplot(122)
    ax2d.set_xlim([0, 360])
    ax2d.set_ylim([0, 90])
    ax2d.set_xlabel("Azimut (°)")
    ax2d.set_ylabel("Elevación (°)")
    ax2d.set_title("Coloca el Sol")
    ax2d.grid(True)
    sun_point_2d, = ax2d.plot([150], [25], 'o', color='gold', markersize=8)

    # Deslizadores
    ax_az = plt.axes([0.25, 0.08, 0.5, 0.03])
    ax_el = plt.axes([0.25, 0.03, 0.5, 0.03])
    slider_az = Slider(ax_az, "Azimut (°)", 0, 360, valinit=150)
    slider_el = Slider(ax_el, "Elevación (°)", 0, 90, valinit=25)

    # Elementos dinámicos
    sun_vec = None       # vector verde
    proj_vec = None      # vector azul
    sun_point_3d = None  # Sol
    vert_line = None     # línea entrecortada
    panel_surface = None # panel

    ax3d.text2D(0.02, 0.92,
                f"Ubicación: Lat {LATITUDE:.4f}, Lon {LONGITUDE:.4f}",
                transform=ax3d.transAxes, fontsize=10)

    def update(val):
        nonlocal sun_vec, proj_vec, sun_point_3d, vert_line, panel_surface

        # Borrar para no dejar rastro
        for artist in (sun_vec, proj_vec, sun_point_3d, vert_line, panel_surface):
            if artist is not None:
                try:
                    artist.remove()
                except Exception:
                    try:
                        artist[0].remove()
                    except Exception:
                        pass
        sun_vec = proj_vec = sun_point_3d = vert_line = panel_surface = None

        az = np.deg2rad(slider_az.val)
        el = np.deg2rad(slider_el.val)

        # Vectores
        s = np.array([
            np.cos(el) * np.sin(az),
            np.cos(el) * np.cos(az),
            np.sin(el)
        ])
        p = np.array([s[0], s[1], 0.0])

        # Sol + línea entrecortada
        sun_point_3d = ax3d.scatter(s[0], s[1], s[2], color='gold', s=120, edgecolor='k', linewidth=0.5)
        vert_line, = ax3d.plot([p[0], s[0]], [p[1], s[1]], [0, s[2]],
                               color='gray', linestyle='--', linewidth=1.5)

        # Vectores principales
        sun_vec = ax3d.quiver(0, 0, 0, s[0], s[1], s[2],
                              color='green', linewidth=3, arrow_length_ratio=0.1)
        proj_vec = ax3d.quiver(0, 0, 0, p[0], p[1], 0.0,
                               color='tab:blue', linewidth=2, arrow_length_ratio=0.1)

        # Panel (usando tus mínimos cuadrados)
        phi, psi = least_squares_roll_pitch(az, el)
        n = np.array([
            -np.sin(phi) * np.cos(psi),
            np.sin(psi),
            np.cos(phi) * np.cos(psi)
        ])
        size = 0.6
        u = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(u, n)) > 0.9:
            u = np.array([0.0, 1.0, 0.0])
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
            Poly3DCollection([corners], alpha=0.45, color='dodgerblue', edgecolor='none')
        )

        # UI 2D y título
        sun_point_2d.set_data([slider_az.val], [slider_el.val])
        ax2d.set_title(f"Coloca el Sol — Az: {slider_az.val:.1f}°, El: {slider_el.val:.1f}°")
        ax3d.set_title(f"Manual — Az: {slider_az.val:.1f}°, El: {slider_el.val:.1f}°")

        fig.canvas.draw_idle()

    slider_az.on_changed(update)
    slider_el.on_changed(update)
    update(None)  # primer render
    plt.tight_layout()
    plt.show()

# ===============================
# INTERFAZ GRÁFICA PRINCIPAL (GUI)
# ===============================
def launch_gui():
    root = tk.Tk()
    root.title("Simulador Seguidor Solar – Mínimos Cuadrados + Posición Real")

    ttk.Label(root, text="Ubicación fija:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
    ttk.Label(root, text=f"Lat {LATITUDE:.4f}, Lon {LONGITUDE:.4f}").grid(row=0, column=1, padx=5, pady=5, sticky="w")

    ttk.Label(root, text="Fecha (YYYY-MM-DD):").grid(row=1, column=0, padx=5, pady=5)
    date_entry = ttk.Entry(root)
    date_entry.insert(0, datetime.now().strftime("%Y-%m-%d"))
    date_entry.grid(row=1, column=1, padx=5, pady=5)

    ttk.Label(root, text="Hora inicio (0–23):").grid(row=2, column=0, padx=5, pady=5)
    start_hour_entry = ttk.Entry(root)
    start_hour_entry.insert(0, "6")
    start_hour_entry.grid(row=2, column=1, padx=5, pady=5)

    ttk.Label(root, text="Duración simulación (horas):").grid(row=3, column=0, padx=5, pady=5)
    duration_entry = ttk.Entry(root)
    duration_entry.insert(0, "10")
    duration_entry.grid(row=3, column=1, padx=5, pady=5)

    def run():
        try:
            date = datetime.strptime(date_entry.get(), "%Y-%m-%d")
            start_hour = int(start_hour_entry.get())
            duration = float(duration_entry.get())
        except Exception:
            messagebox.showerror("Error", "Ingrese valores válidos.")
            return

        az, el, phi, psi, times = run_simulation(date, start_hour, duration)
        animate_simulation(az, el, phi, psi, times)
        plot_separate_graphs(az, el, phi, psi, times)

    def manual():
        interactive_simulation_manual()

    ttk.Button(root, text="▶ Ejecutar Simulación", command=run).grid(row=4, column=0, columnspan=2, pady=10)
    ttk.Button(root, text="🕹 Simulación Manual (Interactiva)", command=manual).grid(row=5, column=0, columnspan=2, pady=10)

    root.mainloop()

# ===============================
# EJECUCIÓN
# ===============================
if __name__ == "__main__":
    launch_gui()