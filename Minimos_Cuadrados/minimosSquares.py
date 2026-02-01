import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.animation import FuncAnimation
from datetime import datetime, timedelta
import tkinter as tk
from tkinter import ttk, messagebox
import csv
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 🔹 Librerías astronómicas
from pysolar.solar import get_altitude, get_azimuth
from pytz import timezone

# Configuración de salida
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

# PARÁMETROS FIJOS (EPN)
LATITUDE = -0.2105367
LONGITUDE = -78.491614
TZ = "America/Guayaquil"
MU = 0.4
TOLERANCIA = 1e-7  # Condición de parada

# ===============================
# LÓGICA MATEMÁTICA (MODIFICADA CON CONDICIÓN DE PARADA)
# ===============================
def least_squares_roll_pitch(az, el):
    s = np.array([np.cos(el)*np.sin(az), np.cos(el)*np.cos(az), np.sin(el)])
    phi, psi = 0.0, 0.0
    
    iters_realizadas = 0
    max_safe_iters = 1000 # Límite para evitar bucles infinitos
    
    for i in range(1, max_safe_iters + 1):
        n = np.array([-np.sin(phi)*np.cos(psi), np.sin(psi), np.cos(phi)*np.cos(psi)])
        e = n - s
        
        # Derivadas parciales (Gradiente)
        dJ_dphi = 2 * np.dot(e, [-np.cos(phi)*np.cos(psi), 0, -np.sin(phi)*np.cos(psi)])
        dJ_dpsi = 2 * np.dot(e, [np.sin(phi)*np.sin(psi), np.cos(psi), -np.cos(phi)*np.sin(psi)])
        
        phi_prev, psi_prev = phi, psi
        
        # Actualización
        phi -= MU * dJ_dphi
        psi -= MU * dJ_dpsi
        
        iters_realizadas = i
        
        # CONDICIÓN DE PARADA: ¿El cambio en ambos ángulos es menor a la tolerancia?
        if abs(phi - phi_prev) < TOLERANCIA and abs(psi - psi_prev) < TOLERANCIA:
            break
            
    return phi, psi, iters_realizadas

def getSolarPosition(latitude, longitude, date):
    return get_azimuth(latitude, longitude, date), get_altitude(latitude, longitude, date)

# ===============================
# RENDERIZADO 3D
# ===============================
def render_3d_scene(ax, az_rad, el_rad, phi_rad, psi_rad, time_str):
    ax.clear()
    ax.set_xlim([-1, 1]); ax.set_ylim([-1, 1]); ax.set_zlim([0, 1])
    ax.set_xlabel("Este (E)"); ax.set_ylabel("Norte (N)"); ax.set_zlabel("Arriba (U)")

    title_text = f"Az: {np.rad2deg(az_rad):.1f}°, El: {np.rad2deg(el_rad):.1f}° | Roll: {np.rad2deg(phi_rad):.1f}°, Pitch: {np.rad2deg(psi_rad):.1f}°"
    ax.text2D(0.5, 1.05, title_text, transform=ax.transAxes, fontsize=11, ha='center', weight='bold')
    ax.text2D(0.02, 0.95, f"Hora: {time_str}", transform=ax.transAxes, fontsize=10)

    ax.add_collection3d(Poly3DCollection([[[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]]], alpha=0.1, color='gray'))

    s = np.array([np.cos(el_rad)*np.sin(az_rad), np.cos(el_rad)*np.cos(az_rad), np.sin(el_rad)])
    p = np.array([s[0], s[1], 0.0])
    ax.quiver(0,0,0, s[0],s[1],s[2], color='green', linewidth=3, arrow_length_ratio=0.1)
    ax.quiver(0,0,0, p[0],p[1],0.0, color='tab:blue', linewidth=2, arrow_length_ratio=0.1)
    ax.scatter(s[0], s[1], s[2], color='gold', s=100, edgecolor='k')
    ax.plot([p[0], s[0]], [p[1], s[1]], [0, s[2]], color='gray', linestyle='--')

    n = np.array([-np.sin(phi_rad)*np.cos(psi_rad), np.sin(psi_rad), np.cos(phi_rad)*np.cos(psi_rad)])
    size = 0.6
    u = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(u, n)) > 0.9: u = np.array([0.0, 1.0, 0.0])
    v = np.cross(n, u); u = np.cross(v, n)
    u /= np.linalg.norm(u); v /= np.linalg.norm(v)
    corners = [-size*u - size*v, -size*u + size*v, size*u + size*v, size*u - size*v]
    ax.add_collection3d(Poly3DCollection([corners], alpha=0.4, facecolor='dodgerblue', edgecolor='royalblue'))

# ===============================
# GRÁFICAS (PANEL 2x2)
# ===============================
def plot_results(t_h, az_h, el_h, phi_h, psi_h):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Análisis de Ángulos - Seguidor Solar EPN', fontsize=14, weight='bold')
    time_labels = [t.strftime('%H:%M') for t in t_h]

    axs[0, 0].plot(time_labels, np.rad2deg(az_h), 'b-o', markersize=2); axs[0, 0].set_title('Azimut (°)')
    axs[0, 1].plot(time_labels, np.rad2deg(el_h), 'r-o', markersize=2); axs[0, 1].set_title('Elevación (°)')
    axs[1, 0].plot(time_labels, np.rad2deg(phi_h), 'g-o', markersize=2); axs[1, 0].set_title('Roll (°)')
    axs[1, 1].plot(time_labels, np.rad2deg(psi_h), 'm-o', markersize=2); axs[1, 1].set_title('Pitch (°)')

    for ax in axs.flat:
        ax.grid(True, alpha=0.3)
        ax.set_xticks(time_labels[::max(1, len(time_labels)//6)])
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# ===============================
# GUI Y CONTROLADORES
# ===============================
def launch_gui():
    root = tk.Tk()
    root.title("Control Solar EPN")
    frame = ttk.Frame(root, padding=20); frame.pack()

    ttk.Label(frame, text="Fecha (YYYY-MM-DD):").grid(row=0, column=0, pady=5)
    date_ent = ttk.Entry(frame); date_ent.insert(0, datetime.now().strftime("%Y-%m-%d")); date_ent.grid(row=0, column=1)

    ttk.Label(frame, text="Hora Inicio (0-23):").grid(row=1, column=0, pady=5)
    hour_ent = ttk.Entry(frame); hour_ent.insert(0, "7"); hour_ent.grid(row=1, column=1)

    ttk.Label(frame, text="Duración (Horas):").grid(row=2, column=0, pady=5)
    dur_ent = ttk.Entry(frame); dur_ent.insert(0, "10"); dur_ent.grid(row=2, column=1)

    def run_automatic():
        try:
            d = datetime.strptime(date_ent.get(), "%Y-%m-%d")
            h, dur = int(hour_ent.get()), float(dur_ent.get())
            tz = timezone(TZ)
            t0 = tz.localize(datetime(d.year, d.month, d.day, h, 0))

            az_h, el_h, phi_h, psi_h, t_h, csv_data = [], [], [], [], [], []
            print(f"\n{'HORA':<8} | {'AZIMUT':<8} | {'ELEV':<8} | {'ROLL':<8} | {'PITCH':<8} | {'ITERS'}")
            print("-" * 65)
            for i in range(int(dur*12)):
                t = t0 + timedelta(minutes=5*i)
                az_d, el_d = getSolarPosition(LATITUDE, LONGITUDE, t)
                if el_d > 0:
                    az_r, el_r = np.deg2rad(az_d), np.deg2rad(el_d)
                    p_r, s_r, iters = least_squares_roll_pitch(az_r, el_r)
                    
                    az_h.append(az_r); el_h.append(el_r); phi_h.append(p_r); psi_h.append(s_r); t_h.append(t)
                    
                    row = [t.strftime('%H:%M'), round(az_d,2), round(el_d,2), round(np.rad2deg(p_r),2), round(np.rad2deg(s_r),2), iters]
                    csv_data.append(row)
                    print(f"{row[0]:<8} | {row[1]:<8.2f} | {row[2]:<8.2f} | {row[3]:<8.2f} | {row[4]:<8.2f} | {row[5]}")
            
            with open('datos_solar.csv','w',newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Hora','Azimut','Elevacion','Roll','Pitch', 'Iters'])
                writer.writerows(csv_data)

            plot_results(t_h, az_h, el_h, phi_h, psi_h)

            fig = plt.figure(figsize=(12,7)); ax = fig.add_subplot(111, projection='3d')
            ani = FuncAnimation(fig, lambda i: render_3d_scene(ax, az_h[i], el_h[i], phi_h[i], psi_h[i], t_h[i].strftime('%H:%M')),
                                frames=len(az_h), interval=150, repeat=True)
            plt.show()
        except Exception as e: messagebox.showerror("Error", str(e))

    def run_manual():
        try:
            d = datetime.strptime(date_ent.get(), "%Y-%m-%d")
            tz = timezone(TZ)
            t_base = tz.localize(datetime(d.year, d.month, d.day, int(hour_ent.get()), 0))

            init_az, init_el, start_time_str = None, None, ""
            for minute in range(0, 1440, 15):
                t_check = t_base + timedelta(minutes=minute)
                az, el = getSolarPosition(LATITUDE, LONGITUDE, t_check)
                if el > 0:
                    init_az, init_el, start_time_str = az, el, t_check.strftime('%H:%M')
                    break
            
            if init_az is None:
                messagebox.showwarning("Aviso","No se encontró sol en la fecha seleccionada."); return

            fig = plt.figure(figsize=(12,7)); ax3d = fig.add_subplot(121, projection='3d'); ax2d = fig.add_subplot(122)
            ax2d.set_xlim([0,360]); ax2d.set_ylim([0,90]); ax2d.grid(True)
            sun_dot, = ax2d.plot([init_az],[init_el],'ro',markersize=8)

            fig.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.22, wspace=0.35)
            s_az = Slider(plt.axes([0.18,0.08,0.64,0.03]), "Azimut (°)", 0, 360, valinit=init_az)
            s_el = Slider(plt.axes([0.18,0.03,0.64,0.03]), "Elevación (°)", 0, 90, valinit=init_el)

            def update_m(val):
                az_r, el_r = np.deg2rad(s_az.val), np.deg2rad(s_el.val)
                p, s, _ = least_squares_roll_pitch(az_r, el_r)
                render_3d_scene(ax3d, az_r, el_r, p, s, f"Manual: {start_time_str}")
                sun_dot.set_data([s_az.val],[s_el.val])
                fig.canvas.draw_idle()

            s_az.on_changed(update_m); s_el.on_changed(update_m)
            update_m(None)
            plt.show()
        except Exception as e: messagebox.showerror("Error", str(e))

    ttk.Button(frame, text="▶ Simulación Automática (Tabla/CSV/Gráficas/GIF)", command=run_automatic).grid(row=4,columnspan=2,pady=10,sticky="ew")
    ttk.Button(frame, text="🕹 Control Manual (Inicia en Amanecer)", command=run_manual).grid(row=5,columnspan=2,pady=5,sticky="ew")
    root.mainloop()

if _name_ == "_main_":
    launch_gui()
