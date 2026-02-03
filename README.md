# Proyecto_MN_B2_PanelSolar

## Descripción

Este proyecto implementa una simulación de un **seguidor solar de dos grados de libertad (2 GDL)** desarrollado en Python.  
El sistema orienta automáticamente un panel solar para que su superficie sea lo más perpendicular posible a la dirección de incidencia del Sol, maximizando la captación de energía.

La simulación integra el cálculo de la posición solar, un modelo matemático de orientación del panel, un método numérico de optimización por mínimos cuadrados y visualizaciones gráficas en 2D y 3D.

---

## Marco Teórico

### Posición solar

La posición del Sol se describe mediante dos ángulos principales:

- **Azimut (α):** ángulo horizontal medido desde el norte.
- **Elevación (θ):** ángulo vertical del Sol respecto al horizonte.

Estos valores dependen de la fecha, hora y ubicación geográfica (latitud y longitud) y se calculan mediante modelos astronómicos.

### Seguidor solar de dos grados de libertad

El seguidor solar se modela como un sistema mecánico con dos grados de libertad:

- **Roll (φ):** rotación alrededor del eje que apunta al norte.
- **Pitch (ψ):** rotación alrededor del eje que apunta al este.

No se considera el ángulo **yaw**, ya que el sistema sigue una configuración fija similar a la empleada en los seguidores solares de la EPN.

### Modelo matemático del panel

La orientación del panel se representa mediante un **vector normal** \( n \), definido en función de los ángulos de roll y pitch.  
La dirección del Sol se representa mediante un **vector unitario** \( s \), calculado a partir del azimut y la elevación.

El objetivo del sistema es alinear el vector normal del panel con el vector solar.

---

## Métodos Utilizados

### Método 1: Cálculo de la posición solar (Modelo físico–astronómico)

Este método permite determinar la posición del Sol a lo largo del día utilizando un modelo astronómico implementado en la librería **PySolar**.

**Entradas:**
- Latitud
- Longitud
- Fecha y hora

**Salidas:**
- Azimut (α)
- Elevación (θ)

Este método no es iterativo, pero es fundamental ya que proporciona los datos de entrada para el método numérico de optimización.

```python
from pysolar.solar import get_altitude, get_azimuth
from datetime import datetime
from pytz import timezone

def getSolarPosition(
    latitude: float = -0.2105367,
    longitude: float = -78.491614,
    date: datetime = datetime.now(tz=timezone("America/Guayaquil")),
):
    az = get_azimuth(latitude, longitude, date)
    el = get_altitude(latitude, longitude, date)
    return az, el

azimuth, elevation = getSolarPosition()
print(f"Azimuth: {azimuth:.2f}°, Elevation: {elevation:.2f}°")
```

---

### Método 2: Optimización por mínimos cuadrados (método numérico)

El segundo método corresponde al **método numérico principal del proyecto**.  
Su objetivo es calcular los ángulos de **roll** y **pitch** que alinean el panel con la dirección del Sol.

Se define la función de costo:

\[
J = \| n - s \|^2
\]

La minimización de esta función se realiza mediante un **proceso iterativo tipo descenso por gradiente**, ajustando progresivamente los ángulos hasta cumplir un criterio de convergencia.

```python
def least_squares_roll_pitch(az, el):
    s = np.array([np.cos(el)*np.sin(az),
                  np.cos(el)*np.cos(az),
                  np.sin(el)])

    phi, psi = 0.0, 0.0
    for _ in range(1000):
        n = np.array([-np.sin(phi)*np.cos(psi),
                       np.sin(psi),
                       np.cos(phi)*np.cos(psi)])

        error = n - s
        phi -= 0.4 * np.dot(error, [-np.cos(phi)*np.cos(psi), 0, -np.sin(phi)*np.cos(psi)])
        psi -= 0.4 * np.dot(error, [np.sin(phi)*np.sin(psi), np.cos(psi), -np.cos(phi)*np.sin(psi)])

    return phi, psi
```

---

## Metodología General

1. Ingreso de parámetros iniciales (fecha, hora y duración).
2. Cálculo de la posición solar mediante el **Método 1**.
3. Construcción del vector solar.
4. Optimización de los ángulos del panel mediante el **Método 2**.
5. Registro de resultados en archivos CSV.
6. Visualización de resultados en gráficas y animaciones.

---

## Resultados

La simulación permite observar:

- La trayectoria diaria del Sol.
- La evolución temporal de los ángulos de roll y pitch.
- La correcta orientación del panel durante el período simulado.

Se generan gráficos 2D, una escena 3D interactiva y una animación que ilustran el comportamiento del sistema.

---

## Conclusiones

El proyecto demuestra que un **seguidor solar de dos grados de libertad** es suficiente para mantener una orientación adecuada del panel durante el día.

El uso del **método de mínimos cuadrados** proporciona una solución numérica estable y eficiente para el cálculo de los ángulos de control, validando el enfoque aplicado en la simulación.

El cálculo de la posición solar mediante un **modelo astronómico** permite obtener valores precisos de azimut y elevación, garantizando que los datos de entrada del sistema representen adecuadamente el movimiento real del Sol y sirvan como base confiable para la simulación.

---

## Tecnologías Utilizadas

- Python 3
- NumPy
- Matplotlib
- Tkinter
- PySolar
