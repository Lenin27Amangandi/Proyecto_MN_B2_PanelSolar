# Métodos Mínimos Cuadrados y L-BFGS-B

## Instalación de las librerías necesarias

A continuación se describen las librerías utilizadas para el correcto funcionamiento de los métodos de **Mínimos Cuadrados** y **L-BFGS-B**, junto con su propósito y forma de instalación.

---

### `import sys`

Proporciona acceso a funciones y variables del intérprete de Python.

* 📦 **Instalación**: No requiere instalación adicional (incluida por defecto en Python).

---

### `import numpy as np`

Permite realizar cálculos numéricos eficientes y trabajar con matrices y arreglos multidimensionales.

* 📦 **Instalación**:

```bash
pip install numpy
```

---

### `import matplotlib.pyplot as plt`

Importa el submódulo `pyplot` de Matplotlib (alias `plt`). Se utiliza para crear gráficos 2D, histogramas y diagramas de dispersión.

* 📦 **Instalación**:

```bash
pip install matplotlib
```

---

### `from matplotlib.widgets import Slider`

Permite crear gráficos interactivos mediante barras deslizables para modificar parámetros en tiempo real.

* 📦 **Instalación**: Incluido en Matplotlib.

---

### `from matplotlib.animation import FuncAnimation`

Herramienta para crear animaciones actualizando repetidamente una función gráfica.

* 📦 **Instalación**: Incluido en Matplotlib.

---

### `from datetime import datetime, timedelta`

Importa clases para manipular fechas, horas y lapsos de tiempo.

* 📦 **Instalación**: Biblioteca estándar de Python.

---

### `import tkinter as tk`

Biblioteca estándar de Python para crear interfaces gráficas de usuario (GUI).

* 📦 **Instalación**: Incluida por defecto en la mayoría de instalaciones de Python.

---

### `from tkinter import ttk, messagebox`

Módulos de Tkinter para widgets avanzados y cuadros de diálogo.

* 📦 **Instalación**: Incluidos por defecto con Python.

---

### `import csv`

Permite leer y escribir archivos CSV (Comma Separated Values).

* 📦 **Instalación**: Biblioteca estándar de Python.

---

### `from mpl_toolkits.mplot3d.art3d import Poly3DCollection`

Se utiliza para dibujar polígonos 3D (caras y superficies) en gráficos tridimensionales.

* 📦 **Instalación**: Incluido en Matplotlib.

---

### `from pysolar.solar import get_altitude, get_azimuth`

Permite calcular la posición del sol (altitud y azimut) para una ubicación y tiempo determinados.

* 📦 **Instalación**:

```bash
pip install pysolar
```

---

### `from pytz import timezone`

Gestiona y convierte zonas horarias utilizando la base de datos Olson.

* 📦 **Instalación**:

```bash
pip install pytz
```

---

### `from scipy.optimize import minimize`

Función de SciPy para encontrar el mínimo de una función escalar, con o sin restricciones (incluye L-BFGS-B).

* 📦 **Instalación**:

```bash
pip install scipy
```

---

### `import plotly.graph_objects as go`

Interfaz de bajo nivel de Plotly para crear visualizaciones interactivas.

* 📦 **Instalación**:

```bash
pip install plotly
```

---

### `import pandas as pd`

Herramienta principal para manipulación, limpieza y análisis de datos estructurados mediante DataFrames.

* 📦 **Instalación**:

```bash
pip install pandas
```

---

### `from mpl_toolkits.mplot3d import Axes3D`

Habilita la creación de gráficos tridimensionales (3D) en Matplotlib.

* 📦 **Instalación**: Incluido por defecto en Matplotlib.

---

### `import matplotlib.animation as animation`

Se utiliza para crear visualizaciones animadas y dinámicas.

* 📦 **Instalación**:

```bash
pip install matplotlib
```

---

### `from IPython.display import HTML, display`

Permite mostrar contenido enriquecido (HTML, animaciones, widgets) en Jupyter Notebook o Google Colab.

* 📦 **Instalación**: Incluido en IPython (instalado por defecto con Jupyter / Anaconda).

---

📌 **Nota**: Para evitar conflictos, se recomienda usar un entorno virtual (`venv` o `conda`) antes de instalar las librerías.
