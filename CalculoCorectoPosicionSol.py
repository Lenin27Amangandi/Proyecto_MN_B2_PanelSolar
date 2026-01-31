from datetime import datetime, timedelta
from pysolar.solar import get_azimuth, get_altitude
from pytz import timezone

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

# Función principal
def run_simulation():
    # Obtener la fecha y hora de inicio
    fecha_str = input("Ingrese la fecha (YYYY-MM-DD): ")
    hora_str = input("Ingrese la hora de inicio (HH:MM): ")

    try:
        # Convertir la fecha y hora ingresadas en un objeto datetime
        start_datetime = datetime.strptime(f"{fecha_str} {hora_str}", "%Y-%m-%d %H:%M")
    except ValueError:
        print("Formato de fecha o hora incorrecto.")
        return

    # Hacer que la hora esté consciente de la zona horaria de Ecuador
    start_datetime = timezone('America/Guayaquil').localize(start_datetime)

    # Obtener el tiempo de simulación
    tiempo_simulacion = float(input("Ingrese el tiempo de simulación en horas: "))

    # Crear el intervalo de tiempo para la simulación (usando la constante definida)
    interval = timedelta(minutes=INTERVALO_MINUTOS)

    # Crear una lista de tiempos para la simulación (usando el intervalo de tiempo fijo)
    times = [start_datetime + interval * i for i in range(int(tiempo_simulacion * 60 / INTERVALO_MINUTOS))]

    # Imprimir los resultados
    print("\nPosición Solar para cada intervalo (Azimut, Elevación):")
    for t in times:
        az, el = get_solar_position(t, latitude, longitude)
        print(f"Hora: {t.strftime('%Y-%m-%d %H:%M:%S')} | Azimut: {az:.2f}° | Elevación: {el:.2f}°")

# Ejecutar la simulación
if __name__ == "__main__":
    run_simulation()
