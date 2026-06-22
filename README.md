# Modelo SIR (Inmunidad Permanente)

Esta rama contiene el código principal del TFG dedicado al modelo clásico de Kermack y McKendrick (Modelo SIR). Aquí se modela la dinámica de propagación de una enfermedad donde los individuos adquieren inmunidad tras recuperarse.

Se incluyen simulaciones que abarcan desde la aproximación macroscópica (determinista) hasta las fluctuaciones microscópicas (estocásticas), así como el efecto de intervenciones sanitarias en la evolución de la epidemia.

---

## Descripción de los Archivos

A continuación, se detalla el propósito de cada uno de los *scripts* de R incluidos en esta carpeta:

* **`SIR_Determinista_Equilibrio_Endemico.R`**: Análisis del sistema determinista con dinámica vital (nacimientos y muertes). Muestra cómo el sistema oscila y converge hacia un punto de equilibrio endémico a largo plazo determinado por el teorema del umbral.
* **`Simulacion_SIR_Estocastica_Determinista.R`**: En este*script* se comparan las trayectorias generadas por el sistema de Ecuaciones Diferenciales Ordinarias (EDOs) frente a las trayectorias de salto discretas obtenidas mediante el algoritmo exacto de Gillespie.
* **`Aproximacion_Van_Kampen_SIR.R`**: Implementación de las Ecuaciones Diferenciales Estocásticas (en inglés, SDE). Utiliza la expansión del tamaño del sistema de Van Kampen para aislar y simular el ruido gaussiano alrededor de la trayectoria determinista.
* **`Visualizacion_3D_Distribucion_Cuasiestacionaria.R`**: Código centrado en la representación gráfica tridimensional de la distribución de probabilidad (cuasi-estacionaria) de los estados del sistema en su fase endémica.
* **`Caso_1Intervencion.R`**: Simulación de un escenario epidemiológico en el que se aplican múltiples medidas de mitigación en una única intervención (por ejemplo, un confinamiento o cuarentena) que alteran las tasas iniciales que hacen que la epidemía perdure, tratando de conseguir la mitigación total. 
* **`Caso_2Intervenciones.R`**: Extensión del escenario anterior aplicando múltiples intervenciones de control sanitario para observar la respuesta del sistema.
