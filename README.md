# Modelo SIS (Sin Inmunidad - Apéndice)

Esta rama contiene el código complementario analizado en el apéndice del TFG. Corresponde al modelo epidemiológico SIS, diseñado para enfermedades infecciosas donde no existe periodo de inmunidad: los individuos infectados se recuperan e inmediatamente vuelven a ser susceptibles.

---

## Descripción del Archivo

* **`CodigoModeloSIS.R`**: *Script* unificado que contiene todo lo necesario para la visualización para este modelo. Al ejecutarse, este archivo se encarga de generar dos gráficas:
  1. **Análisis del sistema determinista con dinámica vital (nacimientos y muertes)**: Muestra cómo las curvas poblacionales de susceptibles e infectados oscilan y convergen hacia un punto de equilibrio endémico a largo plazo determinado por el teorema del umbral.
  2. **Comparativa estocástica frente a determinista**: Superpone las trayectorias generadas por el sistema de Ecuaciones Diferenciales Ordinarias (EDOs) frente a las trayectorias de salto discretas obtenidas mediante el algoritmo exacto de Gillespie, centrando la visualización únicamente en la evolución de la curva de los infectados para analizar el impacto del ruido.
