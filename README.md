# Modelo SIRS (Inmunidad Temporal)

Esta rama contiene el código adaptado al modelo epidemiológico SIRS. A diferencia del modelo clásico, aquí se introduce el concepto de pérdida de inmunidad: los individuos recuperados regresan al compartimento de susceptibles tras un periodo de tiempo, lo que permite el estudio de reinfecciones y dinámicas endémicas a largo plazo.

---

## Descripción de los Archivos

A continuación, se detalla el propósito de cada uno de los *scripts* de R incluidos en esta rama:

* **`SIRS_Determinista_Libre_Enfermedad.R`**: Análisis del sistema determinista en un escenario donde $\mathcal{R}_0 < 1$. Muestra cómo, a pesar de la pérdida de inmunidad, la infección no logra sostenerse y la población converge al equilibrio donde la enfermedad se extingue.
* **`SIRS_Determinista_Equilibrio_Endémico.R`**: Análisis del sistema determinista cuando $\mathcal{R}_0 > 1$. Ilustra la convergencia hacia un atractor o punto de equilibrio endémico, demostrando que la enfermedad persiste en la población de forma indefinida.
* **`Simulacion_SIRS_Estocastica_Determinista.R`**: *Script* comparativo que superpone las trayectorias generadas por las Ecuaciones Diferenciales Ordinarias (EDOs) con las trayectorias probabilísticas obtenidas mediante el algoritmo estocástico de Gillespie.
* **`Caso_Intervencionista.R`**: Simulación de un escenario con la aplicación de medidas de salud pública destinadas mitigar la progresión de la epidemia. Se evalúa si las intervenciones son suficientes para llevar el sistema desde un estado endémico hacia la erradicación de la enfermedad.
