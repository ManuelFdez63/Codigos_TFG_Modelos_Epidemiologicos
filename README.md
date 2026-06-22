# Códigos TFG: Modelos Epidemiológicos

Este repositorio contiene las implementaciones de las las simulaciones numéricas desarrolladas en R para el Trabajo de Fin de Grado (TFG) sobre modelización estocástica de la propagación epidemiológica. 

Para facilitar el desarrollo independiente y el orden de los scripts, los códigos de cada modelo se encuentra organizado en ramas independientes dentro de este repositorio.

---

## Estructuras de las Ramas
 
Para ver los códigos de un escenario concreto hay que seleccionar la rama correspondiente en el menú desplegable de GitHub:

* **Rama `Modelo_SIR`**: Contiene el código del modelo clásico con inmunidad permanente. Incluye las simulaciones deterministas (EDOs), estocásticas exactas (Algoritmo de Gillespie) y aproximadas (SDE).
* **Rama `Modelo_SIRS`**: Contiene el código adaptado para enfermedades donde la inmunidad se pierde con el tiempo, analizando la transición hacia equilibrios endémicos.
* **Rama `Modelo_SIS`**: Contiene los códigos correspondientes a los resultados mostrados en el apéndice del trabajo, donde los individuos no adquieren inmunidad y regresan directamente al estado susceptible al recuperarse.

---

*Nota: En el README de cada rama encontrarás las instrucciones específicas para ejecutar los archivos correspondientes.*
