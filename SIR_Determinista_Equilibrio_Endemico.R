# --- 1. CONFIGURACIÓN DEL ENTORNO Y DEFINICIÓN DE PARÁMETROS ---
library(deSolve)      # Resolución numérica de sistemas de ecuaciones diferenciales ordinarias
library(ggplot2)      # Generación de representaciones gráficas de alta calidad
library(reshape2)     # Manipulación y estructuración de conjuntos de datos
library(gridExtra)    # Gestión de composiciones gráficas en paneles

# Definición de los coeficientes y constantes del sistema
params <- c(
  beta = 0.8,         # Coeficiente de transmisión por contacto efectivo
  gamma = 0.5,        # Tasa de recuperación 
  b = 0.05,           # Tasa de dinámica vital 
  N = 100             # Tamaño de la población total fija
)

# Especificación de las condiciones iniciales y del dominio temporal
initial_state <- c(S = 90, I = 10, R = 0)
times <- seq(0, 160, by = 0.1) # Horizonte temporal para la integración numérica

# --- 2. DETERMINACIÓN ANALÍTICA DE LOS PUNTOS DE EQUILIBRIO ---
# El comportamiento asintótico del sistema está determinado por el número básico de reproducción
R0 <- params["beta"] / (params["b"] + params["gamma"])

# Cálculo de las coordenadas del punto de equilibrio endémico (S*, I*, R*)
S_star <- params["N"] / R0
factor <- (1 - (1/R0)) / (params["b"] + params["gamma"])
I_star <- params["b"] * params["N"] * factor
R_star <- params["gamma"] * params["N"] * factor

# --- 3. MODELADO DETERMINISTA MEDIANTE EDOs ---
# Formulación del sistema de ecuaciones diferenciales que rigen la dinámica SIR
sir_foto_determinista <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Evolución temporal de los compartimentos epidemiológicos
    dS <- - (beta * S * I) / N + b * (I + R)
    dI <- (beta * S * I) / N - (b + gamma) * I
    dR <- (gamma * I) - (b * R)
    return(list(c(dS, dI, dR)))
  })
}

# Ejecución de la integración numérica mediante el algoritmo estándar de deSolve
out_det <- ode(y = initial_state, times = times, func = sir_foto_determinista, parms = params)
df_det <- as.data.frame(out_det)

# --- 4. VISUALIZACIÓN DE LA CONVERGENCIA AL EQUILIBRIO ---
# Configuración de los atributos estéticos para diferenciar trayectorias y horizontes estacionales
mis_colores <- c(
  "Susceptibles (S)" = "steelblue",   "Equilibrio S*" = "steelblue",
  "Infectados (I)"   = "firebrick",   "Equilibrio I*" = "firebrick",
  "Recuperados (R)"  = "forestgreen", "Equilibrio R*" = "forestgreen"
)

mis_tipos <- c(
  "Susceptibles (S)" = "solid",       "Equilibrio S*" = "dotted",
  "Infectados (I)"   = "solid",       "Equilibrio I*" = "dotted",
  "Recuperados (R)"  = "solid",       "Equilibrio R*" = "dotted"
)

# Generación de la figura técnica que muestra la evolución de los compartimentos
g_final <- ggplot(df_det, aes(x = time)) +
  # Representación de los niveles de equilibrio teórico (asíntotas)
  geom_hline(aes(yintercept = S_star, color = "Equilibrio S*", linetype = "Equilibrio S*"), size = 0.8) +
  geom_hline(aes(yintercept = I_star, color = "Equilibrio I*", linetype = "Equilibrio I*"), size = 0.8) +
  geom_hline(aes(yintercept = R_star, color = "Equilibrio R*", linetype = "Equilibrio R*"), size = 0.8) +
  
  # Trazado de las curvas de población obtenidas por integración
  geom_line(aes(y = S, color = "Susceptibles (S)", linetype = "Susceptibles (S)"), size = 1.2) +
  geom_line(aes(y = I, color = "Infectados (I)", linetype = "Infectados (I)"), size = 1.2) +
  geom_line(aes(y = R, color = "Recuperados (R)", linetype = "Recuperados (R)"), size = 1.2) +
  
  labs(title = NULL, subtitle = NULL, y = "Población", x = "Tiempo (t)", color = NULL, linetype = NULL) +
  scale_color_manual(values = mis_colores) +
  scale_linetype_manual(values = mis_tipos) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 113)) +
  
  # Configuración de la leyenda en formato de caja técnica
  theme(
    legend.position = c(0.82, 0.8), 
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm")
  )

print(g_final)