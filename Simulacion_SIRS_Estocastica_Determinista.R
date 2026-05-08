# --- 1. CONFIGURACIÓN DEL ENTORNO Y CARGA DE LIBRERÍAS ---
library(deSolve)      # Resolución numérica de sistemas de ecuaciones diferenciales
library(GillespieSSA) # Implementación de simulaciones estocásticas exactas
library(ggplot2)      # Generación de representaciones gráficas avanzadas
library(reshape2)     # Reestructuración de conjuntos de datos para visualización
library(gridExtra)    # Organización de composiciones gráficas multipanel

# Especificación de coeficientes y constantes del sistema
# Se define un escenario con R0 > 1 para observar el equilibrio endémico propiciado por el SIRS
params <- c(
  beta = 0.8,         # Coeficiente de transmisión
  gamma = 0.1,        # Tasa de recuperación (ajustado para R0 > 1)
  b = 0.05,           # Tasa de dinámica vital (natalidad/mortalidad)
  mu = 0.15,          # Tasa de pérdida de inmunidad 
  N = 100             # Tamaño de la población total
)

# Definición del vector de estado inicial
initial_state <- c(S = 90, I = 10, R = 0)

# Especificación del dominio temporal para la integración numérica
times <- seq(0, 100, by = 0.1)

# --- 2. EVALUACIÓN DEL NÚMERO BÁSICO DE REPRODUCCIÓN (R0) ---
R0 <- params["beta"] / (params["b"] + params["gamma"])
print(paste("El Número Básico de Reproducción (R0) es:", round(R0, 2)))

# --- 3. MODELADO DETERMINISTA MEDIANTE ECUACIONES DIFERENCIALES ---
# Formulación funcional del sistema de EDOs para el modelo SIRS con demografía
sirs_det_func <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Definición de las derivadas temporales de los compartimentos
    dS <- b * N - (beta * S * I) / N - b * S + mu * R
    dI <- (beta * S * I) / N - (b + gamma) * I
    dR <- (gamma * I) - (b * R) - mu * R
    return(list(c(dS, dI, dR)))
  })
}

# Ejecución de la integración numérica del sistema
out_det <- ode(y = initial_state, times = times, func = sirs_det_func, parms = params)
df_det <- as.data.frame(out_det)

# --- 4. MODELADO ESTOCÁSTICO (ALGORITMO DE GILLESPIE) ---
# Matriz de cambio de estado (transiciones) para el modelo SIRS
nu <- matrix(c(
  -1, +1,  0,   # 1. Infección (S -> I)
  0, -1, +1,   # 2. Recuperación (I -> R)
  +1,  0, -1,   # 3. Pérdida de inmunidad (R -> S) 
  +1,  0,  0,   # 4. Nacimiento (-> S)
  -1,  0,  0,   # 5. Mortalidad S
  0, -1,  0,   # 6. Mortalidad I
  0,  0, -1    # 7. Mortalidad R
), nrow = 3, byrow = FALSE)

# Vector de propensiones asociadas a cada evento
a <- c(
  "beta * S * I / N", # Tasa de infección
  "gamma * I",        # Tasa de recuperación
  "mu * R",           # Tasa de pérdida de inmunidad
  "b * N",            # Tasa de nacimiento
  "b * S",            # Tasa de mortalidad susceptible
  "b * I",            # Tasa de mortalidad infectado
  "b * R"             # Tasa de mortalidad recuperado
)

# Ejecución de simulaciones mediante el algoritmo de Gillespie
set.seed(42)
lista_sims <- list()
n_sims <- 20

for(i in 1:n_sims){
  res <- ssa(x0 = initial_state, a = a, nu = nu, parms = params, tf = 100, 
             method = ssa.d(), simName = "SIRS")
  temp <- as.data.frame(res$data)
  if(ncol(temp) >= 4) colnames(temp) <- c("time", "S", "I", "R")
  temp$simulacion <- factor(i)
  lista_sims[[i]] <- temp
}
df_stoch <- do.call(rbind, lista_sims)

# --- 5. REESTRUCTURACIÓN DE DATOS PARA ANÁLISIS ---
# Transformación a formato largo (long format) para procesamiento con ggplot2
df_stoch_long <- melt(df_stoch, id.vars = c("time", "simulacion"), 
                      measure.vars = c("S", "I", "R"),
                      variable.name = "Estado", value.name = "Poblacion")

df_det_long <- melt(df_det, id.vars = "time", 
                    measure.vars = c("S", "I", "R"),
                    variable.name = "Estado", value.name = "Poblacion")

# --- 6. IMPLEMENTACIÓN DE LA FUNCIÓN DE VISUALIZACIÓN ---
# Generación de gráficos comparativos por compartimento epidemiológico
crear_grafico <- function(estado_filtro, color_linea, titulo) {
  
  datos_stoch <- subset(df_stoch_long, Estado == estado_filtro)
  datos_det   <- subset(df_det_long, Estado == estado_filtro)
  
  ggplot() +
    # Representación de trayectorias estocásticas individuales
    geom_step(data = datos_stoch, 
              aes(x = time, y = Poblacion, group = simulacion, 
                  color = "Simulación (Gillespie)", linetype = "Simulación (Gillespie)"), 
              alpha = 0.2, size = 0.3) +
    
    # Representación de la trayectoria determinista de referencia
    geom_line(data = datos_det, 
              aes(x = time, y = Poblacion, 
                  color = "Media Determinista (EDO)", linetype = "Media Determinista (EDO)"), 
              size = 0.8) +
    
    labs(title = titulo, x = "Tiempo (t)", y = "Individuos", color = NULL, linetype = NULL) +
    theme_bw() +
    
    scale_color_manual(values = c("Simulación (Gillespie)" = color_linea, 
                                  "Media Determinista (EDO)" = "black")) +
    
    scale_linetype_manual(values = c("Simulación (Gillespie)" = "solid", 
                                     "Media Determinista (EDO)" = "dashed")) +
    
    # Ajuste del eje temporal y definición del dominio de visualización
    scale_x_continuous(breaks = seq(0, 100, by = 25)) +
    coord_cartesian(xlim = c(0, 100)) +
    
    theme(
      legend.position = c(0.70, 0.9), 
      legend.background = element_rect(fill = "white", color = "black", size = 0.3),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.4, "cm"),
      legend.margin = margin(3,3,3,3),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# --- 7. GENERACIÓN DE PANELES POR COMPARTIMENTO ---
g_S <- crear_grafico("S", "steelblue", "Susceptibles (S)")
g_I <- crear_grafico("I", "firebrick", "Infectados (I)")
g_R <- crear_grafico("R", "forestgreen", "Recuperados (R)")

# --- 8. COMPOSICIÓN GRÁFICA FINAL ---
grid.arrange(g_S, g_I, g_R, ncol = 3)
