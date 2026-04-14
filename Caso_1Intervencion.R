# --- 1. CONFIGURACIÓN DEL ENTORNO Y CARGA DE LIBRERÍAS ---
library(deSolve)      # Resolución de sistemas de ecuaciones diferenciales ordinarias
library(GillespieSSA) # Implementación de simulaciones estocásticas exactas
library(ggplot2)      # Generación de representaciones gráficas avanzadas
library(reshape2)     # Reestructuración y manipulación de conjuntos de datos
library(gridExtra)    # Disposición de múltiples objetos gráficos en paneles

# --- 2. DEFINICIÓN DE PARÁMETROS Y ESCENARIOS DE INTERVENCIÓN ---

# Especificación de los umbrales temporales y tamaño poblacional
t_intervencion <- 150
t_final <- 400
N_pob <- 100

# CONFIGURACIÓN FASE 1: Régimen de propagación libre (R0 > 1)
params_fase1 <- c(beta = 1, gamma = 0.2, b = 0.05, N = N_pob)

# --- ESCENARIOS DE INTERVENCIÓN (FASE 2) ---

# CASO 1: Restricción de movilidad (Reducción del parámetro de transmisión)
ghjsj <- c(beta = 0.2, gamma = 0.2, b = 0.05, N = 100)

# CASO 2: Aislamiento clínico y tratamiento (Incremento de la tasa de recuperación)
jfhj <- c(beta = 1.0, gamma = 1, b = 0.05, N = 100)

# CASO 3: Intervención combinada (Estrategia de supresión)
params_fase2 <- c(beta = 0.4, gamma = 0.5, b = 0.05, N = 100)

# Definición del estado inicial del sistema (Vector de estado x0)
x0 <- c(S = 90, I = 10, R = 0)

# Definición de la matriz estequiométrica y funciones de propensión para Gillespie
nu <- matrix(c(-1, 1, 0,  0, -1, 1,  1, 0, 0,  -1, 0, 0,  0, -1, 0,  0, 0, -1), 
             nrow = 3, byrow = FALSE)
a <- c("beta*S*I/N", "gamma*I", "b*N", "b*S", "b*I", "b*R")

# --- 3. SIMULACIÓN ESTOCÁSTICA MEDIANTE CONCATENACIÓN (STITCHING) ---
set.seed(123)
n_sims <- 20
lista_trayectorias <- list()

for(i in 1:n_sims) {
  
  # Ejecución de la Fase 1: Simulación hasta el instante de intervención
  res1 <- ssa(x0 = x0, a = a, nu = nu, parms = params_fase1, tf = t_intervencion, 
              method = ssa.d(), simName = "Fase1")
  df1 <- as.data.frame(res1$data)
  colnames(df1) <- c("time", "S", "I", "R")
  
  # Extracción del vector de estado final de la Fase 1 como condición inicial para Fase 2
  estado_final_fase1 <- unlist(df1[nrow(df1), c("S", "I", "R")])
  
  # Ejecución de la Fase 2: Simulación desde la intervención hasta el horizonte temporal final
  res2 <- ssa(x0 = estado_final_fase1, a = a, nu = nu, parms = params_fase2, 
              tf = (t_final - t_intervencion), method = ssa.d(), simName = "Fase2")
  df2 <- as.data.frame(res2$data)
  colnames(df2) <- c("time", "S", "I", "R")
  
  # Ajuste del eje temporal para asegurar la continuidad de la trayectoria
  df2$time <- df2$time + t_intervencion
  
  # Unión de los segmentos temporales en una única estructura de datos
  df_total <- rbind(df1, df2)
  df_total$simulacion <- factor(i)
  lista_trayectorias[[i]] <- df_total
}

df_stoch <- do.call(rbind, lista_trayectorias)

# --- 4. INTEGRACIÓN NUMÉRICA DEL MODELO DETERMINISTA HÍBRIDO ---

# Definición del sistema de ecuaciones diferenciales ordinarias (Modelo SIR con demografía)
sir_ode <- function(t, x, parms) {
  with(as.list(c(x, parms)), {
    dS <- b*N - beta*S*I/N - b*S
    dI <- beta*S*I/N - (b+gamma)*I
    dR <- gamma*I - b*R
    list(c(dS, dI, dR))
  })
}

# Integración numérica de la Fase 1
times1 <- seq(0, t_intervencion, by=0.1)
out1 <- ode(y = x0, times = times1, func = sir_ode, parms = params_fase1)
last_state_det <- out1[nrow(out1), -1] 

# Integración numérica de la Fase 2 basada en el estado terminal de la Fase 1
times2 <- seq(t_intervencion, t_final, by=0.1)
out2 <- ode(y = last_state_det, times = times2, func = sir_ode, parms = params_fase2)

df_det <- rbind(as.data.frame(out1), as.data.frame(out2))

# --- 5. REPRESENTACIÓN GRÁFICA Y COMPARATIVA DE RESULTADOS ---
df_stoch_long <- melt(df_stoch, id.vars = c("time", "simulacion"))
df_det_long   <- melt(df_det, id.vars = "time")

# Definición de función para la visualización por compartimentos individuales
plot_compartment <- function(comp, color_stoch, title) {
  ggplot() +
    # Delimitación visual del periodo de intervención (Sombreado)
    annotate("rect", xmin = t_intervencion, xmax = t_final, ymin = -Inf, ymax = Inf, 
             fill = "gray90", alpha = 0.5) +
    
    # Indicador de transición de fase
    geom_vline(xintercept = t_intervencion, linetype="dashed", color="gray40") +
    
    # Representación de las realizaciones estocásticas individuales
    geom_step(data = subset(df_stoch_long, variable == comp),
              aes(x = time, y = value, group = simulacion),
              color = color_stoch, alpha = 0.2, size = 0.3) +
    
    # Representación de la trayectoria determinista media
    geom_line(data = subset(df_det_long, variable == comp),
              aes(x = time, y = value),
              color = "black", linetype = "dashed", size = 0.8) +
    
    # Etiquetado de la fase de intervención
    annotate("text", x = t_intervencion + 2, y = Inf, label = "Intervención\n (Aislamiento + Confinamiento)", 
             hjust = 0, vjust = 1.5, size = 3, fontface="italic", color="gray30") +
    
    labs(title = title, x = "Tiempo (t)", y = "Población") +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 100)) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust=0.5))
}

# Construcción de paneles para S, I y R
gS <- plot_compartment("S", "steelblue", "Susceptibles (S)")
gI <- plot_compartment("I", "firebrick", "Infectados (I)")
gR <- plot_compartment("R", "forestgreen", "Recuperados (R)")

# Disposición final de la comparativa mediante grid.arrange
grid.arrange(gS, gI, gR, ncol = 3)