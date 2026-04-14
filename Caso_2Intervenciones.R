# --- 1. CONFIGURACIÓN Y LIBRERÍAS ---
library(deSolve)      # Resolución de sistemas de EDOs (Determinista)
library(GillespieSSA) # Implementación del algoritmo de Gillespie (Estocástico)
library(ggplot2)      # Visualización de datos
library(reshape2)     # Reestructuración de dataframes
library(gridExtra)    # Disposición de múltiples gráficos

# --- 2. DEFINICIÓN DE ETAPAS Y PARÁMETROS ---

# Configuración temporal de las fases
t_fase2 <- 100  # Inicio de la segunda etapa
t_fase3 <- 250  # Inicio de la tercera etapa
t_final <- 400
N_pob   <- 100  # Tamaño de la población total

# PARÁMETROS FASE 1: Configuración inicial de propagación
params_fase1 <- c(beta = 1.0, gamma = 0.15, b = 0.05, N = N_pob)

# PARÁMETROS FASE 2: Configuración de mitigación intermedia
params_fase2 <- c(beta = 0.8, gamma = 0.15, b = 0.05, N = N_pob)

# PARÁMETROS FASE 3: Configuración final del escenario
params_fase3 <- c(beta = 2.4, gamma = 0.15, b = 0.05, N = N_pob)

# Estado inicial del sistema
x0 <- c(S = 9, I = 10, R = 0)

# Definición de la estructura estocástica (Matriz de cambio y probabilidades de transición)
nu <- matrix(c(-1, 1, 0,  0, -1, 1,  1, 0, 0,  -1, 0, 0,  0, -1, 0,  0, 0, -1), 
             nrow = 3, byrow = FALSE)
a <- c("beta*S*I/N", "gamma*I", "b*N", "b*S", "b*I", "b*R")

# --- 3. SIMULACIÓN ESTOCÁSTICA (ALGORITMO DE GILLESPIE) ---

set.seed(123)
n_sims <- 20
lista_trayectorias <- list()

for(i in 1:n_sims) {
  
  # Ejecución Fase 1
  res1 <- ssa(x0 = x0, a = a, nu = nu, parms = params_fase1, tf = t_fase2, 
              method = ssa.d(), simName = "Fase1")
  df1 <- as.data.frame(res1$data)
  colnames(df1) <- c("time", "S", "I", "R")
  
  # Transferencia de estados para Fase 2
  estado_fin_1 <- unlist(df1[nrow(df1), c("S", "I", "R")])
  
  # Ejecución Fase 2
  duracion_fase2 <- t_fase3 - t_fase2
  res2 <- ssa(x0 = estado_fin_1, a = a, nu = nu, parms = params_fase2, 
              tf = duracion_fase2, method = ssa.d(), simName = "Fase2")
  df2 <- as.data.frame(res2$data)
  colnames(df2) <- c("time", "S", "I", "R")
  df2$time <- df2$time + t_fase2
  
  # Transferencia de estados para Fase 3
  estado_fin_2 <- unlist(df2[nrow(df2), c("S", "I", "R")])
  
  # Ejecución Fase 3
  duracion_fase3 <- t_final - t_fase3
  res3 <- ssa(x0 = estado_fin_2, a = a, nu = nu, parms = params_fase3, 
              tf = duracion_fase3, method = ssa.d(), simName = "Fase3")
  df3 <- as.data.frame(res3$data)
  colnames(df3) <- c("time", "S", "I", "R")
  df3$time <- df3$time + t_fase3
  
  # Consolidación de datos de la simulación i
  df_total <- rbind(df1, df2, df3)
  df_total$simulacion <- factor(i)
  lista_trayectorias[[i]] <- df_total
}

df_stoch <- do.call(rbind, lista_trayectorias)

# --- 4. MODELADO DETERMINISTA ---

sir_ode <- function(t, x, parms) {
  with(as.list(c(x, parms)), {
    dS <- b*N - beta*S*I/N - b*S
    dI <- beta*S*I/N - (b+gamma)*I
    dR <- gamma*I - b*R
    list(c(dS, dI, dR))
  })
}

# Resolución secuencial de las ecuaciones diferenciales
out1 <- ode(y = x0, times = seq(0, t_fase2, by=0.1), func = sir_ode, parms = params_fase1)
last_1 <- out1[nrow(out1), -1]

out2 <- ode(y = last_1, times = seq(t_fase2, t_fase3, by=0.1), func = sir_ode, parms = params_fase2)
last_2 <- out2[nrow(out2), -1]

out3 <- ode(y = last_2, times = seq(t_fase3, t_final, by=0.1), func = sir_ode, parms = params_fase3)

df_det <- rbind(as.data.frame(out1), as.data.frame(out2), as.data.frame(out3))

# --- 5. VISUALIZACIÓN DE RESULTADOS ---

df_stoch_long <- melt(df_stoch, id.vars = c("time", "simulacion"))
df_det_long   <- melt(df_det, id.vars = "time")

# Función para la generación de gráficos por compartimento
plot_compartment <- function(comp, color_stoch, title) {
  ggplot() +
    annotate("rect", xmin = t_fase2, xmax = t_fase3, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.1) +
    annotate("rect", xmin = t_fase3, xmax = t_final, ymin = -Inf, ymax = Inf, fill = "green", alpha = 0.1) +
    
    # Delimitadores de transición
    geom_vline(xintercept = c(t_fase2, t_fase3), linetype="dashed", color="gray40") +
    
    # Capas de datos (Estocástico y Determinista)
    geom_step(data = subset(df_stoch_long, variable == comp),
              aes(x = time, y = value, group = simulacion),
              color = color_stoch, alpha = 0.15, size = 0.3) +
    
    geom_line(data = subset(df_det_long, variable == comp),
              aes(x = time, y = value),
              color = "black", linetype = "dashed", size = 0.8) +
    
    # Anotaciones de las medidas implementadas
    annotate("text", x = t_fase2 + 5, y = Inf, label = "Medidas\nSuaves", 
             hjust = 0, vjust = 1.5, size = 3, fontface="italic", color="gray30") +
    annotate("text", x = t_fase3 + 5, y = Inf, label = "Medidas\nDrásticas", 
             hjust = 0, vjust = 1.5, size = 3, fontface="italic", color="gray30") +
    
    labs(title = title, x = "Tiempo (t)", y = "Población") +
    coord_cartesian(xlim = c(0, t_final), ylim = c(0, 100)) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust=0.5))
}

# Composición final de la figura
gS <- plot_compartment("S", "steelblue", "Susceptibles (S)")
gI <- plot_compartment("I", "firebrick", "Infectados (I)")
gR <- plot_compartment("R", "forestgreen", "Recuperados (R)")

grid.arrange(gS, gI, gR, ncol = 3)
