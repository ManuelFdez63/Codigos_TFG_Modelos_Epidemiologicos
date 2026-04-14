# --- 1. CONFIGURACIÓN DEL ENTORNO Y DEFINICIÓN DE PARÁMETROS ---

library(GillespieSSA) # Implementación del algoritmo de Gillespie
library(ggplot2)      # Visualización avanzada de datos
library(ellipse)      # Generación de elipses de confianza basadas en la covarianza

# Se define un tamaño poblacional elevado (N=5000) para garantizar la convergencia
# hacia el límite y asegurar la validez de la aproximación gaussiana.
params <- c(
  beta = 0.8,  
  gamma = 0.5, 
  b = 0.05,    
  N = 5000     
)

# Extracción de parámetros para asegurar la estabilidad numérica de los cálculos
beta_val  <- unname(params["beta"])
gamma_val <- unname(params["gamma"])
b_val     <- unname(params["b"])
N_val     <- unname(params["N"])

# --- 2. FUNDAMENTOS TEÓRICOS: EQUILIBRIO Y APROXIMACIÓN POR DIFUSIÓN ---
# Cálculo del número básico de reproducción (R0) y equilibrios
R0 <- beta_val / (b_val + gamma_val)
phi_S <- 1 / R0
phi_I <- (b_val / (gamma_val + b_val)) * (1 - 1/R0)

# Coordenadas del equilibrio endémico en términos de individuos
S_eq <- N_val * phi_S
I_eq <- N_val * phi_I

# Cálculo de la matriz de difusión B en el equilibrio 
B_SS <- beta_val * phi_S * phi_I + b_val * (1 - phi_S)
B_SI <- -(beta_val * phi_S * phi_I + b_val * phi_I)
B_II <- beta_val * phi_S * phi_I + (gamma_val + b_val) * phi_I

# Resolución de la ecuación de Lyapunov para la matriz de covarianza 
# Se determinan las fluctuaciones de segundo orden alrededor del equilibrio
Sigma_SI <- -B_II / (2 * beta_val * phi_I)
Sigma_SS <- (B_SS - 2 * (gamma_val + b_val) * Sigma_SI) / (2 * (beta_val * phi_I + b_val))
Sigma_II <- (B_SI + beta_val * phi_I * Sigma_SS - (beta_val * phi_I + b_val) * Sigma_SI) / (gamma_val + b_val)

# Construcción de la matriz de covarianza escalada por el tamaño del sistema
Cov_Matrix <- N_val * matrix(c(Sigma_SS, Sigma_SI, 
                               Sigma_SI, Sigma_II), nrow = 2)

# Generación de la elipse de confianza (nivel de significación del 95%)
elipse_puntos <- as.data.frame(ellipse(x = Cov_Matrix, centre = c(S_eq, I_eq), level = 0.95))
colnames(elipse_puntos) <- c("S", "I")

# --- 3. IMPLEMENTACIÓN COMPUTACIONAL: SIMULACIÓN DE GILLESPIE ---

# Matriz de cambio y probabilidades de transición
nu <- matrix(c(
  -1, +1,  0,   
  0, -1, +1,   
  +1,  0,  0,   
  -1,  0,  0,   
  0, -1,  0,   
  0,  0, -1    
), nrow = 3, byrow = FALSE)

a <- c("beta * S * I / N", "gamma * I", "b * N", "b * S", "b * I", "b * R")

# Definición del estado inicial mediante redondeo al entero más cercano
initial_state_eq <- c(S = round(S_eq), I = round(I_eq), R = round(N_val - S_eq - I_eq))

set.seed(42) 
res <- ssa(x0 = initial_state_eq, a = a, nu = nu, parms = params, 
           tf = 500, method = ssa.d(), simName = "Fase")

# Procesamiento de resultados y filtrado del régimen transitorio
df_fase <- as.data.frame(res$data)
colnames(df_fase) <- c("time", "S", "I", "R")
df_fase_estacionario <- subset(df_fase, time > 100)

# --- 4. REPRESENTACIÓN EN EL PLANO DE FASES (S, I) ---
g_fases <- ggplot() +
  # Representación de las realizaciones estocásticas (nube de puntos)
  geom_point(data = df_fase_estacionario, aes(x = S, y = I, color = "Simulación (Gillespie)"), 
             alpha = 0.3, size = 1.2) +
  # Superposición del contorno teórico de la distribución normal bivariada
  geom_path(data = elipse_puntos, aes(x = S, y = I, color = "Elipse Teórica"), 
            linewidth = 1.2) +
  # Localización del punto de equilibrio determinista
  geom_point(aes(x = S_eq, y = I_eq, color = "Equilibrio Determinista"), 
             shape = 4, size = 4, stroke = 1.5) +
  
  labs(title = NULL, x = "Susceptibles (S)", y = "Infectados (I)", color = NULL) +
  theme_bw() +
  
  # Configuración de escalas cromáticas y personalización de la leyenda
  scale_color_manual(values = c("Simulación (Gillespie)" = "steelblue", 
                                "Elipse Teórica" = "firebrick",
                                "Equilibrio Determinista" = "black")) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.35, "cm"),
    legend.margin = margin(t=3, r=3, b=3, l=3),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

print(g_fases)