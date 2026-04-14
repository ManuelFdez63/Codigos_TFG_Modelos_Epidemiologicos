# ==============================================================================
# SCRIPT: REPRESENTACIÓN TRIDIMENSIONAL DE LA DISTRIBUCIÓN CUASI-ESTACIONARIA
# Análisis de la superficie de densidad de probabilidad mediante la LNA
# ==============================================================================

# --- 1. GESTIÓN DE LIBRERÍAS Y DEPENDENCIAS ---
# En caso de ausencia de los paquetes en el entorno local, se requiere su instalación:
# install.packages("mvtnorm")
# install.packages("rgl")

library(mvtnorm) # Evaluación de la densidad de probabilidad normal bivariada
library(rgl)     # Generación de gráficos tridimensionales interactivos

# --- 2. PARÁMETROS Y CÁLCULO DE LA MATRIZ DE COVARIANZA (VAN KAMPEN) ---
# Definición del sistema epidemiológico bajo la Linear Noise Approximation (LNA)
params <- c(beta = 0.8, gamma = 0.5, b = 0.05, N = 5000)

# Extracción de escalares para el cálculo de los momentos de la distribución
beta_val  <- unname(params["beta"])
gamma_val <- unname(params["gamma"])
b_val     <- unname(params["b"])
N_val     <- unname(params["N"])

# Determinación de los puntos de equilibrio endémico del modelo determinista
R0 <- beta_val / (b_val + gamma_val)
phi_S <- 1 / R0
phi_I <- (b_val / (gamma_val + b_val)) * (1 - 1/R0)

S_eq <- N_val * phi_S
I_eq <- N_val * phi_I

# Cálculo de los elementos de la matriz de difusión B en el estado estacionario
B_SS <- beta_val * phi_S * phi_I + b_val * (1 - phi_S)
B_SI <- -(beta_val * phi_S * phi_I + b_val * phi_I)
B_II <- beta_val * phi_S * phi_I + (gamma_val + b_val) * phi_I

# Resolución de la ecuación de Lyapunov para la obtención de las covarianzas
Sigma_SI <- -B_II / (2 * beta_val * phi_I)
Sigma_SS <- (B_SS - 2 * (gamma_val + b_val) * Sigma_SI) / (2 * (beta_val * phi_I + b_val))
Sigma_II <- (B_SI + beta_val * phi_I * Sigma_SS - (beta_val * phi_I + b_val) * Sigma_SI) / (gamma_val + b_val)

# Construcción de la matriz de covarianza escalada al tamaño poblacional N
Cov_Matrix <- N_val * matrix(c(Sigma_SS, Sigma_SI, 
                               Sigma_SI, Sigma_II), nrow = 2)

# --- 3. GENERACIÓN DEL DOMINIO Y EVALUACIÓN DE LA DENSIDAD ---
# Estimación de las desviaciones típicas para la delimitación del área de influencia
sd_S <- sqrt(Cov_Matrix[1, 1])
sd_I <- sqrt(Cov_Matrix[2, 2])

# Definición de la malla espacial (+/- 3 desviaciones estándares respecto al equilibrio)
S_seq <- seq(S_eq - 3*sd_S, S_eq + 3*sd_S, length.out = 50)
I_seq <- seq(I_eq - 3*sd_I, I_eq + 3*sd_I, length.out = 50)

# Inicialización de la matriz de densidad de probabilidad Z
Z_matrix <- matrix(0, nrow = length(S_seq), ncol = length(I_seq))

# Evaluación de la función de densidad de probabilidad normal multivariante
for (i in 1:length(S_seq)) {
  for (j in 1:length(I_seq)) {
    punto_actual <- c(S_seq[i], I_seq[j])
    Z_matrix[i, j] <- dmvnorm(x = punto_actual, 
                              mean = c(S_eq, I_eq), 
                              sigma = Cov_Matrix)
  }
}

# --- 4. VISUALIZACIÓN TRIDIMENSIONAL (INTERFAZ RGL) ---
# Inicialización del dispositivo gráfico rgl
open3d()

# Generación de la superficie de densidad (Campana de Gauss bivariada)
persp3d(x = S_seq, 
        y = I_seq, 
        z = Z_matrix, 
        col = "lightblue",      
        xlab = "Susceptibles (S)", 
        ylab = "Infectados (I)", 
        zlab = "Probabilidad Pi(S,I)",
        box = TRUE, axes = TRUE)

# Superposición de la malla alámbrica para optimizar la percepción del volumen
surface3d(x = S_seq, y = I_seq, z = Z_matrix, color = "black", front = "lines", back = "lines")