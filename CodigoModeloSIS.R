# --- 1. CONFIGURACIÓN INICIAL ---
library(deSolve)      # Para resolver Ecuaciones Diferenciales (Determinista)
library(GillespieSSA) # Para el algoritmo estocástico (Simulación exacta)
library(ggplot2)      # Para gráficos
library(reshape2)     # Para manipular datos para ggplot

# Parámetros del modelo 
params <- c(
  beta = 1.0,       # Tasa de transmisión
  mu = 0.3,         # Tasa de recuperación
  b = 0.1,          # Tasa de nacimiento/muerte
  N = 100           # Población total fija
)

# Condiciones iniciales
initial_state <- c(S = 95, I = 5)

# Tiempo de simulación e incremento
times <- seq(0, 50, by = 0.1)

# Cálculo del R0
R0 <- params["beta"] / (params["b"] + params["mu"])
print(paste("El Número Básico de Reproducción (R0) es:", round(R0, 2)))



# --- 2. MODELO DETERMINISTA ---

# Definición del sistema de EDOs
sis_determinista <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - (beta * S * I) / N + (mu + b) * I
    dI <- (beta * S * I) / N - (mu + b) * I
    return(list(c(dS, dI)))
  })
}

# Resolución numérica del sistema 
out_det <- ode(y = initial_state, times = times, func = sis_determinista, parms = params)
df_det <- as.data.frame(out_det)


# --- GRÁFICO DETERMINISTA ---
# 1. Preparación de etiquetas 
valor_equilibrio <- params["N"] * (1 - 1/R0)
label_equilibrio <- paste("Equilibrio Teórico =", round(valor_equilibrio, 2), "personas")


mis_colores <- c("Infectados (I)" = "firebrick", "Susceptibles (S)" = "steelblue")
mis_colores[label_equilibrio] <- "gray40" 

mis_tipos <- c("Infectados (I)" = "solid", "Susceptibles (S)" = "solid")
mis_tipos[label_equilibrio] <- "dashed" 


# 2. Generación del Gráfico
g1 <- ggplot(df_det, aes(x = time)) +
  geom_line(aes(y = I, color = "Infectados (I)", linetype = "Infectados (I)"), size = 1.2) +
  geom_line(aes(y = S, color = "Susceptibles (S)", linetype = "Susceptibles (S)"), size = 1.2) +
  
  geom_hline(aes(yintercept = valor_equilibrio, 
                 color = label_equilibrio, 
                 linetype = label_equilibrio), 
             size = 0.8) + 
  labs(title = NULL,
       subtitle = NULL,
       y = "Población", 
       x = "Tiempo (t)", 
       color = NULL,    
       linetype = NULL) +
  
  theme_bw() +
  scale_color_manual(values = mis_colores) +
  scale_linetype_manual(values = mis_tipos) +
  
  theme(
    legend.position = c(0.8, 0.85), 
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.5, "cm"),
    legend.margin = margin(4, 4, 4, 4)
  )

print(g1)




# --- 3. MODELO ESTOCÁSTICO (Algoritmo de Gillespie) ---

# Se definen los vectores de cambio de estado
nu <- matrix(c(-1, +1,   # Cambio en S
               +1, -1),  # Cambio en I
             nrow = 2, byrow = TRUE)

# Definición de las probabilidades de transición
a <- c("beta * S * I / N", "(mu + b) * I")

# Se ejecutan múltiples simulaciones
n_sims <- 20
set.seed(123)

lista_sims <- list()

for(i in 1:n_sims){
  res <- ssa(
    x0 = initial_state,
    a = a,
    nu = nu,
    parms = params,
    tf = 50,
    method = ssa.d(), 
    simName = "SIS"
  )
  
  temp <- as.data.frame(res$data)
  if(ncol(temp) >= 3) colnames(temp) <- c("time", "S", "I")
  
  temp$simulacion <- factor(i) 
  lista_sims[[i]] <- temp
}

df_stoch <- do.call(rbind, lista_sims)

# --- 4. COMPARATIVA Y VISUALIZACIÓN ---

# 1. Etiquetas y colores
val_equil <- params["N"] * (1 - 1/R0)
txt_equil <- paste("Equilibrio Teórico =", round(val_equil, 2), "personas")

# Nombres para las categorías
cat_stoch <- "Simulaciones Estocásticas (Gillespie)"
cat_det   <- "Media Determinista (EDO)"

# Colores y tipos de línea
mis_colores <- c()
mis_colores[cat_stoch] <- "firebrick"
mis_colores[cat_det]   <- "black"
mis_colores[txt_equil] <- "blue"

mis_tipos <- c()
mis_tipos[cat_stoch] <- "solid"
mis_tipos[cat_det]   <- "solid"
mis_tipos[txt_equil] <- "dashed" 


# 2. Generación del gráfico
g2 <- ggplot() +
  # Capa 1: Estocástica 
  geom_step(data = df_stoch, aes(x = time, y = I, group = simulacion, 
                                 color = cat_stoch, linetype = cat_stoch), 
            alpha = 0.3, size = 0.5) +
  
  # Capa 2: Determinista
  geom_line(data = df_det, aes(x = time, y = I, 
                               color = cat_det, linetype = cat_det), 
            size = 1.1) +
  
  # Capa 3: Equilibrio
  geom_hline(aes(yintercept = val_equil, 
                 color = txt_equil, linetype = txt_equil), 
             size = 0.8) +
  
  # Etiquetas
  labs(title = NULL,
       subtitle = NULL,
       x = "Tiempo (t)", 
       y = "Número de Infectados (I(t))",
       color = NULL, linetype = NULL) +
  theme_bw() +
  
  # Configuración de Colores y Líneas
  scale_color_manual(values = mis_colores) +
  scale_linetype_manual(values = mis_tipos) +
  
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1, linewidth = 1))) +
  theme(
    legend.position = c(0.80, 0.25), 
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.5, "cm"),
    legend.margin = margin(4, 4, 4, 4)
  )

print(g2)

