#######################################################################
# Modelo NPZD a desarrollar para lago Llanquihue                      #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#######################################################################


# Hay varias maneras de abordar el tema, cada creativo

# Opcion 1 para resolver el modelo

# Instalar y cargar el paquete "deSolve"
install.packages("deSolve")
library(deSolve)


# Definir parámetros del modelo
g_max <- 0.5   # Tasa máxima de crecimiento del fitoplancton
m <- 0.1       # Tasa de mortalidad del zooplancton
l <- 0.2       # Tasa de excreción del zooplancton
k_n <- 2       # Tasa de remineralización de nutrientes
k_p <- 0.5     # Tasa de asimilación de nutrientes por el fitoplancton
k_z <- 0.1     # Tasa de asimilación de fitoplancton por el zooplancton

# Definir condiciones iniciales
N <- 10   # Concentración inicial de nutrientes
P <- 3    # Concentración inicial de fitoplancton
Z <- 1.5    # Concentración inicial de zooplancton
initial_state <- c(N, P, Z)

# Definir función de derivadas
derivatives <- function(time, state, parameters) {
  N <- state[1]  # Concentración de nutrientes
  P <- state[2]  # Concentración de fitoplancton
  Z <- state[3]  # Concentración de zooplancton
  
  dN <- -parameters$k_n * P + parameters$k_z * parameters$l * Z
  dP <- parameters$g_max * P * N / (parameters$k_n + N) - parameters$m * P - parameters$k_p * Z * P
  dZ <- parameters$k_p * Z * P - parameters$l * Z * (1 - exp(-parameters$k_z * P))
  
  return(list(c(dN, dP, dZ)))
}

# Definir tiempos de simulación
times <- seq(0, 200, by = 0.1)

# Realizar la simulación
output <- ode(y = initial_state, times = times, func = derivatives, parms = list(g_max = g_max, m = m, l = l, k_n = k_n, k_p = k_p, k_z = k_z))


# Graficar los resultados
plot(output, col = c("blue", "green", "red"), xlab = "Tiempo", ylab = "Concentración",
     main = "Simulación NPZ")

# Agregar leyenda
legend("topright", legend = c("Nutrientes", "Fitoplancton", "Zooplancton"), col = c("blue", "green", "red"), lty = 1)


# Convertir el objeto output en un data frame
output_df <- as.data.frame(output)

# Crear el gráfico para la concentración de DIN
din_plot <- ggplot(output_df, aes(x = tiempo, y = DIN, color = "DIN")) +
  geom_line() +
  labs(x = "Tiempo", y = "Concentración", title = "Simulación NPZD") +
  scale_color_manual(values = "blue")

# Crear el gráfico para la concentración de PHYTO
phyto_plot <- ggplot(output_df, aes(x = tiempo, y = PHYTO, color = "Fitoplancton")) +
  geom_line() +
  labs(x = "Tiempo", y = "Concentración", title = "Simulación NPZD") +
  scale_color_manual(values = "green")

# Crear el gráfico para la concentración de ZOO
zoo_plot <- ggplot(output_df, aes(x = tiempo, y = ZOO, color = "Zooplancton")) +
  geom_line() +
  labs(x = "Tiempo", y = "Concentración", title = "Simulación NPZD") +
  scale_color_manual(values = "red")

# Mostrar los gráficos en una única ventana
grid.arrange(din_plot, phyto_plot, zoo_plot, nrow = 2)

grid.arrange(plot_din, plot_phyto, plot_zoo, plot_det, nrow = 2)



#######################################################################
# Modelo NPZD a desarrollar para lago Llanquihue                      #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#######################################################################

renv::update

# Option 2 para resolver el modelo

# OD NPZD model in R using ode function drom deSolve
install.packages('deSolve')
require(deSolve)

#variable inicial de estado  [nmolN/m3]
state <- c(DIN = 10.0, PHYTO = 0.5, ZOO = 0.3, DET = 5.0)

# Parametros
parms <- c(
  maxUptake  = 1.0, #[/day]
  ksPAR      = 140, # [uEinst/m2/s]
  ksDIN      = 1.0, # [nmolN/m3]
  maxGrazing = 1.0, # [/day]
  ksPHYTO    = 1.0,  # [nmolN/m3]
  pFaeces    = 0.3,  # [-]
  excretionRate = 0.1, # [/day]
  mortalityRate = 0.4, # [/nmolN/m3/day]
  mineralizationRate = 0.05, # [/day]
  depth              = 10 # [m]
)


# formulacion del modelo
NPZD <- function(t,state,parameters)
{
  with(as.list(c(state,parameters)),{
    
   # Forcing function: Luz es una funcion Seno
    PAR <- 0.5*(540+440*sin(2*pi*(t-81)/365))*exp(-0.05*depth/2)
    
     # Rate expressions [nmolN/m3/day]
    N_Uptake       <- maxUptake * min(PAR/(PAR +ksPAR),DIN/(DIN+ksDIN)) * PHYTO
    Grazing        <- maxGrazing * PHYTO/(PHYTO+ksPHYTO) * ZOO
    Faeces         <- Grazing * pFaeces
    Excretion      <- excretionRate * ZOO
    Mortality      <- mortalityRate * ZOO * ZOO
    Mineralization <- mineralizationRate * DET
    
    # Balance de masas [nmolN/m3/day]
    dPHYTO <- N_Uptake - Grazing
    dZOO   <- Grazing - Faeces - Excretion  - Mortality
    dDET   <- Faeces + Mortality - Mineralization
    dDIN   <- Excretion + Mineralization - N_Uptake
 
    # Nitrogeno total
    TotalN <- PHYTO + ZOO + DET + DIN
    
    # Output
    return(list(c(dDIN, dPHYTO, dZOO, dDET),
                TotalN = TotalN, PAR = PAR)
           )
    
     })
}

# Run the model

# salida para 2 años
outtimes <-seq(from = 0, to = 2*365,length.out=100)

# Solucion del modelo usando el packete deSolve 
out <- ode(y = state,parms=parms, func=NPZD, times=outtimes)

# visualize output
plot(out)


library(ggplot2)

# Convertir los datos de salida a un dataframe
df <- as.data.frame(out)

# Crear el gráfico utilizando ggplot2
ggplot(df, aes(x = time)) +
  geom_line(aes(y = DIN, color = "DIN")) +
  geom_line(aes(y = PHYTO, color = "Fitoplancton")) +
  geom_line(aes(y = ZOO, color = "Zooplancton")) +
  geom_line(aes(y = DET, color = "Detritos")) +
  labs(x = "Tiempo", y = "Concentración [nmolN/m3]",
       title = "Simulación del modelo NPZD") +
  scale_color_manual(values = c("DIN" = "blue", "Fitoplancton" = "green",
                                "Zooplancton" = "red", "Detritos" = "orange")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = seq(0, 730, by = 100)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2))



library(ggplot2)

# Convertir los datos de salida a un dataframe
df <- as.data.frame(out)

# Gráfico de DIN
ggplot(df, aes(x = time, y = DIN)) +
  geom_line(color = "blue") +
  labs(x = "Tiempo", y = "Concentración DIN [nmolN/m3]",
       title = "Simulación del modelo NPZD - DIN") +
  theme_minimal()

# Gráfico de Fitoplancton
ggplot(df, aes(x = time, y = PHYTO)) +
  geom_line(color = "green") +
  labs(x = "Tiempo", y = "Concentración Fitoplancton [nmolN/m3]",
       title = "Simulación del modelo NPZD - Fitoplancton") +
  theme_minimal()

# Gráfico de Zooplancton
ggplot(df, aes(x = time, y = ZOO)) +
  geom_line(color = "red") +
  labs(x = "Tiempo", y = "Concentración Zooplancton [nmolN/m3]",
       title = "Simulación del modelo NPZD - Zooplancton") +
  theme_minimal()

# Gráfico de Detritos
ggplot(df, aes(x = time, y = DET)) +
  geom_line(color = "orange") +
  labs(x = "Tiempo", y = "Concentración Detritos [nmolN/m3]",
       title = "Simulación del modelo NPZD - Detritos") +
  theme_minimal()


library(ggplot2)
library(gridExtra)

# Convertir los datos de salida a un dataframe
df <- as.data.frame(out)

# Gráfico de DIN
plot_din <- ggplot(df, aes(x = time, y = DIN)) +
  geom_line(color = "blue") +
  labs(x = "Tiempo", y = "Concentración DIN [nmolN/m3]",
       title = "Simulación del modelo NPZD - DIN") +
  theme_minimal()

# Gráfico de Fitoplancton
plot_phyto <- ggplot(df, aes(x = time, y = PHYTO)) +
  geom_line(color = "green") +
  labs(x = "Tiempo", y = "Concentración Fitoplancton [nmolN/m3]",
       title = "Simulación del modelo NPZD - Fitoplancton") +
  theme_minimal()

# Gráfico de Zooplancton
plot_zoo <- ggplot(df, aes(x = time, y = ZOO)) +
  geom_line(color = "red") +
  labs(x = "Tiempo", y = "Concentración Zooplancton [nmolN/m3]",
       title = "Simulación del modelo NPZD - Zooplancton") +
  theme_minimal()

# Gráfico de Detritos
plot_det <- ggplot(df, aes(x = time, y = DET)) +
  geom_line(color = "orange") +
  labs(x = "Tiempo", y = "Concentración Detritos [nmolN/m3]",
       title = "Simulación del modelo NPZD - Detritos") +
  theme_minimal()

# Mostrar todos los gráficos en una misma ventana
grid.arrange(plot_din, plot_phyto, plot_zoo, plot_det, nrow = 2)



