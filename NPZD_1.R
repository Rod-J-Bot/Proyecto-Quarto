#######################################################################
# Modelo NPZD a desarrollar para lago Llanquihue                      #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#######################################################################

# Option 2 para resolver el modelo

# OD NPZD model in R using ode function drom deSolve
install.packages('deSolve')
require(deSolve)

#variable inicial de estado  [nmolN/m3]
state <- c(DIN = 10.0, PHYTO = 0.5, ZOO = 0.3, DET = 5.0)

# Tabla de datos
tabla_datos <- data.frame(
  Nombre = c("Dato1", "Dato2", "Dato3"),
  Valor = c(1.0, 2.0, 3.0)
)

# Parametros
parms <- c(
  maxUptake  = 1.0, 1.5,  #[/day]
  ksPAR      = 140, 155,  # [uEinst/m2/s]
  ksDIN      = 1.0, 1.5,  # [nmolN/m3]
  maxGrazing = 1.0, 1.5,  # [/day]
  ksPHYTO    = 1.0, 1.5,   # [nmolN/m3]
  pFaeces    = 0.3, 0.6,  # [-]
  excretionRate = 0.1, 0.2, # [/day]
  mortalityRate = 0.4, 0.8,  # [/nmolN/m3/day]
  mineralizationRate = 0.05, 0.08, # [/day]
  depth              = 10, 15,  # [m]
  K                  = 20, 22, # [ug/l] P mortality half saturation constant
  u                 =  0.028, 0.035, # [d^-1) P mortality rate
  tablaDatos = tabla_datos  # Agregar la tabla de datos como un parámetro
)


# formulacion del modelo
NPZD <- function(t,state,parameters)
{
  with(as.list(c(state,parameters)),{
    
    # Forcing function: Luz es una funcion Seno
    PAR <- 0.5*(540+440*sin(2*pi*(t-81)/365))*exp(-0.05*depth/2)
    
    # Rate expressions [nmolN/m3/day]
    N_Uptake       <- maxUptake * min(PAR/(PAR +ksPAR),DIN/(DIN+ksDIN)) * PHYTO
    Mortality_P    <- ((u*PHYTO)/(K+PHYTO)) * PHYTO
    Grazing        <- maxGrazing * PHYTO/(PHYTO+ksPHYTO) * ZOO
    Faeces         <- Grazing * pFaeces
    Excretion      <- excretionRate * ZOO
    Mortality      <- mortalityRate * ZOO * ZOO
    Mineralization <- mineralizationRate * DET
    
    # Balance de masas [nmolN/m3/day]
    dPHYTO <- N_Uptake - Grazing - Mortality_P
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
outtimes <-seq(from = 0, to = 2*365,length.out=30)

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