## Option 3 para resolver el modelo

renv::deactivate()

# OD NPZD model in R using ode function from deSolve
install.packages("deSolve")
require(deSolve)
renv::deactivate()

renv::install("deSolve")

library(deSolve)

#variable inicial de estado  [nmolN/m3]
state <- c(DIN = 10.0, PHYTO = 0.5, ZOO = 0.3, DET = 5.0)

# Parámetros
# Definimos la tabla de parámetros
params_table <- data.frame(
  Parameter = c("maxUptake", "ksPAR", "ksDIN", "maxGrazing", "ksPHYTO", "pFaeces", "excretionRate", "mortalityRate", "mineralizationRate", "depth" , "K", "u"),
  Sitio1    = c(1.0, 140, 1.0, 1.0, 1.0, 0.3, 0.1, 0.4, 0.05, 10, 20, 0.028),
  Sitio2    = c(1.5, 155, 1.5, 1.5, 1.5, 0.6, 0.2, 0.8, 0.08, 15, 22, 0.035),
  Sitio3    = c(1.0, 140, 1.0, 1.0, 1.0, 0.3, 0.1, 0.4, 0.05, 10, 20, 0.028),
  depth     = c(10, 15, 20)  # Profundidades correspondientes a cada sitio
)

print(params_table)

# Formulacion del modelo
NPZD <- function(t, state, params_table)
{
  with(as.list(c(state, params_table)), {
    
    # Forcing function: Luz es una funcion Seno
    PAR <- 0.5 * (540 + 440 * sin(2 * pi * (t - 81) / 365)) * exp(-0.05 * depth / 2)
    
    # Rate expressions [nmolN/m3/day]
    N_Uptake       <- maxUptake * min(PAR / (PAR + ksPAR), DIN / (DIN + ksDIN)) * PHYTO
    Mortality_P    <- ((u * PHYTO) / (K + PHYTO)) * PHYTO
    Grazing        <- maxGrazing * PHYTO / (PHYTO + ksPHYTO) * ZOO
    Faeces         <- Grazing * pFaeces
    Excretion      <- excretionRate * ZOO
    Mortality      <- mortalityRate * ZOO * ZOO
    Mineralization <- mineralizationRate * DET
    
    # Balance de masas [nmolN/m3/day]
    dPHYTO <- N_Uptake - Grazing - Mortality_P
    dZOO   <- Grazing - Faeces - Excretion - Mortality
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

NPZD

# Run the model

# Salida para 2 años
outtimes <- seq(from = 0, to = 2 * 365, length.out = 30)

# Solucion del modelo usando el paquete deSolve 
out <- ode(y = state, parms = params_table["Sitio1"], func = NPZD, times = outtimes)
out <- ode(y = state, parms = as.list(params_table["Sitio1", ]), func = NPZD, times = outtimes)


params_sitio1 <- params_table[, c("maxUptake", "ksPAR", "ksDIN", "maxGrazing", "ksPHYTO", "pFaeces", "excretionRate", "mortalityRate", "mineralizationRate", "depth", "K", "u")]
out <- ode(y = state, parms = params_sitio1, func = NPZD, times = outtimes)



out <- ode(y = state, parms = as.list(params_sitio1), func = NPZD, times = outtimes)



# Visualizar la salida
plot(out)
