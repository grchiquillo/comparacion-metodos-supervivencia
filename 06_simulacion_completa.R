############################################################################

# 6. Simulación completa covariables + tiempo supervivencia + censura

#se incluye como riesgo basal el **beta0**, el tiempo al evento sigue una 
#Weibull al igual que la censura



simular_dataset_supervivencia <- function(n = 1000,
                                          seed = 123,
                                          tipo_relacion = c("lineal", "no_lineal"),
                                          censura_objetivo = 0.50,
                                          gamma_evento = 1,
                                          beta0,
                                          gamma_cens = 1,
                                          lambda_cens = NULL,
                                          admin_time = NULL) {
  
  tipo_relacion <- match.arg(tipo_relacion)
  set.seed(seed)
  
  # 1) covariables
  sim <- simular_covariables(n = n, seed = seed)
  
  # 2) tiempo verdadero al evento
  sim_evento <- simular_tiempo_evento(
    data = sim,
    tipo_relacion = tipo_relacion,
    gamma = gamma_evento,
    beta0 = beta0
  )
  
  T_true <- sim_evento$t_evento
  eta <- sim_evento$eta
  lambda_i <- sim_evento$lambda_i
  
  # 3) calibrar lambda de censura si no viene dado
  if (is.null(lambda_cens)) {
    lambda_cens <- calibrar_lambda_censura(
      data = sim,
      tipo_relacion = tipo_relacion,
      gamma_evento = gamma_evento,
      beta0 = beta0,
      gamma_cens = gamma_cens,
      censura_objetivo = censura_objetivo,
      admin_time = admin_time,
      seed = seed
    )
  }
  
  # 4) tiempo de censura Weibull
  C <- rweibull(
    n = n,
    shape = gamma_cens,
    scale = lambda_cens
  )
  
  # 5) tiempo observado y estado
  if (!is.null(admin_time)) {
    Y <- pmin(T_true, C, admin_time)
    delta <- as.integer(T_true <= C & T_true <= admin_time)
  } else {
    Y <- pmin(T_true, C)
    delta <- as.integer(T_true <= C)
  }
  
  sim$true_event_time <- T_true
  sim$censoring_time <- C
  sim$time <- Y
  sim$event <- delta
  
  sim$overall_survival_months <- Y
  sim$overall_survival_num <- delta
  sim$overall_survival <- factor(
    ifelse(delta == 1, "Muerto", "Vivo"),
    levels = c("Muerto", "Vivo")
  )
  
  attr(sim, "censura_lograda") <- mean(delta == 0)
  attr(sim, "tipo_relacion") <- tipo_relacion
  attr(sim, "beta0") <- beta0
  attr(sim, "lambda_cens") <- lambda_cens
  
  return(sim)
}


## 6.1 Prueba de la simulación

### 6.1.1 Simulación escenario n500_c0.5_lineal:

#- Tamaño de muestra: 500
#- Tipo de relación: Lineal
#- % censura: 50%

set.seed(123)

sim_superv_50 <- simular_dataset_supervivencia(
  n = 500,
  seed = 123,
  tipo_relacion = "lineal",
  censura_objetivo = 0.50,
  gamma_evento = 1,
  beta0 = beta0_calibrado,
  gamma_cens = 1,
  admin_time = max_obs
)

# Ver atributos guardados
attr(sim_superv_50, "censura_lograda")
attr(sim_superv_50, "tipo_relacion")
attr(sim_superv_50, "beta0")
attr(sim_superv_50, "lambda_cens")


### 6.1.2 Descripción general del conjunto de datos simulado

# Resumen numérico
summary(sim_superv_50[, sapply(sim_superv_50, is.numeric)])

# Resumen categórico
lapply(sim_superv_50[, sapply(sim_superv_50, is.factor)], table)


