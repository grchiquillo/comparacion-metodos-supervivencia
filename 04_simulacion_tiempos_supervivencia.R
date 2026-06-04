############################################################################
# 4. Simulación tiempos de supervivencia

## 4.1 Estandarización de las variables continuas usando parámetros de los datos reales

edad_mean  <- mean(datos1$age_at_diagnosis, na.rm = TRUE)
edad_sd    <- sd(datos1$age_at_diagnosis, na.rm = TRUE)

nodes_mean <- mean(datos1$lymph_nodes_examined_positive, na.rm = TRUE)
nodes_sd   <- sd(datos1$lymph_nodes_examined_positive, na.rm = TRUE)

size_mean  <- mean(datos1$tumor_size, na.rm = TRUE)
size_sd    <- sd(datos1$tumor_size, na.rm = TRUE)


## 4.2 Simulación del tiempo bajo una distribución Weibull
### 4.2.1 Supervivencia en la METABRIC original

#Para saber que parámetro lamda usar y que tenga coherencia con el 
#comportamiento de la base de datos real METABRIC se evalua el 
#comportamiento real de la supervivencia


metabric_surv <- datos1 %>%
  transmute(
    time = as.numeric(overall_survival_months),
  
    status = ifelse(as.integer(overall_survival) ==1, 1, 0)
  ) %>%
  filter(!is.na(time), !is.na(status), time > 0)


#### 4.2.1.1 calibración la Weibull basal con la mediana Kaplan-Meier

km_fit <- survfit(Surv(time, status) ~ 1, data = metabric_surv)

resumen_km <- summary(km_fit)$table
mediana_km <- as.numeric(resumen_km["median"])
max_obs <- max(metabric_surv$time, na.rm = TRUE)

mediana_km
max_obs

##### 4.2.1.1.1 Gráfica de supervivencia en la cohorte METABRIC

dir.create("figuras_anexos", showWarnings = FALSE, recursive = TRUE)

base_km <- datos1 %>%
  mutate(
    tiempo = as.numeric(overall_survival_months),
    
    evento = as.numeric(overall_survival)
  ) %>%
  filter(!is.na(tiempo), !is.na(evento), tiempo > 0)

# Ajuste Kaplan-Meier
km_fit <- survfit(Surv(tiempo, evento) ~ 1, data = base_km)

# Mediana de supervivencia observada
mediana_supervivencia <- summary(km_fit)$table["median"]

mediana_supervivencia


# Extraer información de la curva Kaplan-Meier
km_sum <- summary(km_fit)

df_km <- data.frame(
  tiempo = km_sum$time,
  supervivencia = km_sum$surv,
  lower = km_sum$lower,
  upper = km_sum$upper
)

# Figura S1
figura_s1 <- ggplot(df_km, aes(x = tiempo, y = supervivencia)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = mediana_supervivencia, linetype = "dashed") +
  annotate(
    "text",
    x = mediana_supervivencia,
    y = 0.55,
    label = paste0("Mediana = ", round(mediana_supervivencia, 1), " meses"),
    hjust = -0.05,
    size = 3.5
  ) +
  labs(
    title = "Distribución de supervivencia en la cohorte METABRIC",
    x = "Tiempo de supervivencia global (meses)",
    y = "Probabilidad de supervivencia estimada"
  ) +
  theme_minimal(base_size = 12)

# Mostrar figura
figura_s1

# Guardar figura para anexos
ggsave(
  filename = "figuras_anexos/Figura_S1_KM_METABRIC.png",
  plot = figura_s1,
  width = 8,
  height = 5,
  dpi = 300
)



#### 4.2.1.2 lamda cero y beta cero calibrados

gamma_evento <- 1

lambda0_calibrada <- mediana_km / (log(2))^(1 / gamma_evento)
beta0_calibrado <- -gamma_evento * log(lambda0_calibrada)

lambda0_calibrada
beta0_calibrado



### 4.2.2 Función para construir el predictor lineal del riesgo (de las covariables)

#Una vez simuladas y estandarizadas las covariables, se procede a 
#la construcción del predictor lineal del riesgo, el cual resume el efecto 
#conjunto de las covariables sobre el tiempo al evento.

#La función construir_fx devuelve la parte explicada por las variables.


construir_fx <- function(data, tipo_relacion = c("lineal", "no_lineal")) {
  tipo_relacion <- match.arg(tipo_relacion)
  
  edad_z  <- (data$age_at_diagnosis - mean(data$age_at_diagnosis)) / sd(data$age_at_diagnosis)
  nodes_z <- (data$lymph_nodes_examined_positive - mean(data$lymph_nodes_examined_positive)) / sd(data$lymph_nodes_examined_positive)
  size_z  <- (data$tumor_size - mean(data$tumor_size)) / sd(data$tumor_size)
  
  her2_num  <- ifelse(data$her2_status == "Positivo", 1, 0)
  chemo_num <- ifelse(data$chemotherapy == "Si", 1, 0)
  surg_num  <- ifelse(data$type_of_breast_surgery == "Mastectomia", 1, 0)
  
  grade2 <- ifelse(data$neoplasm_histologic_grade == "Grado 2", 1, 0)
  grade3 <- ifelse(data$neoplasm_histologic_grade == "Grado 3", 1, 0)
  
  if (tipo_relacion == "lineal") {
    fx <-
      0.25 * edad_z +
      0.35 * nodes_z +
      0.30 * her2_num +
      0.20 * size_z +
      0.25 * grade2 +
      0.45 * grade3 +
      0.15 * chemo_num +
      0.10 * surg_num
  }
  
  if (tipo_relacion == "no_lineal") {
    fx <-
      0.20 * edad_z +
      0.15 * (edad_z^2) +
      0.05 * (edad_z^3) +
      0.30 * log(data$lymph_nodes_examined_positive + 1) +
      0.30 * her2_num +
      0.15 * size_z +
      0.20 * (size_z^2) +
      0.25 * grade2 +
      0.45 * grade3 +
      0.10 * chemo_num +
      0.08 * surg_num +
      0.18 * her2_num * size_z
  }
  
  return(fx)
}



### 4.2.3 Construcción del predictor lineal completo

construir_eta <- function(data,
                          tipo_relacion = c("lineal", "no_lineal"),
                          beta0) {
  tipo_relacion <- match.arg(tipo_relacion)
  fx <- construir_fx(data, tipo_relacion = tipo_relacion)
  eta <- beta0 + fx
  return(eta)
}



## 4.3 Construcción función para simular tiempo al evento

#A continuación se construye la función simular_tiempo_evento, 
#la cual genera para individuo un tiempo al evento a partir de una 
#distribución weibull.

#Esta función incluye el tipo de relacion (lineal o no lineal), 
#el gamma (parámetro de forma de la weibull) y el lambda (escala basal).


simular_tiempo_evento <- function(data,
                                  tipo_relacion = c("lineal", "no_lineal"),
                                  gamma = 1,
                                  beta0) {
  
  tipo_relacion <- match.arg(tipo_relacion)
  
  eta <- construir_eta(
    data = data,
    tipo_relacion = tipo_relacion,
    beta0 = beta0
  )
  
  lambda_i <- exp(-eta / gamma)
  
  t_evento <- rweibull(
    n = nrow(data),
    shape = gamma,
    scale = lambda_i
  )
  
  return(list(
    t_evento = t_evento,
    eta = eta,
    lambda_i = lambda_i
  ))
}



