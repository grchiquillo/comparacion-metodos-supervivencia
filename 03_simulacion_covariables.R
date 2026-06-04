############################################################

# 3. Simulación Covariables

#Se simulan las variables en cadena:
  
#age_at_diagnosis ~ Normal truncada.

#lymph_nodes_examined_positive ~ Regresión Poisson o binomial negativa 
#si hay sobredispersión con age_at_diagnosis.

#her2_status ~ Regresión logística binaria con age_at_diagnosis + lymph_nodes_examined_positive.

#tumor_size ~ Regresión lineal con age_at_diagnosis +
# lymph_nodes_examined_positive + her2_status.

#neoplasm_histologic_grade ~ Multinomial con age_at_diagnosis + 
#lymph_nodes_examined_positive + her2_status + tumor_size.

#chemotheraphy ~ Regresión logística binaria con age_at_diagnosis + 
#lymph_nodes_examined_positive + her2_status + tumor_size + neoplasm_histologic_grade

#type_of_breast_surgery ~ Regresión logística binaria con age_at_diagnosis + 
#lymph_nodes_examined_positive + her2_status + tumor_size + neoplasm_histologic_grade
#+ chemotheraphy.


## 3.1 Modelos

# age_at_diagnosis ~ Normal truncada
age_mu  <- mean(datos1$age_at_diagnosis, na.rm = TRUE)
age_sd  <- sd(datos1$age_at_diagnosis, na.rm = TRUE)
age_min <- min(datos1$age_at_diagnosis, na.rm = TRUE)
age_max <- max(datos1$age_at_diagnosis, na.rm = TRUE)


# lymph_nodes_examined_positive ~ lineal con age
m_nodes_pois <- glm(
  lymph_nodes_examined_positive ~ age_at_diagnosis,
  data = datos1,
  family = poisson(link = "log")
)

# Verificación simple de sobredispersión
dispersion_nodes <- sum(residuals(m_nodes_pois, type = "pearson")^2) / m_nodes_pois$df.residual
dispersion_nodes

# Si hay sobredispersión importante, usar binomial negativa
if (dispersion_nodes > 1.5) {
  m_nodes <- glm.nb(
    lymph_nodes_examined_positive ~ age_at_diagnosis,
    data = datos1
  )
  tipo_modelo_nodes <- "negbin"
} else {
  m_nodes <- m_nodes_pois
  tipo_modelo_nodes <- "poisson"
}




# her2_status
m_her2 <- glm(her2_status ~ age_at_diagnosis + lymph_nodes_examined_positive,data = datos1, family = binomial)

# tumor_size
m_size <- lm(
  tumor_size ~ age_at_diagnosis + lymph_nodes_examined_positive + her2_status,
  data = datos1
)

# neoplasm_histologic_grade
m_grade <- multinom(
  neoplasm_histologic_grade ~ age_at_diagnosis + lymph_nodes_examined_positive + her2_status + tumor_size,
  data = datos1,
  trace = FALSE
)

# chemotherapy
m_chemo <- glm(
  chemotherapy ~ age_at_diagnosis + lymph_nodes_examined_positive + her2_status +
    tumor_size + neoplasm_histologic_grade,
  data = datos1,
  family = binomial
)

# type_of_breast_surgery
m_surg <- glm(
  type_of_breast_surgery ~ age_at_diagnosis + lymph_nodes_examined_positive +
    her2_status + tumor_size + neoplasm_histologic_grade + chemotherapy,
  data = datos1,
  family = binomial
)


### 3.1.1 Modelo para variable número de nodos

m_nodes_pois
exp(coef(m_nodes_pois))

### 3.1.2 Modelo para variable her2

m_her2
exp(coef(m_her2))

### 3.1.3 Modelo para variable tamaño del tumor

m_size

### 3.1.4 Modelo para variable grado histológico
m_grade
exp(coef(m_grade))


### 3.1.5 Modelo para variable quimioterapia

m_chemo
exp(coef(m_chemo))

### 3.1.6 Modelo para variable tipo de cirugía

m_surg
exp(coef(m_surg))

levels(model.frame(m_surg)$type_of_breast_surgery) #Conocer niveles ref de esta variable



## 3.2 Validación de los modelos usados en la simulación

#Se realiza la validación de los modelos ajustados en el punto anterior 
#para verificar si dicho modelo representa adecuadamente la relación entre 
#las variables y si cumple los supuestos necesarios para su uso.

### 3.2.1 Funciones de validación

#Función para validar modelo poisson
validar_poisson <- function(modelo) {
  cat("Resumen del modelo\n")
  print(summary(modelo))
  
  cat("\nAIC:\n")
  print(AIC(modelo))
  
  cat("\nRazón Deviance / gl residual:\n")
  overdisp <- modelo$deviance / modelo$df.residual
  print(overdisp)
  
  cat("\nPrueba de sobredispersión:\n")
  print(dispersiontest(modelo))
  
  par(mfrow = c(1, 2))
  plot(fitted(modelo), residuals(modelo, type = "pearson"),
       xlab = "Valores ajustados", ylab = "Residuos de Pearson",
       main = "Poisson: ajustados vs residuos")
  abline(h = 0, lty = 2, col = "red")
  
  qqnorm(residuals(modelo, type = "deviance"),
         main = "QQ-plot residuos deviance")
  qqline(residuals(modelo, type = "deviance"), col = "red")
}


#Función para validar modelos logisticos binarios
validar_logistico <- function(modelo, y_real, positivo, corte = 0.5) {
  y_real <- factor(y_real)
  negativo <- setdiff(levels(y_real), positivo)[1]
  
  p <- predict(modelo, type = "response")
  pred <- ifelse(p >= corte, positivo, negativo)
  pred <- factor(pred, levels = levels(y_real))
  
  y_bin <- ifelse(y_real == positivo, 1, 0)
  
  cat("Resumen del modelo\n")
  print(summary(modelo))
  
  cat("\nAIC:\n")
  print(AIC(modelo))
  
  cat("\nHosmer-Lemeshow:\n")
  print(hoslem.test(y_bin, p, g = 10))
  
  cat("\nMatriz de confusión:\n")
  print(table(Real = y_real, Predicho = pred))
  
  cat("\nAccuracy:\n")
  print(mean(pred == y_real))
  
  roc_obj <- roc(y_bin, p)
  cat("\nAUC:\n")
  print(auc(roc_obj))
  
  plot(roc_obj, main = "Curva ROC")
}



#Función para validar modelo lineal
validar_lineal <- function(modelo) {
  cat("Resumen del modelo\n")
  print(summary(modelo))
  
  cat("\nAIC:\n")
  print(AIC(modelo))
  
  cat("\nVIF:\n")
  print(vif(modelo))
  
  cat("\nRMSE:\n")
  rmse <- sqrt(mean(residuals(modelo)^2))
  print(rmse)
  
  par(mfrow = c(2, 2))
  plot(modelo)
}



#Función para validar modelo multinomial
validar_multinom <- function(modelo, data, yvar) {
  real <- factor(data[[yvar]])
  pred <- predict(modelo, type = "class")
  
  cat("Resumen del modelo\n")
  print(summary(modelo))
  
  cat("\nAIC:\n")
  print(AIC(modelo))
  
  cat("\nMatriz de confusión:\n")
  print(table(Real = real, Predicho = pred))
  
  cat("\nAccuracy:\n")
  print(mean(pred == real))
  
  modelo_nulo <- multinom(as.formula(paste(yvar, "~ 1")),
                          data = data, trace = FALSE)
  
  pseudoR2 <- 1 - as.numeric(logLik(modelo) / logLik(modelo_nulo))
  
  cat("\nPseudo R2 de McFadden:\n")
  print(pseudoR2)
}



### 3.2.2 Validación de cada modelo

#### 3.2.2.1 Validación - número de nodos

validar_poisson(m_nodes_pois)
dispersiontest(m_nodes_pois) #prueba sobredispersión

#Prueba y comparación modelo Poisson vs Binomial negativo
m_nodes_nb <- glm.nb(lymph_nodes_examined_positive ~ age_at_diagnosis,
                     data = datos1)
summary(m_nodes_nb)

AIC(m_nodes_pois, m_nodes_nb)

#Al comparar el AIC del modelo de poisson con el modelo de regresión binomial
#negativa, se observa que el binomial negativo presenta un valor mucho menor
#lo que indica un mejor ajuste relativo del modelo a los datos. 
#Por lo anterior se cambia el modelo utilizado para la varaible numero
#de nodos positivos a una REGRESIÓN BINOMIAL NEGATIVA.



#### 3.2.2.2 Validación para modelo variable her2

validar_logistico(m_her2, datos1$her2_status, positivo = "Positivo")

#### 3.2.2.3 Validación para modelo variable tamaño tumor

validar_lineal(m_size)

#### 3.2.2.4 Validación para modelo variable grado histológico

validar_multinom(m_grade, datos1, "neoplasm_histologic_grade")


#### 3.2.2.5 Validación para modelo variable quimioterapia

validar_logistico(m_chemo, datos1$chemotherapy, positivo = "Si")


#### 3.2.2.6 Validación para modelo variable tipo cirugia

validar_logistico(m_surg, datos1$type_of_breast_surgery, positivo = "Mastectomia")



## 3.3 Funciones auxiliares para simular variables

#Se incluyen las siguientes funciones auxiliares para simular las variables 
#con el fin de replicar de mejor forma la estructura estadistica de 
#los datos originales.


# Funciones auxiliares para simular variables a partir de modelos ajustados

# Binaria desde glm binomial
sim_bin_from_glm <- function(modelo, newdata, niveles = NULL) {
  p <- predict(modelo, newdata = newdata, type = "response")
  y <- rbinom(n = nrow(newdata), size = 1, prob = p)
  
  if (!is.null(niveles)) {
    y <- factor(y, levels = c(0, 1), labels = niveles)
  }
  
  return(y)
}


# Continua desde lm
sim_cont_from_lm <- function(modelo, newdata) {
  mu <- predict(modelo, newdata = newdata)
  sigma <- summary(modelo)$sigma
  y <- rnorm(n = nrow(newdata), mean = mu, sd = sigma)
  return(y)
}


# Multinomial desde multinom
sim_multinom <- function(modelo, newdata, niveles = NULL) {
  probs <- predict(modelo, newdata = newdata, type = "probs")
  
  if (is.vector(probs)) {
    probs <- cbind(1 - probs, probs)
  }
  
  y <- apply(probs, 1, function(p) sample(seq_along(p), size = 1, prob = p))
  
  if (!is.null(niveles)) {
    y <- factor(y, levels = seq_along(niveles), labels = niveles)
  }
  
  return(y)
}


# Conteo desde Poisson o binomial negativa
sim_count_from_model <- function(modelo, newdata, tipo = c("poisson", "negbin")) {
  tipo <- match.arg(tipo)
  
  mu <- predict(modelo, newdata = newdata, type = "response")
  
  if (tipo == "poisson") {
    y <- rpois(n = nrow(newdata), lambda = mu)
  }
  
  if (tipo == "negbin") {
    # En glm.nb del paquete MASS, theta controla la dispersión
    theta <- modelo$theta
    y <- rnbinom(n = nrow(newdata), mu = mu, size = theta)
  }
  
  return(y)
}



## 3.4 Función para simular variables

### 3.4.1 Guardar niveles originales de las variables categóricas
niveles_her2  <- levels(datos1$her2_status)
niveles_grade <- levels(datos1$neoplasm_histologic_grade)
niveles_chemo <- levels(datos1$chemotherapy)
niveles_surg  <- levels(datos1$type_of_breast_surgery)


### 3.4.2 Función de simulación de covariables

simular_covariables <- function(n = 1000, seed = 123) {
  set.seed(seed)
  
  sim <- data.frame(
    age_at_diagnosis = rtruncnorm(
      n = n,
      a = age_min,
      b = age_max,
      mean = age_mu,
      sd = age_sd
    )
  )
  
  # lymph_nodes_examined_positive
  sim$lymph_nodes_examined_positive <- sim_count_from_model(
    modelo = m_nodes_nb,
    newdata = sim,
    tipo = tipo_modelo_nodes
  )
  
  # her2_status
  sim$her2_status <- sim_bin_from_glm(
    modelo = m_her2,
    newdata = sim,
    niveles = niveles_her2
  )
  
  # tumor_size
  sim$tumor_size <- sim_cont_from_lm(
    modelo = m_size,
    newdata = sim
  )
  sim$tumor_size <- pmax(0.1, sim$tumor_size)
  
  # neoplasm_histologic_grade
  sim$neoplasm_histologic_grade <- sim_multinom(
    modelo = m_grade,
    newdata = sim,
    niveles = niveles_grade
  )
  
  # chemotherapy
  sim$chemotherapy <- sim_bin_from_glm(
    modelo = m_chemo,
    newdata = sim,
    niveles = niveles_chemo
  )
  
  # type_of_breast_surgery
  sim$type_of_breast_surgery <- sim_bin_from_glm(
    modelo = m_surg,
    newdata = sim,
    niveles = niveles_surg
  )
  
  return(sim)
}


### 3.4.3 Prueba de la función
sim_datos <- simular_covariables(n = 500, seed = 123)

head(sim_datos)
str(sim_datos)
summary(sim_datos)


### 3.4.4 Comparación entre datos reales vs datos simulados en la prueba
summary(datos1$age_at_diagnosis)
summary(sim_datos$age_at_diagnosis)

summary(datos1$lymph_nodes_examined_positive)
summary(sim_datos$lymph_nodes_examined_positive)

summary(datos1$tumor_size)
summary(sim_datos$tumor_size)


table(datos1$her2_status)
table(sim_datos$her2_status)

table(datos1$neoplasm_histologic_grade)
table(sim_datos$neoplasm_histologic_grade)

table(datos1$chemotherapy)
table(sim_datos$chemotherapy)

table(datos1$type_of_breast_surgery)
table(sim_datos$type_of_breast_surgery)


