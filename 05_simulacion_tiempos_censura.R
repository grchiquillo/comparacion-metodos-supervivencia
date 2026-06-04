#########################################################################

# 5. Simulación tiempos de censura

#En el paso anterior los tiempos al evento se simularon a partir de 
#una distribución Weibull, cuya escala depende de las covariables a través 
#de un predictor de riesgo. Ahora se incorpora la censura mediante una variable 
#Bernoulli, la cual permite controlar la proporción de observaciones 
#censuradas en el conjunto de datos. 

#Por último se construyó el tiempo observado como el mínimo entre el tiempo 
#al evento y el tiempo de censura asignado.



## 5.1 Función para calibrar la escala de censura

#Esta función busca el valor de lambda_c que haga que la censura lograda se 
#acerque a 20%, 50%, 75% o 95% (los escenarios que quiero simular)

calibrar_lambda_censura <- function(data,
                                    tipo_relacion = c("lineal", "no_lineal"),
                                    gamma_evento = 1,
                                    beta0,
                                    gamma_cens = 1,
                                    censura_objetivo = 0.50,
                                    admin_time = NULL,
                                    n_piloto = 5000,
                                    seed = 123,
                                    lower_log = log(1e-4),
                                    upper_log = log(1e8)) {
  
  tipo_relacion <- match.arg(tipo_relacion)
  set.seed(seed)
  
  # Tomar muestra piloto de covariables
  if (nrow(data) > n_piloto) {
    idx <- sample(seq_len(nrow(data)), n_piloto)
    data_piloto <- data[idx, , drop = FALSE]
  } else {
    data_piloto <- data
  }
  
  # Simular una sola vez los tiempos verdaderos al evento
  sim_evento <- simular_tiempo_evento(
    data = data_piloto,
    tipo_relacion = tipo_relacion,
    gamma = gamma_evento,
    beta0 = beta0
  )
  
  T_true <- sim_evento$t_evento
  
  # Fijar aleatoriedad de censura para que la función objetivo sea estable
  set.seed(seed + 999)
  u_cens <- runif(length(T_true))
  
  f_obj <- function(log_lambda_c) {
    
    lambda_c <- exp(log_lambda_c)
    
    C <- qweibull(
      p = u_cens,
      shape = gamma_cens,
      scale = lambda_c
    )
    
    if (!is.null(admin_time)) {
      Y <- pmin(T_true, C, admin_time)
      delta <- as.integer(T_true <= C & T_true <= admin_time)
    } else {
      Y <- pmin(T_true, C)
      delta <- as.integer(T_true <= C)
    }
    
    censura_lograda <- mean(delta == 0)
    censura_lograda - censura_objetivo
  }
  
  f_lower <- f_obj(lower_log)
  f_upper <- f_obj(upper_log)
  
  # Caso normal: hay cambio de signo y se puede usar uniroot
  if (!is.na(f_lower) && !is.na(f_upper) && f_lower * f_upper <= 0) {
    
    raiz <- uniroot(
      f = f_obj,
      lower = lower_log,
      upper = upper_log,
      tol = 1e-8
    )
    
    lambda_c <- exp(raiz$root)
    return(lambda_c)
  }
  
  # Caso alternativo: no hay raíz exacta dentro del intervalo
  # Se busca el lambda que más se acerque a la censura objetivo
  
  grid_log <- seq(lower_log, upper_log, length.out = 200)
  valores <- sapply(grid_log, f_obj)
  
  idx_mejor <- which.min(abs(valores))
  lambda_c <- exp(grid_log[idx_mejor])
  
  censura_aprox <- valores[idx_mejor] + censura_objetivo
  
  warning(
    paste0(
      "No se encontró raíz exacta para censura_objetivo = ",
      censura_objetivo,
      ". Se usó lambda_c aproximado. Censura esperada aproximada = ",
      round(censura_aprox, 4)
    )
  )
  
  return(lambda_c)
}


