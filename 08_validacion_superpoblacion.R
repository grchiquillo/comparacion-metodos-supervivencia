###############################################################################

# 8. Validación en superpoblación de n=10.000

## 8.1. Preparación de rutas y objetos

ruta_super <- "prueba_superpoblacion_n10000"

dir.create(ruta_super, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "superpoblaciones"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "hiperparametros_srf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "hiperparametros_gbs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "validacion_en_linea"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "validacion_en_linea", "cox"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "validacion_en_linea", "srf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "validacion_en_linea", "gbs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "resultados_consolidados"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_super, "figuras"), showWarnings = FALSE, recursive = TRUE)

stopifnot(exists("tabla_24_escenarios"))
stopifnot(nrow(tabla_24_escenarios) == 24)

nombres_escenarios <- tabla_24_escenarios$escenario_nombre

length(nombres_escenarios)
head(nombres_escenarios)



## 8.2 Generación de superpoblaciones n=10.000 por escenario

generar_superpoblacion_escenario <- function(escenario_row,
                                             n_super = 10000,
                                             seed_base = 700000) {
  
  set.seed(seed_base + escenario_row$escenario_id)
  
  super_data <- simular_dataset_supervivencia(
    n = n_super,
    seed = seed_base + escenario_row$escenario_id,
    tipo_relacion = escenario_row$tipo_relacion,
    censura_objetivo = escenario_row$censura_objetivo,
    gamma_evento = 1,
    beta0 = beta0_calibrado,
    gamma_cens = 1,
    admin_time = 120
  )
  
  attr(super_data, "escenario_nombre") <- escenario_row$escenario_nombre
  attr(super_data, "escenario_id") <- escenario_row$escenario_id
  attr(super_data, "n_entrenamiento_original") <- escenario_row$n
  attr(super_data, "n_superpoblacion") <- n_super
  attr(super_data, "censura_objetivo") <- escenario_row$censura_objetivo
  attr(super_data, "tipo_relacion") <- escenario_row$tipo_relacion
  
  return(super_data)
}


tiempo_superpoblaciones <- system.time({
  
  for (i in seq_len(nrow(tabla_24_escenarios))) {
    
    esc_i <- tabla_24_escenarios[i, ]
    
    nombre_archivo <- file.path(
      ruta_super,
      "superpoblaciones",
      paste0(esc_i$escenario_nombre, "_superpob_n10000.rds")
    )
    
    cat("\nGenerando superpoblación:", esc_i$escenario_nombre, "\n")
    
    super_i <- generar_superpoblacion_escenario(
      escenario_row = esc_i,
      n_super = 10000,
      seed_base = 700000
    )
    
    saveRDS(super_i, nombre_archivo)
    
    cat(
      "Censura lograda:",
      round(mean(super_i$event == 0), 4),
      "| Eventos:",
      sum(super_i$event == 1),
      "| Censurados:",
      sum(super_i$event == 0),
      "\n"
    )
  }
})

tiempo_superpoblaciones


## 8.3 Cargar superpoblaciones y verificar correspondencia

superpoblaciones_10000 <- vector("list", nrow(tabla_24_escenarios))
names(superpoblaciones_10000) <- tabla_24_escenarios$escenario_nombre

for (i in seq_len(nrow(tabla_24_escenarios))) {
  
  nombre_esc <- tabla_24_escenarios$escenario_nombre[i]
  
  archivo_i <- file.path(
    ruta_super,
    "superpoblaciones",
    paste0(nombre_esc, "_superpob_n10000.rds")
  )
  
  superpoblaciones_10000[[nombre_esc]] <- readRDS(archivo_i)
}

# Verificación de correspondencia
verificacion_superpoblaciones <- purrr::map_dfr(
  names(superpoblaciones_10000),
  function(nombre_esc) {
    
    lista_esc <- escenarios_24_listas[[nombre_esc]]
    super_i <- superpoblaciones_10000[[nombre_esc]]
    
    tibble::tibble(
      escenario = nombre_esc,
      n_entrenamiento = attr(lista_esc, "n"),
      n_superpoblacion = nrow(super_i),
      censura_entrenamiento = attr(lista_esc, "censura_objetivo"),
      censura_super = attr(super_i, "censura_objetivo"),
      tipo_entrenamiento = attr(lista_esc, "tipo_relacion"),
      tipo_super = attr(super_i, "tipo_relacion"),
      censura_lograda_super = mean(super_i$event == 0),
      correspondencia = (
        attr(lista_esc, "censura_objetivo") == attr(super_i, "censura_objetivo") &
          attr(lista_esc, "tipo_relacion") == attr(super_i, "tipo_relacion") &
          nrow(super_i) == 10000
      )
    )
  }
)

verificacion_superpoblaciones

stopifnot(all(verificacion_superpoblaciones$correspondencia))


### 8.3.1 Funciones auxiliares generales

# Funciones base para ajuste, optimización y validación

formula_cox <- Surv(time, event) ~
  age_at_diagnosis +
  type_of_breast_surgery +
  chemotherapy +
  neoplasm_histologic_grade +
  her2_status +
  lymph_nodes_examined_positive +
  tumor_size

formula_srf <- formula_cox


alinear_factores_modelo <- function(train_data, test_data) {
  
  train_data <- as.data.frame(train_data)
  test_data  <- as.data.frame(test_data)
  
  vars_factor <- names(train_data)[sapply(train_data, is.factor)]
  
  for (v in vars_factor) {
    test_data[[v]] <- factor(
      as.character(test_data[[v]]),
      levels = levels(train_data[[v]])
    )
  }
  
  list(
    train = train_data,
    test = test_data
  )
}


crear_xy_xgb <- function(data) {
  
  data <- as.data.frame(data)
  
  terms_x <- stats::delete.response(stats::terms(formula_cox))
  
  x <- model.matrix(
    terms_x,
    data = data
  )
  
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  
  y <- ifelse(
    data$event == 1,
    data$time,
    -data$time
  )
  
  list(
    x = x,
    y = y
  )
}


alinear_columnas_xgb_por_nombres <- function(x_new, columnas_ref) {
  
  x_new <- as.matrix(x_new)
  
  cols_faltantes <- setdiff(columnas_ref, colnames(x_new))
  
  if (length(cols_faltantes) > 0) {
    for (cc in cols_faltantes) {
      x_new <- cbind(x_new, rep(0, nrow(x_new)))
      colnames(x_new)[ncol(x_new)] <- cc
    }
  }
  
  cols_sobrantes <- setdiff(colnames(x_new), columnas_ref)
  
  if (length(cols_sobrantes) > 0) {
    x_new <- x_new[, !colnames(x_new) %in% cols_sobrantes, drop = FALSE]
  }
  
  x_new <- x_new[, columnas_ref, drop = FALSE]
  
  x_new
}


ajustar_breslow_desde_lp <- function(data, lp) {
  
  data <- data %>%
    dplyr::mutate(lp = lp) %>%
    dplyr::arrange(time)
  
  tiempos_evento <- sort(unique(data$time[data$event == 1]))
  
  if (length(tiempos_evento) == 0) {
    return(NULL)
  }
  
  exp_lp <- exp(data$lp)
  
  incrementos <- sapply(tiempos_evento, function(tt) {
    
    d_t <- sum(data$time == tt & data$event == 1)
    riesgo_t <- sum(exp_lp[data$time >= tt])
    
    if (riesgo_t <= 0) {
      return(NA_real_)
    }
    
    d_t / riesgo_t
  })
  
  if (any(is.na(incrementos))) {
    return(NULL)
  }
  
  tibble::tibble(
    time = tiempos_evento,
    cumhaz = cumsum(incrementos)
  )
}


predecir_supervivencia_desde_lp <- function(lp,
                                            bh,
                                            times) {
  
  if (is.null(bh)) {
    return(NULL)
  }
  
  H0_t <- sapply(times, function(tt) {
    
    idx <- which(bh$time <= tt)
    
    if (length(idx) == 0) {
      return(0)
    } else {
      return(bh$cumhaz[max(idx)])
    }
  })
  
  surv_mat <- outer(
    exp(lp),
    H0_t,
    function(riesgo, h0) exp(-riesgo * h0)
  )
  
  surv_mat
}


### 8.3.2 Función homogénea para calcular IBS 


crear_G_censura_km <- function(data_eval) {
  
  fit_cens <- survival::survfit(
    survival::Surv(time, 1 - event) ~ 1,
    data = data_eval
  )
  
  tiempos_cens <- fit_cens$time
  surv_cens <- fit_cens$surv
  
  G_fun <- function(t) {
    
    idx <- findInterval(t, tiempos_cens)
    
    out <- ifelse(
      idx == 0,
      1,
      surv_cens[idx]
    )
    
    pmax(out, 1e-06)
  }
  
  return(G_fun)
}


brier_score_ipcw_graf <- function(data_eval,
                                  surv_prob_mat,
                                  times) {
  
  surv_prob_mat <- as.matrix(surv_prob_mat)
  
  idx_ok <- which(
    !is.na(data_eval$time) &
      !is.na(data_eval$event) &
      data_eval$time > 0
  )
  
  data_eval <- data_eval[idx_ok, ]
  
  if (nrow(surv_prob_mat) >= max(idx_ok)) {
    surv_prob_mat <- surv_prob_mat[idx_ok, , drop = FALSE]
  }
  
  if (nrow(surv_prob_mat) != nrow(data_eval) &&
      ncol(surv_prob_mat) == nrow(data_eval)) {
    surv_prob_mat <- t(surv_prob_mat)
  }
  
  if (nrow(surv_prob_mat) != nrow(data_eval)) {
    stop("La matriz de supervivencia no coincide con el número de sujetos.")
  }
  
  if (ncol(surv_prob_mat) != length(times)) {
    stop("La matriz de supervivencia no coincide con la grilla de tiempos.")
  }
  
  G_fun <- crear_G_censura_km(data_eval)
  
  y <- data_eval$time
  delta <- data_eval$event
  
  bs_t <- numeric(length(times))
  
  for (j in seq_along(times)) {
    
    tt <- times[j]
    
    y_t <- as.numeric(y > tt)
    
    G_t <- pmax(G_fun(tt), 1e-06)
    G_y <- pmax(G_fun(pmin(y, tt)), 1e-06)
    
    w <- ifelse(
      y <= tt & delta == 1,
      1 / G_y,
      ifelse(y > tt, 1 / G_t, 0)
    )
    
    bs_t[j] <- mean(
      w * (y_t - surv_prob_mat[, j])^2,
      na.rm = TRUE
    )
  }
  
  return(bs_t)
}


ibs_graf_ipcw <- function(data_eval,
                          surv_prob_mat,
                          times) {
  
  bs_t <- brier_score_ipcw_graf(
    data_eval = data_eval,
    surv_prob_mat = surv_prob_mat,
    times = times
  )
  
  idx_ok <- which(
    !is.na(bs_t) &
      is.finite(bs_t) &
      !is.na(times) &
      is.finite(times)
  )
  
  if (length(idx_ok) < 2) {
    return(list(
      ibs = NA_real_,
      times = times,
      brier = bs_t
    ))
  }
  
  times_ok <- times[idx_ok]
  bs_ok <- bs_t[idx_ok]
  
  ord <- order(times_ok)
  times_ok <- times_ok[ord]
  bs_ok <- bs_ok[ord]
  
  ibs <- sum(
    diff(times_ok) *
      (head(bs_ok, -1) + tail(bs_ok, -1)) / 2
  ) / (max(times_ok) - min(times_ok))
  
  return(list(
    ibs = ibs,
    times = times_ok,
    brier = bs_ok
  ))
}


### 8.3.3 Funciones de predicción de supervivencia

predecir_supervivencia_cox <- function(modelo_cox,
                                       newdata,
                                       times) {
  
  sf <- survival::survfit(
    modelo_cox,
    newdata = newdata
  )
  
  resumen <- summary(sf, times = times, extend = TRUE)
  surv_raw <- resumen$surv
  
  n_new <- nrow(newdata)
  n_times <- length(times)
  
  if (length(surv_raw) != n_new * n_times) {
    return(NULL)
  }
  
  surv_mat <- matrix(
    surv_raw,
    nrow = n_times,
    ncol = n_new
  )
  
  surv_mat <- t(surv_mat)
  
  return(surv_mat)
}


predecir_supervivencia_srf <- function(modelo_srf,
                                       newdata,
                                       times) {
  
  pred <- tryCatch(
    predict(modelo_srf, newdata = newdata),
    error = function(e) NULL
  )
  
  if (is.null(pred)) {
    return(NULL)
  }
  
  tiempos_modelo <- pred$time.interest
  surv_raw <- pred$survival
  
  if (is.null(surv_raw)) {
    return(NULL)
  }
  
  surv_raw <- as.matrix(surv_raw)
  
  if (nrow(surv_raw) != nrow(newdata) &&
      ncol(surv_raw) == nrow(newdata)) {
    surv_raw <- t(surv_raw)
  }
  
  if (nrow(surv_raw) != nrow(newdata)) {
    return(NULL)
  }
  
  surv_mat <- matrix(
    NA_real_,
    nrow = nrow(newdata),
    ncol = length(times)
  )
  
  for (j in seq_along(times)) {
    
    tt <- times[j]
    idx <- which(tiempos_modelo <= tt)
    
    if (length(idx) == 0) {
      surv_mat[, j] <- 1
    } else {
      surv_mat[, j] <- surv_raw[, max(idx)]
    }
  }
  
  risk_score <- pred$predicted
  
  if (is.null(risk_score) || length(risk_score) != nrow(newdata)) {
    risk_score <- -surv_mat[, ncol(surv_mat)]
  }
  
  return(list(
    surv = surv_mat,
    risk = risk_score
  ))
}


### 8.3.4 Funciones de optimización de hiperparámetros en GBS

# Funciones auxiliares para GBS
crear_xy_xgb <- function(data) {
  
  data <- droplevels(data)
  
  x <- model.matrix(
    formula_cox,
    data = data
  )
  
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  
  y <- ifelse(
    data$event == 1,
    data$time,
    -data$time
  )
  
  return(list(
    x = x,
    y = y
  ))
}


ajustar_breslow_desde_lp <- function(data, lp) {
  
  data <- data %>%
    dplyr::mutate(lp = lp) %>%
    dplyr::arrange(time)
  
  tiempos_evento <- sort(unique(data$time[data$event == 1]))
  
  if (length(tiempos_evento) == 0) {
    return(NULL)
  }
  
  exp_lp <- exp(data$lp)
  
  incrementos <- sapply(tiempos_evento, function(tt) {
    
    d_t <- sum(data$time == tt & data$event == 1)
    riesgo_t <- sum(exp_lp[data$time >= tt])
    
    if (riesgo_t <= 0) {
      return(NA_real_)
    }
    
    d_t / riesgo_t
  })
  
  if (any(is.na(incrementos))) {
    return(NULL)
  }
  
  tibble::tibble(
    time = tiempos_evento,
    cumhaz = cumsum(incrementos)
  )
}


predecir_supervivencia_desde_lp <- function(lp,
                                            bh,
                                            times) {
  
  H0_t <- sapply(times, function(tt) {
    
    idx <- which(bh$time <= tt)
    
    if (length(idx) == 0) {
      return(0)
    } else {
      return(bh$cumhaz[max(idx)])
    }
  })
  
  surv_mat <- outer(
    exp(lp),
    H0_t,
    function(riesgo, h0) exp(-riesgo * h0)
  )
  
  return(surv_mat)
}


## 8.4 Funciones de ajuste de los modelos de supervivencia

### 8.4.1 Ajuste CoxPH

ajustar_cox_modelo <- function(train_data) {
  
  train_data <- droplevels(train_data)
  
  modelo_cox <- tryCatch(
    suppressWarnings(
      coxph(
        formula_cox,
        data = train_data,
        ties = "breslow",
        x = TRUE,
        y = TRUE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(modelo_cox)) {
    return(list(
      modelo = NULL,
      estado_modelo = "fallo_ajuste"
    ))
  }
  
  coefs <- coef(modelo_cox)
  
  if (length(coefs) == 0 || any(is.na(coefs)) || any(!is.finite(coefs))) {
    return(list(
      modelo = modelo_cox,
      estado_modelo = "no_convergente"
    ))
  }
  
  return(list(
    modelo = modelo_cox,
    estado_modelo = "ok"
  ))
}


### 8.4.2 Ajuste SRF

ajustar_srf_modelo_fijo <- function(train_data, mejor_grid) {
  
  train_data <- droplevels(train_data)
  
  if (is.null(mejor_grid) || nrow(mejor_grid) == 0) {
    return(list(
      modelo = NULL,
      estado_modelo = "sin_hiperparametros"
    ))
  }
  
  modelo_srf <- tryCatch(
    suppressWarnings(
      randomForestSRC::rfsrc(
        formula_srf,
        data = train_data,
        ntree = mejor_grid$ntree[1],
        mtry = mejor_grid$mtry[1],
        nodesize = mejor_grid$nodesize[1],
        splitrule = as.character(mejor_grid$splitrule[1]),
        importance = "none",
        forest = TRUE,
        save.memory = TRUE,
        perf.type = "none",
        statistics = FALSE,
        proximity = FALSE,
        distance = FALSE,
        forest.wt = FALSE,
        membership = FALSE,
        var.used = FALSE,
        split.depth = FALSE
      )
    ),
    error = function(e) {
      message("Error en ajuste SRF: ", e$message)
      NULL
    }
  )
  
  if (is.null(modelo_srf)) {
    return(list(
      modelo = NULL,
      estado_modelo = "fallo_ajuste"
    ))
  }
  
  return(list(
    modelo = modelo_srf,
    estado_modelo = "ok",
    mejor_grid = mejor_grid
  ))
}



### 8.4.3 Ajuste GBS

alinear_columnas_xgb_por_nombres <- function(x_new, columnas_ref) {
  
  x_new <- as.matrix(x_new)
  
  cols_faltantes <- setdiff(columnas_ref, colnames(x_new))
  
  if (length(cols_faltantes) > 0) {
    for (cc in cols_faltantes) {
      x_new <- cbind(x_new, rep(0, nrow(x_new)))
      colnames(x_new)[ncol(x_new)] <- cc
    }
  }
  
  cols_sobrantes <- setdiff(colnames(x_new), columnas_ref)
  
  if (length(cols_sobrantes) > 0) {
    x_new <- x_new[, !colnames(x_new) %in% cols_sobrantes, drop = FALSE]
  }
  
  x_new <- x_new[, columnas_ref, drop = FALSE]
  
  return(x_new)
}


ajustar_gbs_modelo_fijo <- function(train_data, mejor_grid) {
  
  train_data <- droplevels(train_data)
  
  if (is.null(mejor_grid) || nrow(mejor_grid) == 0) {
    return(list(
      modelo = NULL,
      estado_modelo = "sin_hiperparametros"
    ))
  }
  
  xy_train <- crear_xy_xgb(train_data)
  x_train <- xy_train$x
  
  dtrain <- xgb.DMatrix(
    data = x_train,
    label = xy_train$y
  )
  
  nrounds_usar <- if ("best_nrounds" %in% names(mejor_grid)) {
    mejor_grid$best_nrounds[1]
  } else {
    mejor_grid$nrounds[1]
  }
  
  if (is.na(nrounds_usar) || !is.finite(nrounds_usar) || nrounds_usar < 1) {
    nrounds_usar <- mejor_grid$nrounds[1]
  }
  
  modelo_gbs <- tryCatch(
    xgb.train(
      params = list(
        objective = "survival:cox",
        eval_metric = "cox-nloglik",
        eta = mejor_grid$eta[1],
        max_depth = mejor_grid$max_depth[1],
        subsample = mejor_grid$subsample[1],
        colsample_bytree = mejor_grid$colsample_bytree[1]
      ),
      data = dtrain,
      nrounds = nrounds_usar,
      verbose = 0
    ),
    error = function(e) NULL
  )
  
  if (is.null(modelo_gbs)) {
    return(list(
      modelo = NULL,
      estado_modelo = "fallo_ajuste"
    ))
  }
  
  lp_train <- predict(modelo_gbs, dtrain)
  
  bh <- tryCatch(
    ajustar_breslow_desde_lp(train_data, lp_train),
    error = function(e) NULL
  )
  
  if (is.null(bh)) {
    return(list(
      modelo = modelo_gbs,
      estado_modelo = "fallo_breslow"
    ))
  }
  
  return(list(
    modelo = modelo_gbs,
    estado_modelo = "ok",
    mejor_grid = mejor_grid,
    bh = bh,
    columnas_x = colnames(x_train),
    nrounds_usado = nrounds_usar
  ))
}



## 8.5 Optimización de hiperparámetros para SRF y GBS

### 8.5.0 Definición grillas de hiperparámetros

dir.create(file.path(ruta_super, "grillas_hiperparametros"), showWarnings = FALSE)
dir.create(file.path(ruta_super, "hiperparametros_srf"), showWarnings = FALSE)
dir.create(file.path(ruta_super, "hiperparametros_gbs"), showWarnings = FALSE)

grid_srf <- expand.grid(
  mtry = c(2, 3, 4),
  nodesize = c(3, 5, 10),
  ntree = c(250, 500),
  splitrule = c("logrank"),
  stringsAsFactors = FALSE
)

grid_gbs <- expand.grid(
  eta = c(0.01, 0.1),
  max_depth = c(3, 7),
  nrounds = c(100, 500),
  subsample = c(0.7, 1),
  colsample_bytree = c(0.7, 1),
  stringsAsFactors = FALSE
)

grid_srf
nrow(grid_srf)

grid_gbs
nrow(grid_gbs)

write.csv(
  grid_srf,
  file.path(ruta_super, "grillas_hiperparametros", "grid_srf_superpoblacion.csv"),
  row.names = FALSE
)

write.csv(
  grid_gbs,
  file.path(ruta_super, "grillas_hiperparametros", "grid_gbs_superpoblacion.csv"),
  row.names = FALSE
)


### 8.5.1 Función auxiliar para validación cruzada interna en la optimización de hiperparámetros

# Funciones auxiliares para validación cruzada interna

crear_folds_estratificados <- function(event, k = 5, seed = 123) {
  
  set.seed(seed)
  
  idx_evento <- sample(which(event == 1))
  idx_cens   <- sample(which(event == 0))
  
  k_efectivo <- min(k, length(idx_evento), length(idx_cens))
  
  if (k_efectivo < 2) {
    stop("No hay suficientes eventos y censurados para construir al menos 2 folds.")
  }
  
  split_evento <- split(
    idx_evento,
    rep(seq_len(k_efectivo), length.out = length(idx_evento))
  )
  
  split_cens <- split(
    idx_cens,
    rep(seq_len(k_efectivo), length.out = length(idx_cens))
  )
  
  folds <- vector("list", k_efectivo)
  
  for (i in seq_len(k_efectivo)) {
    folds[[i]] <- sort(c(split_evento[[i]], split_cens[[i]]))
  }
  
  folds
}


crear_time_grid_cv <- function(data,
                               n_puntos = 30,
                               q_min = 0.05,
                               q_max = 0.80) {
  
  tiempos_evento <- data$time[data$event == 1]
  
  if (length(tiempos_evento) >= 10) {
    t_min <- as.numeric(stats::quantile(tiempos_evento, q_min, na.rm = TRUE))
    t_max <- as.numeric(stats::quantile(tiempos_evento, q_max, na.rm = TRUE))
  } else {
    t_min <- max(0.01, min(data$time, na.rm = TRUE))
    t_max <- as.numeric(stats::quantile(data$time, q_max, na.rm = TRUE))
  }
  
  if (!is.finite(t_min) || !is.finite(t_max) || t_max <= t_min) {
    t_min <- max(0.01, min(data$time, na.rm = TRUE))
    t_max <- as.numeric(stats::quantile(data$time, 0.80, na.rm = TRUE))
  }
  
  seq(t_min, t_max, length.out = n_puntos)
}


crear_funcion_censura_cv <- function(data_ref) {
  
  fit_cens <- survival::survfit(
    survival::Surv(time, 1 - event) ~ 1,
    data = data_ref
  )
  
  tiempos <- fit_cens$time
  surv    <- fit_cens$surv
  
  G_fun <- function(t) {
    
    sapply(t, function(tt) {
      
      idx <- which(tiempos <= tt)
      
      if (length(idx) == 0) {
        return(1)
      } else {
        return(max(surv[max(idx)], 1e-06))
      }
    })
  }
  
  G_fun
}


ibs_ipcw_cv <- function(data_ref,
                        data_eval,
                        surv_prob_mat,
                        times) {
  
  if (is.null(surv_prob_mat) || length(times) < 2) {
    return(NA_real_)
  }
  
  surv_prob_mat <- as.matrix(surv_prob_mat)
  
  if (nrow(surv_prob_mat) != nrow(data_eval) &&
      ncol(surv_prob_mat) == nrow(data_eval)) {
    surv_prob_mat <- t(surv_prob_mat)
  }
  
  if (nrow(surv_prob_mat) != nrow(data_eval) ||
      ncol(surv_prob_mat) != length(times)) {
    return(NA_real_)
  }
  
  G_fun <- crear_funcion_censura_cv(data_ref)
  
  y <- data_eval$time
  delta <- data_eval$event
  
  G_y <- G_fun(pmax(y - 1e-08, 0))
  G_y[G_y < 1e-06] <- 1e-06
  
  bs_t <- rep(NA_real_, length(times))
  
  for (j in seq_along(times)) {
    
    tt <- times[j]
    G_t <- max(G_fun(tt), 1e-06)
    
    term_evento <- as.numeric(y <= tt & delta == 1) *
      ((0 - surv_prob_mat[, j])^2) / G_y
    
    term_superv <- as.numeric(y > tt) *
      ((1 - surv_prob_mat[, j])^2) / G_t
    
    bs_t[j] <- mean(term_evento + term_superv, na.rm = TRUE)
  }
  
  idx_ok <- which(!is.na(bs_t) & is.finite(bs_t))
  
  if (length(idx_ok) < 2) {
    return(NA_real_)
  }
  
  times_ok <- times[idx_ok]
  bs_ok <- bs_t[idx_ok]
  
  sum(
    diff(times_ok) *
      (head(bs_ok, -1) + tail(bs_ok, -1)) / 2
  ) / (max(times_ok) - min(times_ok))
}


### 8.5.2 Funciones de predicción para SRF en validación cruzada
# Predicción interna para SRF

predecir_supervivencia_srf_cv <- function(modelo_srf,
                                          newdata,
                                          times) {
  
  pred <- tryCatch(
    predict(modelo_srf, newdata = newdata),
    error = function(e) NULL
  )
  
  if (is.null(pred)) {
    return(NULL)
  }
  
  tiempos_modelo <- pred$time.interest
  surv_raw <- pred$survival
  
  if (is.null(surv_raw)) {
    return(NULL)
  }
  
  surv_raw <- as.matrix(surv_raw)
  
  if (nrow(surv_raw) != nrow(newdata) &&
      ncol(surv_raw) == nrow(newdata)) {
    surv_raw <- t(surv_raw)
  }
  
  if (nrow(surv_raw) != nrow(newdata)) {
    return(NULL)
  }
  
  surv_mat <- matrix(
    NA_real_,
    nrow = nrow(newdata),
    ncol = length(times)
  )
  
  for (j in seq_along(times)) {
    
    idx <- which(tiempos_modelo <= times[j])
    
    if (length(idx) == 0) {
      surv_mat[, j] <- 1
    } else {
      surv_mat[, j] <- surv_raw[, max(idx)]
    }
  }
  
  surv_mat
}




### 8.5.3 Evaluar grilla SRF por validación cruzada

evaluar_grid_srf_cv <- function(train_data,
                                grid_srf,
                                k = 5,
                                seed = 123) {
  
  train_data <- as.data.frame(train_data)
  
  folds <- crear_folds_estratificados(
    event = train_data$event,
    k = k,
    seed = seed
  )
  
  tiempos_eval <- crear_time_grid_cv(
    data = train_data,
    n_puntos = 30
  )
  
  resultados_grid <- vector("list", nrow(grid_srf))
  
  for (g in seq_len(nrow(grid_srf))) {
    
    pars <- grid_srf[g, ]
    ibs_folds <- rep(NA_real_, length(folds))
    
    for (f in seq_along(folds)) {
      
      idx_val <- folds[[f]]
      
      dat_tr  <- train_data[-idx_val, , drop = FALSE]
      dat_val <- train_data[idx_val, , drop = FALSE]
      
      modelo_srf <- tryCatch(
        suppressWarnings(
          randomForestSRC::rfsrc(
            formula_srf,
            data = dat_tr,
            ntree = as.integer(pars$ntree),
            mtry = as.integer(pars$mtry),
            nodesize = as.integer(pars$nodesize),
            splitrule = as.character(pars$splitrule),
            importance = "none",
            forest = TRUE,
            save.memory = TRUE,
            perf.type = "none",
            statistics = FALSE,
            proximity = FALSE,
            distance = FALSE,
            forest.wt = FALSE,
            membership = FALSE,
            var.used = FALSE,
            split.depth = FALSE
          )
        ),
        error = function(e) NULL
      )
      
      if (is.null(modelo_srf)) next
      
      surv_val <- predecir_supervivencia_srf_cv(
        modelo_srf = modelo_srf,
        newdata = dat_val,
        times = tiempos_eval
      )
      
      if (is.null(surv_val)) next
      
      ibs_folds[f] <- ibs_ipcw_cv(
        data_ref = dat_tr,
        data_eval = dat_val,
        surv_prob_mat = surv_val,
        times = tiempos_eval
      )
    }
    
    resultados_grid[[g]] <- tibble::tibble(
      mtry = pars$mtry,
      nodesize = pars$nodesize,
      ntree = pars$ntree,
      splitrule = as.character(pars$splitrule),
      ibs_cv = mean(ibs_folds, na.rm = TRUE),
      sd_ibs_cv = sd(ibs_folds, na.rm = TRUE),
      folds_validos = sum(!is.na(ibs_folds))
    )
  }
  
  dplyr::bind_rows(resultados_grid) %>%
    dplyr::filter(
      folds_validos > 0,
      !is.na(ibs_cv),
      is.finite(ibs_cv)
    ) %>%
    dplyr::arrange(ibs_cv)
}



### 8.5.4 Evaluar grilla GBS por validación cruzada

evaluar_grid_gbs_cv <- function(train_data,
                                grid_gbs,
                                k = 5,
                                seed = 123,
                                early_stopping_rounds = 50) {
  
  train_data <- as.data.frame(train_data)
  
  folds <- crear_folds_estratificados(
    event = train_data$event,
    k = k,
    seed = seed
  )
  
  tiempos_eval <- crear_time_grid_cv(
    data = train_data,
    n_puntos = 30
  )
  
  resultados_grid <- vector("list", nrow(grid_gbs))
  
  for (g in seq_len(nrow(grid_gbs))) {
    
    pars <- grid_gbs[g, ]
    
    ibs_folds <- rep(NA_real_, length(folds))
    best_iter_folds <- rep(NA_real_, length(folds))
    
    for (f in seq_along(folds)) {
      
      idx_val <- folds[[f]]
      
      dat_tr  <- train_data[-idx_val, , drop = FALSE]
      dat_val <- train_data[idx_val, , drop = FALSE]
      
      xy_tr <- crear_xy_xgb(dat_tr)
      xy_val <- crear_xy_xgb(dat_val)
      
      x_val <- alinear_columnas_xgb_por_nombres(
        x_new = xy_val$x,
        columnas_ref = colnames(xy_tr$x)
      )
      
      dtrain <- xgboost::xgb.DMatrix(
        data = xy_tr$x,
        label = xy_tr$y
      )
      
      dval <- xgboost::xgb.DMatrix(
        data = x_val,
        label = xy_val$y
      )
      
      modelo_gbs <- tryCatch(
        xgboost::xgb.train(
          params = list(
            objective = "survival:cox",
            eval_metric = "cox-nloglik",
            eta = pars$eta,
            max_depth = pars$max_depth,
            subsample = pars$subsample,
            colsample_bytree = pars$colsample_bytree
          ),
          data = dtrain,
          nrounds = pars$nrounds,
          watchlist = list(train = dtrain, valid = dval),
          early_stopping_rounds = early_stopping_rounds,
          verbose = 0
        ),
        error = function(e) NULL
      )
      
      if (is.null(modelo_gbs)) next
      
      lp_tr <- tryCatch(predict(modelo_gbs, dtrain), error = function(e) NULL)
      lp_val <- tryCatch(predict(modelo_gbs, dval), error = function(e) NULL)
      
      if (is.null(lp_tr) || is.null(lp_val)) next
      
      bh <- ajustar_breslow_desde_lp(
        data = dat_tr,
        lp = lp_tr
      )
      
      if (is.null(bh)) next
      
      surv_val <- predecir_supervivencia_desde_lp(
        lp = lp_val,
        bh = bh,
        times = tiempos_eval
      )
      
      if (is.null(surv_val)) next
      
      ibs_folds[f] <- ibs_ipcw_cv(
        data_ref = dat_tr,
        data_eval = dat_val,
        surv_prob_mat = surv_val,
        times = tiempos_eval
      )
      
      best_iter_folds[f] <- modelo_gbs$best_iteration
    }
    
    best_nrounds <- round(mean(best_iter_folds, na.rm = TRUE))
    
    if (is.na(best_nrounds) || !is.finite(best_nrounds) || best_nrounds < 1) {
      best_nrounds <- pars$nrounds
    }
    
    resultados_grid[[g]] <- tibble::tibble(
      eta = pars$eta,
      max_depth = pars$max_depth,
      nrounds = pars$nrounds,
      subsample = pars$subsample,
      colsample_bytree = pars$colsample_bytree,
      best_nrounds = best_nrounds,
      ibs_cv = mean(ibs_folds, na.rm = TRUE),
      sd_ibs_cv = sd(ibs_folds, na.rm = TRUE),
      folds_validos = sum(!is.na(ibs_folds))
    )
  }
  
  dplyr::bind_rows(resultados_grid) %>%
    dplyr::filter(
      folds_validos > 0,
      !is.na(ibs_cv),
      is.finite(ibs_cv)
    ) %>%
    dplyr::arrange(ibs_cv)
}



### 8.5.5 Selección de repetición piloto viable y de hiperparámetros


diagnosticar_repeticiones_escenario <- function(lista_escenario) {
  
  purrr::map_dfr(
    seq_along(lista_escenario),
    function(r) {
      
      dat <- lista_escenario[[r]]
      
      eventos <- sum(dat$event == 1, na.rm = TRUE)
      censurados <- sum(dat$event == 0, na.rm = TRUE)
      
      tibble::tibble(
        rep = r,
        eventos = eventos,
        censurados = censurados,
        minimo_clase = min(eventos, censurados)
      )
    }
  ) %>%
    dplyr::arrange(dplyr::desc(minimo_clase), rep)
}


seleccionar_repeticion_piloto_viable <- function(lista_escenario,
                                                 rep_piloto = 1,
                                                 k = 5,
                                                 min_clase = 2) {
  
  diagnostico <- diagnosticar_repeticiones_escenario(lista_escenario)
  
  candidatos <- diagnostico %>%
    dplyr::filter(eventos >= min_clase, censurados >= min_clase) %>%
    dplyr::mutate(
      prioridad = ifelse(rep == rep_piloto, 0, 1)
    ) %>%
    dplyr::arrange(prioridad, dplyr::desc(minimo_clase), rep)
  
  if (nrow(candidatos) == 0) {
    stop("No existe ninguna repetición con suficientes eventos y censurados.")
  }
  
  candidato <- candidatos[1, ]
  
  list(
    rep_piloto_usada = candidato$rep,
    k_usado = min(k, candidato$minimo_clase),
    eventos_piloto = candidato$eventos,
    censurados_piloto = candidato$censurados,
    diagnostico = diagnostico
  )
}


seleccionar_hiperparametros_srf_super <- function(lista_escenario,
                                                  grid_srf,
                                                  rep_piloto = 1,
                                                  k = 5,
                                                  seed = 123) {
  
  piloto <- seleccionar_repeticion_piloto_viable(
    lista_escenario = lista_escenario,
    rep_piloto = rep_piloto,
    k = k,
    min_clase = 2
  )
  
  train_data <- lista_escenario[[piloto$rep_piloto_usada]]
  
  grid_cv <- evaluar_grid_srf_cv(
    train_data = train_data,
    grid_srf = grid_srf,
    k = piloto$k_usado,
    seed = seed
  )
  
  if (is.null(grid_cv) || nrow(grid_cv) == 0) {
    stop("La optimización SRF no produjo combinaciones válidas.")
  }
  
  mejor_grid <- grid_cv[1, ] %>%
    dplyr::mutate(
      rep_piloto_usada = piloto$rep_piloto_usada,
      k_cv_usado = piloto$k_usado,
      eventos_piloto = piloto$eventos_piloto,
      censurados_piloto = piloto$censurados_piloto
    )
  
  list(
    mejor_grid = mejor_grid,
    grid_cv = grid_cv,
    diagnostico_repeticiones = piloto$diagnostico
  )
}


seleccionar_hiperparametros_gbs_super <- function(lista_escenario,
                                                  grid_gbs,
                                                  rep_piloto = 1,
                                                  k = 5,
                                                  seed = 123) {
  
  piloto <- seleccionar_repeticion_piloto_viable(
    lista_escenario = lista_escenario,
    rep_piloto = rep_piloto,
    k = k,
    min_clase = 2
  )
  
  train_data <- lista_escenario[[piloto$rep_piloto_usada]]
  
  grid_cv <- evaluar_grid_gbs_cv(
    train_data = train_data,
    grid_gbs = grid_gbs,
    k = piloto$k_usado,
    seed = seed
  )
  
  if (is.null(grid_cv) || nrow(grid_cv) == 0) {
    stop("La optimización GBS no produjo combinaciones válidas.")
  }
  
  mejor_grid <- grid_cv[1, ] %>%
    dplyr::mutate(
      rep_piloto_usada = piloto$rep_piloto_usada,
      k_cv_usado = piloto$k_usado,
      eventos_piloto = piloto$eventos_piloto,
      censurados_piloto = piloto$censurados_piloto
    )
  
  list(
    mejor_grid = mejor_grid,
    grid_cv = grid_cv,
    diagnostico_repeticiones = piloto$diagnostico
  )
}




### 8.5.6 Seleccionar hiperparámetros SRF y GBS por escenario

seleccionar_hiperparametros_srf_super <- function(lista_escenario,
                                                  grid_srf,
                                                  rep_piloto = 1,
                                                  k = 5,
                                                  seed = 123) {
  
  piloto <- seleccionar_repeticion_piloto_viable(
    lista_escenario = lista_escenario,
    rep_piloto = rep_piloto,
    k = k,
    min_clase = 2
  )
  
  if (piloto$rep_piloto_usada != rep_piloto) {
    cat(
      "La repetición piloto", rep_piloto,
      "no fue viable. Se usó rep_piloto =",
      piloto$rep_piloto_usada,
      "| eventos =", piloto$eventos_piloto,
      "| censurados =", piloto$censurados_piloto,
      "| k usado =", piloto$k_usado,
      "\n"
    )
  }
  
  train_data <- lista_escenario[[piloto$rep_piloto_usada]]
  train_data <- droplevels(train_data)
  
  grid_cv <- evaluar_grid_srf_cv(
    train_data = train_data,
    grid_srf = grid_srf,
    k = piloto$k_usado,
    seed = seed
  )
  
  if (is.null(grid_cv) || nrow(grid_cv) == 0) {
    stop("La optimización SRF no produjo combinaciones válidas.")
  }
  
  mejor_grid <- grid_cv[1, ] %>%
    dplyr::mutate(
      rep_piloto_usada = piloto$rep_piloto_usada,
      k_cv_usado = piloto$k_usado,
      eventos_piloto = piloto$eventos_piloto,
      censurados_piloto = piloto$censurados_piloto
    )
  
  list(
    mejor_grid = mejor_grid,
    grid_cv = grid_cv,
    diagnostico_repeticiones = piloto$diagnostico
  )
}


seleccionar_hiperparametros_gbs_super <- function(lista_escenario,
                                                  grid_gbs,
                                                  rep_piloto = 1,
                                                  k = 5,
                                                  seed = 123,
                                                  reps_candidatas = NULL,
                                                  k_candidatos = c(5, 3, 2),
                                                  max_reps_candidatas = 20) {
  
  diagnostico_reps <- diagnosticar_repeticiones_escenario(
    lista_escenario = lista_escenario
  )
  
  if (is.null(reps_candidatas)) {
    
    reps_viables <- diagnostico_reps %>%
      dplyr::filter(
        eventos >= 2,
        censurados >= 2
      ) %>%
      dplyr::arrange(
        dplyr::desc(minimo_clase),
        rep
      ) %>%
      dplyr::pull(rep)
    
    reps_candidatas <- unique(
      c(
        rep_piloto,
        reps_viables
      )
    )
    
    reps_candidatas <- reps_candidatas[
      reps_candidatas %in% seq_along(lista_escenario)
    ]
    
    reps_candidatas <- head(
      reps_candidatas,
      max_reps_candidatas
    )
  }
  
  k_candidatos <- unique(
    c(
      k,
      k_candidatos
    )
  )
  
  k_candidatos <- k_candidatos[
    !is.na(k_candidatos) &
      is.finite(k_candidatos) &
      k_candidatos >= 2
  ]
  
  diagnostico_intentos <- list()
  contador_intentos <- 0
  
  for (k_i in k_candidatos) {
    
    for (rep_i in reps_candidatas) {
      
      cat(
        "\nIntento GBS | rep_piloto =",
        rep_i,
        "| k =",
        k_i,
        "\n"
      )
      
      piloto <- tryCatch(
        seleccionar_repeticion_piloto_viable(
          lista_escenario = lista_escenario,
          rep_piloto = rep_i,
          k = k_i,
          min_clase = 2
        ),
        error = function(e) {
          list(
            error = TRUE,
            mensaje = conditionMessage(e)
          )
        }
      )
      
      if (!is.null(piloto$error) && piloto$error) {
        
        contador_intentos <- contador_intentos + 1
        
        diagnostico_intentos[[contador_intentos]] <- tibble::tibble(
          rep_piloto_solicitada = rep_i,
          rep_piloto_usada = NA_integer_,
          k_solicitado = k_i,
          k_cv_usado = NA_integer_,
          eventos_piloto = NA_integer_,
          censurados_piloto = NA_integer_,
          estado = "fallo_piloto",
          mensaje = piloto$mensaje
        )
        
        next
      }
      
      if (piloto$rep_piloto_usada != rep_i) {
        cat(
          "La repetición piloto", rep_i,
          "no fue viable. Se usó rep_piloto =",
          piloto$rep_piloto_usada,
          "| eventos =", piloto$eventos_piloto,
          "| censurados =", piloto$censurados_piloto,
          "| k usado =", piloto$k_usado,
          "\n"
        )
      }
      
      train_data <- lista_escenario[[piloto$rep_piloto_usada]]
      train_data <- droplevels(train_data)
      
      grid_cv <- tryCatch(
        evaluar_grid_gbs_cv(
          train_data = train_data,
          grid_gbs = grid_gbs,
          k = piloto$k_usado,
          seed = seed + rep_i + k_i
        ),
        error = function(e) {
          tibble::tibble(
            error = TRUE,
            mensaje = conditionMessage(e)
          )
        }
      )
      
      if (
        is.data.frame(grid_cv) &&
        "error" %in% names(grid_cv) &&
        isTRUE(grid_cv$error[1])
      ) {
        
        contador_intentos <- contador_intentos + 1
        
        diagnostico_intentos[[contador_intentos]] <- tibble::tibble(
          rep_piloto_solicitada = rep_i,
          rep_piloto_usada = piloto$rep_piloto_usada,
          k_solicitado = k_i,
          k_cv_usado = piloto$k_usado,
          eventos_piloto = piloto$eventos_piloto,
          censurados_piloto = piloto$censurados_piloto,
          estado = "fallo_grid",
          mensaje = grid_cv$mensaje[1]
        )
        
        next
      }
      
      if (is.null(grid_cv) || nrow(grid_cv) == 0) {
        
        contador_intentos <- contador_intentos + 1
        
        diagnostico_intentos[[contador_intentos]] <- tibble::tibble(
          rep_piloto_solicitada = rep_i,
          rep_piloto_usada = piloto$rep_piloto_usada,
          k_solicitado = k_i,
          k_cv_usado = piloto$k_usado,
          eventos_piloto = piloto$eventos_piloto,
          censurados_piloto = piloto$censurados_piloto,
          estado = "fallo_sin_combinaciones_validas",
          mensaje = "La optimización GBS no produjo combinaciones válidas."
        )
        
        next
      }
      
      mejor_grid <- grid_cv[1, ] %>%
        dplyr::mutate(
          rep_piloto_solicitada = rep_i,
          rep_piloto_usada = piloto$rep_piloto_usada,
          k_solicitado = k_i,
          k_cv_usado = piloto$k_usado,
          eventos_piloto = piloto$eventos_piloto,
          censurados_piloto = piloto$censurados_piloto
        )
      
      contador_intentos <- contador_intentos + 1
      
      diagnostico_intentos[[contador_intentos]] <- tibble::tibble(
        rep_piloto_solicitada = rep_i,
        rep_piloto_usada = piloto$rep_piloto_usada,
        k_solicitado = k_i,
        k_cv_usado = piloto$k_usado,
        eventos_piloto = piloto$eventos_piloto,
        censurados_piloto = piloto$censurados_piloto,
        estado = "ok",
        mensaje = NA_character_
      )
      
      return(
        list(
          mejor_grid = mejor_grid,
          grid_cv = grid_cv,
          diagnostico_repeticiones = diagnostico_reps,
          diagnostico_intentos = dplyr::bind_rows(diagnostico_intentos)
        )
      )
    }
  }
  
  diagnostico_intentos_df <- dplyr::bind_rows(diagnostico_intentos)
  
  stop(
    paste0(
      "La optimización GBS no produjo combinaciones válidas después de ",
      nrow(diagnostico_intentos_df),
      " intentos."
    )
  )
}



### 8.5.7 Ejecutar optimización SRF para los 24 escenarios

registro_hp_srf_super <- tibble::tibble(
  escenario = names(escenarios_24_listas),
  estado = NA_character_,
  rep_piloto_usada = NA_integer_,
  k_cv_usado = NA_integer_,
  eventos_piloto = NA_integer_,
  censurados_piloto = NA_integer_,
  mensaje_error = NA_character_
)

tiempo_hp_srf_super <- system.time({
  
  for (i in seq_along(escenarios_24_listas)) {
    
    nombre_esc <- names(escenarios_24_listas)[i]
    lista_esc <- escenarios_24_listas[[i]]
    
    cat("\n====================================\n")
    cat("Optimizando SRF:", nombre_esc, "\n")
    cat("Escenario", i, "de", length(escenarios_24_listas), "\n")
    cat("====================================\n")
    
    seleccion_hp <- tryCatch(
      seleccionar_hiperparametros_srf_super(
        lista_escenario = lista_esc,
        grid_srf = grid_srf,
        rep_piloto = 1,
        k = 5,
        seed = 123
      ),
      error = function(e) {
        list(error = TRUE, mensaje = conditionMessage(e))
      }
    )
    
    if (!is.null(seleccion_hp$error) && seleccion_hp$error) {
      
      registro_hp_srf_super$estado[i] <- "fallo"
      registro_hp_srf_super$mensaje_error[i] <- seleccion_hp$mensaje
      
      write.csv(
        registro_hp_srf_super,
        file.path(ruta_super, "hiperparametros_srf", "registro_hp_srf_super.csv"),
        row.names = FALSE
      )
      
      next
    }
    
    mejor_grid <- seleccion_hp$mejor_grid
    grid_cv <- seleccion_hp$grid_cv
    
    saveRDS(
      mejor_grid,
      file.path(
        ruta_super,
        "hiperparametros_srf",
        paste0(nombre_esc, "_mejor_grid_srf.rds")
      )
    )
    
    write.csv(
      mejor_grid,
      file.path(
        ruta_super,
        "hiperparametros_srf",
        paste0(nombre_esc, "_mejor_grid_srf.csv")
      ),
      row.names = FALSE
    )
    
    write.csv(
      grid_cv,
      file.path(
        ruta_super,
        "hiperparametros_srf",
        paste0(nombre_esc, "_grid_cv_srf.csv")
      ),
      row.names = FALSE
    )
    
    write.csv(
      seleccion_hp$diagnostico_repeticiones,
      file.path(
        ruta_super,
        "hiperparametros_srf",
        paste0(nombre_esc, "_diagnostico_reps_srf.csv")
      ),
      row.names = FALSE
    )
    
    registro_hp_srf_super$estado[i] <- "ok"
    registro_hp_srf_super$rep_piloto_usada[i] <- mejor_grid$rep_piloto_usada[1]
    registro_hp_srf_super$k_cv_usado[i] <- mejor_grid$k_cv_usado[1]
    registro_hp_srf_super$eventos_piloto[i] <- mejor_grid$eventos_piloto[1]
    registro_hp_srf_super$censurados_piloto[i] <- mejor_grid$censurados_piloto[1]
    
    write.csv(
      registro_hp_srf_super,
      file.path(ruta_super, "hiperparametros_srf", "registro_hp_srf_super.csv"),
      row.names = FALSE
    )
  }
})

tiempo_hp_srf_super

registro_hp_srf_super
table(registro_hp_srf_super$estado, useNA = "ifany")



### 8.5.8 Ejecutar optimización GBS para los 24 escenarios

dir.create(
  file.path(ruta_super, "hiperparametros_gbs"),
  showWarnings = FALSE,
  recursive = TRUE
)

sobrescribir_hp_gbs <- TRUE

registro_hp_gbs_super <- tibble::tibble(
  escenario = tabla_24_escenarios$escenario_nombre,
  estado = NA_character_,
  rep_piloto_solicitada = NA_integer_,
  rep_piloto_usada = NA_integer_,
  k_solicitado = NA_integer_,
  k_cv_usado = NA_integer_,
  eventos_piloto = NA_integer_,
  censurados_piloto = NA_integer_,
  mensaje_error = NA_character_
)

tiempo_hp_gbs_super <- system.time({
  
  for (i in seq_len(nrow(tabla_24_escenarios))) {
    
    nombre_esc <- tabla_24_escenarios$escenario_nombre[i]
    lista_esc <- escenarios_24_listas[[nombre_esc]]
    
    cat("\n====================================\n")
    cat("Optimizando GBS:", nombre_esc, "\n")
    cat("Escenario", i, "de", nrow(tabla_24_escenarios), "\n")
    cat("====================================\n")
    
    archivo_mejor_rds <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_mejor_grid_gbs.rds")
    )
    
    archivo_mejor_csv <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_mejor_grid_gbs.csv")
    )
    
    archivo_grid_csv <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_grid_cv_gbs.csv")
    )
    
    archivo_diag_reps <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_diagnostico_reps_gbs.csv")
    )
    
    archivo_diag_intentos <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_diagnostico_intentos_gbs.csv")
    )
    
    if (file.exists(archivo_mejor_rds) && !sobrescribir_hp_gbs) {
      
      mejor_grid_existente <- readRDS(archivo_mejor_rds)
      
      registro_hp_gbs_super$estado[i] <- "ya_existia"
      registro_hp_gbs_super$rep_piloto_solicitada[i] <- if ("rep_piloto_solicitada" %in% names(mejor_grid_existente)) mejor_grid_existente$rep_piloto_solicitada[1] else NA_integer_
      registro_hp_gbs_super$rep_piloto_usada[i] <- mejor_grid_existente$rep_piloto_usada[1]
      registro_hp_gbs_super$k_solicitado[i] <- if ("k_solicitado" %in% names(mejor_grid_existente)) mejor_grid_existente$k_solicitado[1] else NA_integer_
      registro_hp_gbs_super$k_cv_usado[i] <- mejor_grid_existente$k_cv_usado[1]
      registro_hp_gbs_super$eventos_piloto[i] <- mejor_grid_existente$eventos_piloto[1]
      registro_hp_gbs_super$censurados_piloto[i] <- mejor_grid_existente$censurados_piloto[1]
      
      write.csv(
        registro_hp_gbs_super,
        file.path(
          ruta_super,
          "hiperparametros_gbs",
          "registro_hp_gbs_super.csv"
        ),
        row.names = FALSE
      )
      
      cat("Hiperparámetros GBS ya existentes. Se omite:", nombre_esc, "\n")
      
      next
    }
    
    if (is.null(lista_esc)) {
      
      registro_hp_gbs_super$estado[i] <- "fallo"
      registro_hp_gbs_super$mensaje_error[i] <- "No existe el escenario en escenarios_24_listas."
      
      write.csv(
        registro_hp_gbs_super,
        file.path(
          ruta_super,
          "hiperparametros_gbs",
          "registro_hp_gbs_super.csv"
        ),
        row.names = FALSE
      )
      
      next
    }
    
    seleccion_hp <- tryCatch(
      seleccionar_hiperparametros_gbs_super(
        lista_escenario = lista_esc,
        grid_gbs = grid_gbs,
        rep_piloto = 1,
        k = 5,
        seed = 123,
        reps_candidatas = NULL,
        k_candidatos = c(5, 3, 2),
        max_reps_candidatas = 20
      ),
      error = function(e) {
        list(
          error = TRUE,
          mensaje = conditionMessage(e)
        )
      }
    )
    
    if (!is.null(seleccion_hp$error) && seleccion_hp$error) {
      
      registro_hp_gbs_super$estado[i] <- "fallo"
      registro_hp_gbs_super$mensaje_error[i] <- seleccion_hp$mensaje
      
      write.csv(
        registro_hp_gbs_super,
        file.path(
          ruta_super,
          "hiperparametros_gbs",
          "registro_hp_gbs_super.csv"
        ),
        row.names = FALSE
      )
      
      cat(
        "Falló optimización GBS en",
        nombre_esc,
        ":",
        seleccion_hp$mensaje,
        "\n"
      )
      
      next
    }
    
    mejor_grid <- seleccion_hp$mejor_grid
    grid_cv <- seleccion_hp$grid_cv
    
    saveRDS(
      mejor_grid,
      archivo_mejor_rds
    )
    
    write.csv(
      mejor_grid,
      archivo_mejor_csv,
      row.names = FALSE
    )
    
    write.csv(
      grid_cv,
      archivo_grid_csv,
      row.names = FALSE
    )
    
    if (!is.null(seleccion_hp$diagnostico_repeticiones)) {
      write.csv(
        seleccion_hp$diagnostico_repeticiones,
        archivo_diag_reps,
        row.names = FALSE
      )
    }
    
    if (!is.null(seleccion_hp$diagnostico_intentos)) {
      write.csv(
        seleccion_hp$diagnostico_intentos,
        archivo_diag_intentos,
        row.names = FALSE
      )
    }
    
    registro_hp_gbs_super$estado[i] <- "ok"
    registro_hp_gbs_super$rep_piloto_solicitada[i] <- mejor_grid$rep_piloto_solicitada[1]
    registro_hp_gbs_super$rep_piloto_usada[i] <- mejor_grid$rep_piloto_usada[1]
    registro_hp_gbs_super$k_solicitado[i] <- mejor_grid$k_solicitado[1]
    registro_hp_gbs_super$k_cv_usado[i] <- mejor_grid$k_cv_usado[1]
    registro_hp_gbs_super$eventos_piloto[i] <- mejor_grid$eventos_piloto[1]
    registro_hp_gbs_super$censurados_piloto[i] <- mejor_grid$censurados_piloto[1]
    
    write.csv(
      registro_hp_gbs_super,
      file.path(
        ruta_super,
        "hiperparametros_gbs",
        "registro_hp_gbs_super.csv"
      ),
      row.names = FALSE
    )
    
    cat("Hiperparámetros GBS guardados para:", nombre_esc, "\n")
    
    rm(mejor_grid, grid_cv, seleccion_hp)
    invisible(gc())
  }
})

tiempo_hp_gbs_super

registro_hp_gbs_super %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

registro_hp_gbs_super %>%
  dplyr::count(estado) %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

registro_hp_gbs_super %>%
  dplyr::filter(estado == "fallo") %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)



### 8.5.9 Verificación final hiperparámetros y frecuencia de hiperparámetros en cada escenario

#### 8.5.9.1 Verificar hiperparámetros seleccionados

leer_tabla_hiperparametros_modelo <- function(ruta_hp, sufijo_archivo, modelo_nombre) {
  
  archivos <- list.files(
    path = ruta_hp,
    pattern = paste0(sufijo_archivo, "\\.rds$"),
    full.names = TRUE
  )
  
  purrr::map_dfr(
    archivos,
    function(archivo) {
      
      nombre_esc <- basename(archivo) %>%
        stringr::str_remove(paste0(sufijo_archivo, "\\.rds$"))
      
      readRDS(archivo) %>%
        dplyr::mutate(
          escenario = nombre_esc,
          modelo = modelo_nombre,
          .before = 1
        )
    }
  )
}

hiperparametros_srf_super <- leer_tabla_hiperparametros_modelo(
  ruta_hp = file.path(ruta_super, "hiperparametros_srf"),
  sufijo_archivo = "_mejor_grid_srf",
  modelo_nombre = "SRF"
)

hiperparametros_gbs_super <- leer_tabla_hiperparametros_modelo(
  ruta_hp = file.path(ruta_super, "hiperparametros_gbs"),
  sufijo_archivo = "_mejor_grid_gbs",
  modelo_nombre = "GBS"
)

hiperparametros_srf_super
hiperparametros_gbs_super

dim(hiperparametros_srf_super)
dim(hiperparametros_gbs_super)

escenarios_hp_esperados <- tabla_24_escenarios$escenario_nombre

escenarios_hp_srf <- unique(hiperparametros_srf_super$escenario)
escenarios_hp_gbs <- unique(hiperparametros_gbs_super$escenario)

escenarios_srf_faltantes <- setdiff(
  escenarios_hp_esperados,
  escenarios_hp_srf
)

escenarios_gbs_faltantes <- setdiff(
  escenarios_hp_esperados,
  escenarios_hp_gbs
)

cat("\nEscenarios SRF con hiperparámetros:", length(escenarios_hp_srf), "\n")
cat("Escenarios GBS con hiperparámetros:", length(escenarios_hp_gbs), "\n")

if (length(escenarios_srf_faltantes) > 0) {
  cat("\nEscenarios SRF faltantes:\n")
  print(escenarios_srf_faltantes)
}

if (length(escenarios_gbs_faltantes) > 0) {
  cat("\nEscenarios GBS faltantes:\n")
  print(escenarios_gbs_faltantes)
}

stopifnot(length(escenarios_hp_srf) == 24)
stopifnot(length(escenarios_hp_gbs) == 24)
stopifnot(length(escenarios_srf_faltantes) == 0)
stopifnot(length(escenarios_gbs_faltantes) == 0)


#### 8.5.9.2 Frecuencia de hiperparámetros seleccionados en SRF

tabla_hp_srf_super <- hiperparametros_srf_super %>%
  dplyr::count(
    mtry,
    nodesize,
    ntree,
    splitrule,
    sort = TRUE
  )

tabla_hp_srf_super


#### 8.5.9.3 Frecuencia de hiperparámetros seleccionados en GBS

tabla_hp_gbs_super <- hiperparametros_gbs_super %>%
  dplyr::count(
    eta,
    max_depth,
    nrounds,
    best_nrounds,
    subsample,
    colsample_bytree,
    sort = TRUE
  )

tabla_hp_gbs_super


#### 8.5.9.4 Guardar tablas de hiperparámetros seleccionados

write.csv(
  hiperparametros_srf_super,
  file.path(ruta_super, "hiperparametros_srf", "hiperparametros_srf_super_todos.csv"),
  row.names = FALSE
)

write.csv(
  hiperparametros_gbs_super,
  file.path(ruta_super, "hiperparametros_gbs", "hiperparametros_gbs_super_todos.csv"),
  row.names = FALSE
)

write.csv(
  tabla_hp_srf_super,
  file.path(ruta_super, "hiperparametros_srf", "frecuencia_hiperparametros_srf_super.csv"),
  row.names = FALSE
)

write.csv(
  tabla_hp_gbs_super,
  file.path(ruta_super, "hiperparametros_gbs", "frecuencia_hiperparametros_gbs_super.csv"),
  row.names = FALSE
)




### 8.5.10 Funciones para validación externa en línea

crear_time_grid_super <- function(super_data,
                                  n_puntos = 40,
                                  t_min = 5,
                                  admin_time = 120) {
  
  # Grilla global fija para garantizar comparabilidad del IBS
  # entre escenarios con diferentes proporciones de censura.
  
  
  if (!is.finite(t_min) || !is.finite(admin_time) || admin_time <= t_min) {
    stop("La grilla global no es válida: revise t_min y admin_time.")
  }
  
  seq(
    from = t_min,
    to = admin_time,
    length.out = n_puntos
  )
}


predecir_supervivencia_cox_eficiente <- function(modelo_cox,
                                                 newdata,
                                                 times) {
  
  lp <- predict(
    modelo_cox,
    newdata = newdata,
    type = "lp",
    reference = "zero"
  )
  
  bh <- survival::basehaz(
    modelo_cox,
    centered = FALSE
  )
  
  H0_t <- sapply(times, function(tt) {
    
    idx <- which(bh$time <= tt)
    
    if (length(idx) == 0) {
      return(0)
    } else {
      return(bh$hazard[max(idx)])
    }
  })
  
  surv_mat <- outer(
    exp(lp),
    H0_t,
    function(riesgo, h0) exp(-riesgo * h0)
  )
  
  surv_mat <- pmin(pmax(surv_mat, 0), 1)
  
  list(
    surv = surv_mat,
    risk = lp
  )
}






### 8.5.11 Funciones para calcular IBS con pec() en la validación externa


limpiar_mensaje_error_validacion <- function(e) {
  
  x <- conditionMessage(e)
  x <- gsub("\033\\[[0-9;]*[A-Za-z]", "", x)
  x <- gsub("\\\\u001b\\[[0-9;]*[A-Za-z]", "", x)
  
  x
}


normalizar_evento_01_validacion <- function(x) {
  
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  if (is.logical(x)) {
    return(as.integer(x))
  }
  
  if (is.character(x)) {
    
    x <- dplyr::case_when(
      x %in% c("1", "evento", "Evento", "muerto", "Muerto", "dead", "Dead", "TRUE", "true") ~ 1L,
      x %in% c("0", "censura", "Censura", "censurado", "Censurado", "vivo", "Vivo", "alive", "Alive", "FALSE", "false") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(x))
    )
    
    return(as.integer(x))
  }
  
  x <- as.numeric(x)
  
  if (all(na.omit(unique(x)) %in% c(0, 1))) {
    return(as.integer(x))
  }
  
  if (all(na.omit(unique(x)) %in% c(1, 2))) {
    x <- ifelse(x == 2, 1, 0)
    return(as.integer(x))
  }
  
  as.integer(x)
}


estado_metricas_validacion <- function(ibs, c_index) {
  
  if (is.na(ibs) || !is.finite(ibs)) {
    return("fallo_ibs")
  }
  
  if (is.na(c_index) || !is.finite(c_index)) {
    return("fallo_cindex")
  }
  
  return("ok")
}



#### 8.5.11.1 Objeto auxiliar para usar pec() con matrices de supervivencia precalculadas

predictSurvProb.survprob_precalculada_pec <- function(object,
                                                      newdata,
                                                      times,
                                                      ...) {
  
  surv_ref <- as.matrix(object$surv)
  times_ref <- as.numeric(object$times)
  
  if (nrow(surv_ref) != nrow(newdata)) {
    stop("El número de filas de la matriz de supervivencia no coincide con newdata.")
  }
  
  if (ncol(surv_ref) != length(times_ref)) {
    stop("El número de columnas de la matriz de supervivencia no coincide con times_ref.")
  }
  
  idx <- findInterval(
    x = as.numeric(times),
    vec = times_ref
  )
  
  idx[idx < 1] <- 1
  idx[idx > length(times_ref)] <- length(times_ref)
  
  surv_out <- surv_ref[, idx, drop = FALSE]
  surv_out <- pmin(pmax(surv_out, 0), 1)
  
  surv_out
}


extraer_ibs_crps_validacion <- function(crps_obj,
                                        nombre_modelo = "Modelo") {
  
  val <- NA_real_
  
  if (is.matrix(crps_obj) || is.data.frame(crps_obj)) {
    
    if (nombre_modelo %in% colnames(crps_obj)) {
      val <- as.numeric(crps_obj[nrow(crps_obj), nombre_modelo])
    } else if (nombre_modelo %in% rownames(crps_obj)) {
      val <- as.numeric(crps_obj[nombre_modelo, ncol(crps_obj)])
    } else {
      val <- suppressWarnings(as.numeric(crps_obj[length(crps_obj)]))
    }
    
  } else {
    
    if (nombre_modelo %in% names(crps_obj)) {
      val <- suppressWarnings(as.numeric(crps_obj[[nombre_modelo]]))
    } else {
      val <- suppressWarnings(as.numeric(crps_obj[length(crps_obj)]))
    }
  }
  
  val
}


calcular_ibs_pec_validacion <- function(data_eval,
                                        surv_prob_mat,
                                        times_eval) {
  
  if (!requireNamespace("pec", quietly = TRUE)) {
    return(
      list(
        ibs = NA_real_,
        estado_ibs = "pec_no_disponible",
        mensaje_ibs = "El paquete pec no está instalado o no carga."
      )
    )
  }
  
  if (!requireNamespace("prodlim", quietly = TRUE)) {
    return(
      list(
        ibs = NA_real_,
        estado_ibs = "prodlim_no_disponible",
        mensaje_ibs = "El paquete prodlim no está instalado o no carga."
      )
    )
  }
  
  data_eval <- as.data.frame(data_eval)
  data_eval$event <- normalizar_evento_01_validacion(data_eval$event)
  data_eval$time <- as.numeric(data_eval$time)
  
  idx_ok <- which(
    !is.na(data_eval$time) &
      !is.na(data_eval$event) &
      data_eval$time > 0 &
      data_eval$event %in% c(0, 1)
  )
  
  data_eval <- data_eval[idx_ok, , drop = FALSE]
  surv_prob_mat <- as.matrix(surv_prob_mat)
  
  if (nrow(surv_prob_mat) >= max(idx_ok)) {
    surv_prob_mat <- surv_prob_mat[idx_ok, , drop = FALSE]
  }
  
  if (nrow(surv_prob_mat) != nrow(data_eval) &&
      ncol(surv_prob_mat) == nrow(data_eval)) {
    surv_prob_mat <- t(surv_prob_mat)
  }
  
  if (nrow(surv_prob_mat) != nrow(data_eval)) {
    return(
      list(
        ibs = NA_real_,
        estado_ibs = "fallo_dimension_surv",
        mensaje_ibs = "La matriz de supervivencia no coincide con el número de sujetos."
      )
    )
  }
  
  if (ncol(surv_prob_mat) != length(times_eval)) {
    return(
      list(
        ibs = NA_real_,
        estado_ibs = "fallo_dimension_times",
        mensaje_ibs = "La matriz de supervivencia no coincide con la grilla de tiempos."
      )
    )
  }
  
  objeto_pec <- structure(
    list(
      surv = surv_prob_mat,
      times = as.numeric(times_eval)
    ),
    class = "survprob_precalculada_pec"
  )
  
  avisos_pec <- character()
  
  out <- tryCatch(
    {
      pec_obj <- withCallingHandlers(
        pec::pec(
          object = list(Modelo = objeto_pec),
          formula = Hist(time, event) ~ 1,
          data = data_eval,
          times = as.numeric(times_eval),
          cens.model = "marginal",
          splitMethod = "none",
          exact = FALSE,
          verbose = FALSE
        ),
        warning = function(w) {
          avisos_pec <<- c(avisos_pec, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      
      crps_obj <- withCallingHandlers(
        pec::crps(
          pec_obj,
          times = max(times_eval, na.rm = TRUE)
        ),
        warning = function(w) {
          avisos_pec <<- c(avisos_pec, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      
      ibs_val <- extraer_ibs_crps_validacion(
        crps_obj = crps_obj,
        nombre_modelo = "Modelo"
      )
      
      mensaje_avisos <- ifelse(
        length(avisos_pec) == 0,
        NA_character_,
        paste(unique(avisos_pec), collapse = " | ")
      )
      
      list(
        ibs = ibs_val,
        estado_ibs = ifelse(is.na(ibs_val), "fallo_ibs_pec", "ok"),
        mensaje_ibs = mensaje_avisos
      )
    },
    error = function(e) {
      list(
        ibs = NA_real_,
        estado_ibs = "fallo_ibs_pec",
        mensaje_ibs = limpiar_mensaje_error_validacion(e)
      )
    }
  )
  
  out
}


#### 8.5.11.2 Validación externa CoxPH

validar_cox_superpob <- function(obj_modelo,
                                 train_data,
                                 super_data) {
  
  if (is.null(obj_modelo$modelo) || obj_modelo$estado_modelo != "ok") {
    return(tibble::tibble(
      modelo = "CoxPH",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = obj_modelo$estado_modelo,
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  datos_ok <- alinear_factores_modelo(
    train_data = train_data,
    test_data = super_data
  )
  
  test_ok <- datos_ok$test
  test_ok$event <- normalizar_evento_01_validacion(test_ok$event)
  
  times_eval <- crear_time_grid_super(
    super_data = test_ok,
    n_puntos = 40
  )
  
  pred_test <- tryCatch(
    predecir_supervivencia_cox_eficiente(
      modelo_cox = obj_modelo$modelo,
      newdata = test_ok,
      times = times_eval
    ),
    error = function(e) NULL
  )
  
  if (is.null(pred_test)) {
    return(tibble::tibble(
      modelo = "CoxPH",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = "fallo_prediccion",
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  ibs_out <- calcular_ibs_pec_validacion(
    data_eval = test_ok,
    surv_prob_mat = pred_test$surv,
    times_eval = times_eval
  )
  
  c_index <- tryCatch(
    survival::concordance(
      survival::Surv(test_ok$time, test_ok$event) ~ pred_test$risk,
      reverse = TRUE
    )$concordance,
    error = function(e) NA_real_
  )
  
  tibble::tibble(
    modelo = "CoxPH",
    ibs = ibs_out$ibs,
    c_index = c_index,
    estado_ajuste = estado_metricas_validacion(
      ibs = ibs_out$ibs,
      c_index = c_index
    ),
    metodo_ibs = "pec_crps",
    mensaje_ibs = ibs_out$mensaje_ibs
  )
}


#### 8.5.11.3 Validación externa SRF

validar_srf_superpob <- function(obj_modelo,
                                 train_data,
                                 super_data) {
  
  if (is.null(obj_modelo$modelo) || obj_modelo$estado_modelo != "ok") {
    return(tibble::tibble(
      modelo = "SRF",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = obj_modelo$estado_modelo,
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  datos_ok <- alinear_factores_modelo(
    train_data = train_data,
    test_data = super_data
  )
  
  test_ok <- datos_ok$test
  test_ok$event <- normalizar_evento_01_validacion(test_ok$event)
  
  times_eval <- crear_time_grid_super(
    super_data = test_ok,
    n_puntos = 40
  )
  
  pred_test <- tryCatch(
    predecir_supervivencia_srf(
      modelo_srf = obj_modelo$modelo,
      newdata = test_ok,
      times = times_eval
    ),
    error = function(e) NULL
  )
  
  if (is.null(pred_test)) {
    return(tibble::tibble(
      modelo = "SRF",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = "fallo_prediccion",
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  ibs_out <- calcular_ibs_pec_validacion(
    data_eval = test_ok,
    surv_prob_mat = pred_test$surv,
    times_eval = times_eval
  )
  
  c_index <- tryCatch(
    survival::concordance(
      survival::Surv(test_ok$time, test_ok$event) ~ pred_test$risk,
      reverse = TRUE
    )$concordance,
    error = function(e) NA_real_
  )
  
  tibble::tibble(
    modelo = "SRF",
    ibs = ibs_out$ibs,
    c_index = c_index,
    estado_ajuste = estado_metricas_validacion(
      ibs = ibs_out$ibs,
      c_index = c_index
    ),
    metodo_ibs = "pec_crps",
    mensaje_ibs = ibs_out$mensaje_ibs
  )
}



#### 8.5.11.4 Validación externa GBS

validar_gbs_superpob <- function(obj_modelo,
                                 train_data,
                                 super_data) {
  
  if (is.null(obj_modelo$modelo) || obj_modelo$estado_modelo != "ok") {
    return(tibble::tibble(
      modelo = "GBS",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = obj_modelo$estado_modelo,
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  datos_ok <- alinear_factores_modelo(
    train_data = train_data,
    test_data = super_data
  )
  
  test_ok <- datos_ok$test
  test_ok$event <- normalizar_evento_01_validacion(test_ok$event)
  
  xy_test <- crear_xy_xgb(test_ok)
  
  x_test <- alinear_columnas_xgb_por_nombres(
    x_new = xy_test$x,
    columnas_ref = obj_modelo$columnas_x
  )
  
  dtest <- xgboost::xgb.DMatrix(
    data = x_test,
    label = xy_test$y
  )
  
  lp_test <- tryCatch(
    predict(obj_modelo$modelo, dtest),
    error = function(e) NULL
  )
  
  if (is.null(lp_test)) {
    return(tibble::tibble(
      modelo = "GBS",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = "fallo_prediccion",
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  times_eval <- crear_time_grid_super(
    super_data = test_ok,
    n_puntos = 40
  )
  
  surv_test <- tryCatch(
    predecir_supervivencia_desde_lp(
      lp = lp_test,
      bh = obj_modelo$bh,
      times = times_eval
    ),
    error = function(e) NULL
  )
  
  if (is.null(surv_test)) {
    return(tibble::tibble(
      modelo = "GBS",
      ibs = NA_real_,
      c_index = NA_real_,
      estado_ajuste = "fallo_supervivencia",
      metodo_ibs = "pec_crps",
      mensaje_ibs = NA_character_
    ))
  }
  
  ibs_out <- calcular_ibs_pec_validacion(
    data_eval = test_ok,
    surv_prob_mat = surv_test,
    times_eval = times_eval
  )
  
  c_index <- tryCatch(
    survival::concordance(
      survival::Surv(test_ok$time, test_ok$event) ~ lp_test,
      reverse = TRUE
    )$concordance,
    error = function(e) NA_real_
  )
  
  tibble::tibble(
    modelo = "GBS",
    ibs = ibs_out$ibs,
    c_index = c_index,
    estado_ajuste = estado_metricas_validacion(
      ibs = ibs_out$ibs,
      c_index = c_index
    ),
    metodo_ibs = "pec_crps",
    mensaje_ibs = ibs_out$mensaje_ibs
  )
}



