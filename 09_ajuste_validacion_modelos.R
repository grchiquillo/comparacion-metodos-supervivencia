################################################################################

# 9. Ajuste y validación de modelos

#Se ajustan los modelos entrenando con el 100% de cada repetición simulada y 
#validando con la población externa de n=10000 que se generó para cada escenario.


carpeta_validacion_en_linea <- file.path(
  ruta_super,
  "validacion_en_linea_24_ibs_pec"
)

carpeta_resultados_consolidados <- file.path(
  ruta_super,
  "resultados_consolidados_24_ibs_pec"
)

dir.create(carpeta_validacion_en_linea, showWarnings = FALSE, recursive = TRUE)
dir.create(carpeta_resultados_consolidados, showWarnings = FALSE, recursive = TRUE)


valor_1 <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA)
  }
  x[1]
}


extraer_hiperparametros_usados <- function(nombre_modelo,
                                           mejor_grid = NULL,
                                           ajuste = NULL) {
  
  if (nombre_modelo == "cox") {
    return(tibble::tibble(
      mtry = NA_real_,
      nodesize = NA_real_,
      ntree = NA_real_,
      splitrule = NA_character_,
      eta = NA_real_,
      max_depth = NA_real_,
      nrounds = NA_real_,
      subsample = NA_real_,
      colsample_bytree = NA_real_,
      best_nrounds = NA_real_,
      nrounds_usado = NA_real_
    ))
  }
  
  if (nombre_modelo == "srf") {
    return(tibble::tibble(
      mtry = valor_1(mejor_grid$mtry),
      nodesize = valor_1(mejor_grid$nodesize),
      ntree = valor_1(mejor_grid$ntree),
      splitrule = as.character(valor_1(mejor_grid$splitrule)),
      eta = NA_real_,
      max_depth = NA_real_,
      nrounds = NA_real_,
      subsample = NA_real_,
      colsample_bytree = NA_real_,
      best_nrounds = NA_real_,
      nrounds_usado = NA_real_
    ))
  }
  
  if (nombre_modelo == "gbs") {
    return(tibble::tibble(
      mtry = NA_real_,
      nodesize = NA_real_,
      ntree = NA_real_,
      splitrule = NA_character_,
      eta = valor_1(mejor_grid$eta),
      max_depth = valor_1(mejor_grid$max_depth),
      nrounds = valor_1(mejor_grid$nrounds),
      subsample = valor_1(mejor_grid$subsample),
      colsample_bytree = valor_1(mejor_grid$colsample_bytree),
      best_nrounds = if ("best_nrounds" %in% names(mejor_grid)) {
        valor_1(mejor_grid$best_nrounds)
      } else {
        NA_real_
      },
      nrounds_usado = if (!is.null(ajuste$nrounds_usado)) {
        ajuste$nrounds_usado
      } else {
        NA_real_
      }
    ))
  }
  
  stop("Modelo no reconocido.")
}


leer_hiperparametros_modelo <- function(nombre_modelo,
                                        nombre_esc) {
  
  if (nombre_modelo == "cox") {
    return(NULL)
  }
  
  if (nombre_modelo == "srf") {
    archivo_hp <- file.path(
      ruta_super,
      "hiperparametros_srf",
      paste0(nombre_esc, "_mejor_grid_srf.rds")
    )
  }
  
  if (nombre_modelo == "gbs") {
    archivo_hp <- file.path(
      ruta_super,
      "hiperparametros_gbs",
      paste0(nombre_esc, "_mejor_grid_gbs.rds")
    )
  }
  
  if (!file.exists(archivo_hp)) {
    stop("No existe el archivo de hiperparámetros: ", archivo_hp)
  }
  
  readRDS(archivo_hp)
}


ajustar_y_validar_repeticion <- function(nombre_modelo,
                                         train_data,
                                         super_data,
                                         mejor_grid,
                                         rep_id,
                                         escenario_id) {
  
  etiqueta_modelo <- dplyr::case_when(
    nombre_modelo == "cox" ~ "CoxPH",
    nombre_modelo == "srf" ~ "SRF",
    nombre_modelo == "gbs" ~ "GBS",
    TRUE ~ nombre_modelo
  )
  
  set.seed(900000 + escenario_id * 1000 + rep_id)
  
  out <- tryCatch(
    {
      if (nombre_modelo == "cox") {
        
        ajuste <- ajustar_cox_modelo(train_data)
        
        obj_modelo <- list(
          rep = rep_id,
          modelo = ajuste$modelo,
          estado_modelo = ajuste$estado_modelo
        )
        
        metricas <- validar_cox_superpob(
          obj_modelo = obj_modelo,
          train_data = train_data,
          super_data = super_data
        )
      }
      
      if (nombre_modelo == "srf") {
        
        ajuste <- ajustar_srf_modelo_fijo(
          train_data = train_data,
          mejor_grid = mejor_grid
        )
        
        obj_modelo <- list(
          rep = rep_id,
          modelo = ajuste$modelo,
          estado_modelo = ajuste$estado_modelo,
          mejor_grid = mejor_grid
        )
        
        metricas <- validar_srf_superpob(
          obj_modelo = obj_modelo,
          train_data = train_data,
          super_data = super_data
        )
      }
      
      if (nombre_modelo == "gbs") {
        
        ajuste <- ajustar_gbs_modelo_fijo(
          train_data = train_data,
          mejor_grid = mejor_grid
        )
        
        obj_modelo <- list(
          rep = rep_id,
          modelo = ajuste$modelo,
          estado_modelo = ajuste$estado_modelo,
          mejor_grid = mejor_grid,
          bh = ajuste$bh,
          columnas_x = ajuste$columnas_x,
          nrounds_usado = ajuste$nrounds_usado
        )
        
        metricas <- validar_gbs_superpob(
          obj_modelo = obj_modelo,
          train_data = train_data,
          super_data = super_data
        )
      }
      
      hiper <- extraer_hiperparametros_usados(
        nombre_modelo = nombre_modelo,
        mejor_grid = mejor_grid,
        ajuste = ajuste
      )
      
      dplyr::bind_cols(metricas, hiper)
    },
    error = function(e) {
      
      hiper <- extraer_hiperparametros_usados(
        nombre_modelo = nombre_modelo,
        mejor_grid = mejor_grid,
        ajuste = NULL
      )
      
      dplyr::bind_cols(
        tibble::tibble(
          modelo = etiqueta_modelo,
          ibs = NA_real_,
          c_index = NA_real_,
          estado_ajuste = "fallo_repeticion",
          metodo_ibs = "pec_crps",
          mensaje_ibs = NA_character_,
          mensaje_error = conditionMessage(e)
        ),
        hiper
      )
    }
  )
  
  if (exists("ajuste")) rm(ajuste)
  if (exists("obj_modelo")) rm(obj_modelo)
  invisible(gc())
  
  out
}


ajustar_y_validar_escenario_modelo <- function(nombre_modelo,
                                               nombre_esc,
                                               sobrescribir = FALSE,
                                               guardar_cada = 10) {
  
  cat("\n====================================\n")
  cat("Modelo:", toupper(nombre_modelo), "\n")
  cat("Escenario:", nombre_esc, "\n")
  cat("====================================\n")
  
  carpeta_val <- file.path(
    carpeta_validacion_en_linea,
    nombre_modelo
  )
  
  dir.create(carpeta_val, showWarnings = FALSE, recursive = TRUE)
  
  archivo_resultado_rds <- file.path(
    carpeta_val,
    paste0(nombre_esc, "_", nombre_modelo, "_validacion_en_linea.rds")
  )
  
  archivo_resultado_csv <- file.path(
    carpeta_val,
    paste0(nombre_esc, "_", nombre_modelo, "_validacion_en_linea.csv")
  )
  
  resultados_previos <- tibble::tibble()
  reps_ya <- integer(0)
  
  if (file.exists(archivo_resultado_rds) && !sobrescribir) {
    
    resultados_previos <- readRDS(archivo_resultado_rds)
    reps_ya <- unique(resultados_previos$rep)
    
    if (length(reps_ya) >= 200) {
      cat("El escenario ya tiene 200 repeticiones. Se omite.\n")
      return(resultados_previos)
    }
  }
  
  archivo_escenario <- file.path(
    "escenarios_agrupados",
    paste0(nombre_esc, ".rds")
  )
  
  archivo_super <- file.path(
    ruta_super,
    "superpoblaciones",
    paste0(nombre_esc, "_superpob_n10000.rds")
  )
  
  if (!file.exists(archivo_escenario)) {
    stop("No existe el archivo del escenario: ", archivo_escenario)
  }
  
  if (!file.exists(archivo_super)) {
    stop("No existe la superpoblación: ", archivo_super)
  }
  
  lista_esc <- readRDS(archivo_escenario)
  super_data <- readRDS(archivo_super)
  
  mejor_grid <- leer_hiperparametros_modelo(
    nombre_modelo = nombre_modelo,
    nombre_esc = nombre_esc
  )
  
  esc_row <- tabla_24_escenarios[
    tabla_24_escenarios$escenario_nombre == nombre_esc,
  ]
  
  escenario_id <- esc_row$escenario_id[1]
  n_esc <- esc_row$n[1]
  censura_esc <- esc_row$censura_objetivo[1]
  tipo_esc <- esc_row$tipo_relacion[1]
  
  resultados_nuevos <- list()
  contador_nuevos <- 0
  
  for (j in seq_along(lista_esc)) {
    
    nombre_rep <- names(lista_esc)[j]
    
    rep_id <- suppressWarnings(
      as.integer(stringr::str_extract(nombre_rep, "\\d+"))
    )
    
    if (is.na(rep_id)) {
      rep_id <- j
    }
    
    if (rep_id %in% reps_ya && !sobrescribir) {
      cat("Saltando repetición", rep_id, "porque ya existe resultado.\n")
      next
    }
    
    cat(
      "Modelo:",
      toupper(nombre_modelo),
      "| escenario:",
      nombre_esc,
      "| repetición:",
      rep_id,
      "\n"
    )
    
    train_data <- lista_esc[[j]]
    
    tiempo_rep <- system.time({
      
      metricas_rep <- ajustar_y_validar_repeticion(
        nombre_modelo = nombre_modelo,
        train_data = train_data,
        super_data = super_data,
        mejor_grid = mejor_grid,
        rep_id = rep_id,
        escenario_id = escenario_id
      )
    })
    
    resultado_rep <- dplyr::bind_cols(
      tibble::tibble(
        escenario = nombre_esc,
        n = n_esc,
        censura_objetivo = censura_esc,
        tipo_relacion = tipo_esc,
        rep = rep_id,
        n_validacion = nrow(super_data),
        tiempo_ajuste_validacion_segundos = as.numeric(tiempo_rep["elapsed"])
      ),
      metricas_rep
    )
    
    contador_nuevos <- contador_nuevos + 1
    resultados_nuevos[[contador_nuevos]] <- resultado_rep
    
    rm(train_data, metricas_rep, resultado_rep)
    invisible(gc())
    
    if (contador_nuevos %% guardar_cada == 0) {
      
      resultados_parciales <- dplyr::bind_rows(
        resultados_previos,
        dplyr::bind_rows(resultados_nuevos)
      )
      
      saveRDS(resultados_parciales, archivo_resultado_rds)
      
      write.csv(
        resultados_parciales,
        archivo_resultado_csv,
        row.names = FALSE
      )
      
      cat("Guardado parcial:", contador_nuevos, "nuevas repeticiones.\n")
    }
  }
  
  resultados_esc <- dplyr::bind_rows(
    resultados_previos,
    dplyr::bind_rows(resultados_nuevos)
  )
  
  saveRDS(resultados_esc, archivo_resultado_rds)
  
  write.csv(
    resultados_esc,
    archivo_resultado_csv,
    row.names = FALSE
  )
  
  rm(lista_esc, super_data, mejor_grid, resultados_previos, resultados_nuevos)
  invisible(gc())
  
  return(resultados_esc)
}


ejecutar_ajuste_validacion_modelo <- function(nombre_modelo,
                                              sobrescribir = FALSE) {
  
  nombres_escenarios <- tabla_24_escenarios$escenario_nombre
  
  resultados_modelo <- vector("list", length(nombres_escenarios))
  names(resultados_modelo) <- nombres_escenarios
  
  for (i in seq_along(nombres_escenarios)) {
    
    nombre_esc <- nombres_escenarios[i]
    
    cat("\n####################################\n")
    cat("Modelo:", toupper(nombre_modelo), "\n")
    cat("Escenario", i, "de", length(nombres_escenarios), "\n")
    cat("Nombre:", nombre_esc, "\n")
    cat("####################################\n")
    
    resultados_modelo[[i]] <- tryCatch(
      ajustar_y_validar_escenario_modelo(
        nombre_modelo = nombre_modelo,
        nombre_esc = nombre_esc,
        sobrescribir = sobrescribir,
        guardar_cada = 10
      ),
      error = function(e) {
        
        message("Falló ", toupper(nombre_modelo), " en ", nombre_esc, ": ", e$message)
        
        esc_row <- tabla_24_escenarios[
          tabla_24_escenarios$escenario_nombre == nombre_esc,
        ]
        
        tibble::tibble(
          escenario = nombre_esc,
          n = esc_row$n[1],
          censura_objetivo = esc_row$censura_objetivo[1],
          tipo_relacion = esc_row$tipo_relacion[1],
          rep = NA_integer_,
          n_validacion = NA_integer_,
          tiempo_ajuste_validacion_segundos = NA_real_,
          modelo = toupper(nombre_modelo),
          ibs = NA_real_,
          c_index = NA_real_,
          estado_ajuste = "fallo_escenario",
          metodo_ibs = "pec_crps",
          mensaje_ibs = NA_character_,
          mensaje_error = e$message
        )
      }
    )
    
    invisible(gc())
  }
  
  resultados_modelo_df <- dplyr::bind_rows(resultados_modelo)
  
  saveRDS(
    resultados_modelo_df,
    file.path(
      carpeta_resultados_consolidados,
      paste0("resultados_", nombre_modelo, "_superpob_10000_en_linea.rds")
    )
  )
  
  write.csv(
    resultados_modelo_df,
    file.path(
      carpeta_resultados_consolidados,
      paste0("resultados_", nombre_modelo, "_superpob_10000_en_linea.csv")
    ),
    row.names = FALSE
  )
  
  return(resultados_modelo_df)
}



## 9.1 Ajuste y validación CoxPH



tiempo_cox_en_linea <- system.time({
  resultados_cox_super <- ejecutar_ajuste_validacion_modelo(
    nombre_modelo = "cox",
    sobrescribir = FALSE
  )
})

tiempo_cox_en_linea




## 9.2 Ajuste y validación SRF

tiempo_srf_en_linea <- system.time({
  resultados_srf_super <- ejecutar_ajuste_validacion_modelo(
    nombre_modelo = "srf",
    sobrescribir = FALSE
  )
})

tiempo_srf_en_linea



## 9.3 Ajuste y validación GBS

tiempo_gbs_en_linea <- system.time({
  resultados_gbs_super <- ejecutar_ajuste_validacion_modelo(
    nombre_modelo = "gbs",
    sobrescribir = FALSE
  )
})

tiempo_gbs_en_linea



### 9.4 Verificación del avance de ajuste + validación

verificar_avance_en_linea <- function(nombre_modelo) {
  
  carpeta_val <- file.path(
    carpeta_validacion_en_linea,
    nombre_modelo
  )
  
  purrr::map_dfr(
    tabla_24_escenarios$escenario_nombre,
    function(nombre_esc) {
      
      archivo_esc <- file.path(
        carpeta_val,
        paste0(nombre_esc, "_", nombre_modelo, "_validacion_en_linea.rds")
      )
      
      if (!file.exists(archivo_esc)) {
        return(tibble::tibble(
          modelo = nombre_modelo,
          escenario = nombre_esc,
          n_reps = 0,
          n_ok = 0,
          n_ibs = 0,
          n_c_index = 0,
          completo = FALSE
        ))
      }
      
      res <- readRDS(archivo_esc)
      
      tibble::tibble(
        modelo = nombre_modelo,
        escenario = nombre_esc,
        n_reps = dplyr::n_distinct(res$rep),
        n_ok = sum(res$estado_ajuste == "ok", na.rm = TRUE),
        n_ibs = sum(!is.na(res$ibs)),
        n_c_index = sum(!is.na(res$c_index)),
        completo = dplyr::n_distinct(res$rep) == 200
      )
    }
  )
}

avance_cox <- verificar_avance_en_linea("cox")
avance_srf <- verificar_avance_en_linea("srf")
avance_gbs <- verificar_avance_en_linea("gbs")

avance_general <- dplyr::bind_rows(
  avance_cox,
  avance_srf,
  avance_gbs
)

avance_general

avance_general %>%
  dplyr::filter(!completo)

write.csv(
  avance_general,
  file.path(
    carpeta_resultados_consolidados,
    "avance_ajuste_validacion_en_linea.csv"
  ),
  row.names = FALSE
)



## 9.5 Consolidación resultados

leer_resultados_en_linea <- function(nombre_modelo) {
  
  carpeta_val <- file.path(
    carpeta_validacion_en_linea,
    nombre_modelo
  )
  
  archivos <- list.files(
    carpeta_val,
    pattern = paste0("_", nombre_modelo, "_validacion_en_linea\\.rds$"),
    full.names = TRUE
  )
  
  if (length(archivos) == 0) {
    warning("No se encontraron archivos para el modelo: ", nombre_modelo)
    return(tibble::tibble())
  }
  
  purrr::map_dfr(archivos, readRDS) %>%
    dplyr::filter(
      escenario %in% tabla_24_escenarios$escenario_nombre,
      censura_objetivo %in% c(0.20, 0.50, 0.75)
    )
}


resultados_cox_super <- leer_resultados_en_linea("cox")
resultados_srf_super <- leer_resultados_en_linea("srf")
resultados_gbs_super <- leer_resultados_en_linea("gbs")

resultados_modelos_super <- dplyr::bind_rows(
  resultados_cox_super,
  resultados_srf_super,
  resultados_gbs_super
) %>%
  dplyr::mutate(
    modelo = factor(
      modelo,
      levels = c("CoxPH", "SRF", "GBS")
    )
  )

dim(resultados_modelos_super)

table(resultados_modelos_super$modelo, useNA = "ifany")
table(resultados_modelos_super$estado_ajuste, useNA = "ifany")
table(resultados_modelos_super$censura_objetivo, useNA = "ifany")
table(resultados_modelos_super$tipo_relacion, useNA = "ifany")

stopifnot(!("error_estimacion" %in% names(resultados_modelos_super)))
stopifnot(!any(resultados_modelos_super$censura_objetivo == 0.95, na.rm = TRUE))

saveRDS(
  resultados_modelos_super,
  file.path(
    carpeta_resultados_consolidados,
    "resultados_modelos_superpob_10000_ibs_pec.rds"
  )
)

write.csv(
  resultados_modelos_super,
  file.path(
    carpeta_resultados_consolidados,
    "resultados_modelos_superpob_10000_ibs_pec.csv"
  ),
  row.names = FALSE
)



## 9.6 Tabla resumen por escenario

media_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

sd_na <- function(x) {
  if (sum(!is.na(x)) <= 1) {
    return(NA_real_)
  }
  sd(x, na.rm = TRUE)
}

resultados_modelos_super_validos <- resultados_modelos_super %>%
  dplyr::filter(
    estado_ajuste == "ok",
    !is.na(ibs),
    !is.na(c_index),
    is.finite(ibs),
    is.finite(c_index)
  )

resultados_super_por_escenario <- resultados_modelos_super_validos %>%
  dplyr::group_by(
    escenario,
    n,
    censura_objetivo,
    tipo_relacion,
    n_validacion,
    modelo
  ) %>%
  dplyr::summarise(
    n_repeticiones_validas = dplyr::n(),
    
    media_ibs = media_na(ibs),
    sd_ibs = sd_na(ibs),
    mediana_ibs = median(ibs, na.rm = TRUE),
    p25_ibs = quantile(ibs, 0.25, na.rm = TRUE),
    p75_ibs = quantile(ibs, 0.75, na.rm = TRUE),
    
    media_c_index = media_na(c_index),
    sd_c_index = sd_na(c_index),
    mediana_c_index = median(c_index, na.rm = TRUE),
    p25_c_index = quantile(c_index, 0.25, na.rm = TRUE),
    p75_c_index = quantile(c_index, 0.75, na.rm = TRUE),
    
    metodo_ibs = {
      x <- unique(na.omit(as.character(metodo_ibs)))
      if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
    },
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    tipo_relacion,
    censura_objetivo,
    n,
    modelo
  )

resultados_super_por_escenario

write.csv(
  resultados_super_por_escenario,
  file.path(
    carpeta_resultados_consolidados,
    "resultados_superpob_10000_por_escenario_ibs_pec.csv"
  ),
  row.names = FALSE
)

saveRDS(
  resultados_super_por_escenario,
  file.path(
    carpeta_resultados_consolidados,
    "resultados_superpob_10000_por_escenario_ibs_pec.rds"
  )
)




## 9.7 Creación de Boxplots

preparar_factores_grafico <- function(data) {
  
  data %>%
    dplyr::mutate(
      modelo = factor(
        modelo,
        levels = c("CoxPH", "SRF", "GBS")
      ),
      tipo_relacion_f = factor(
        tipo_relacion,
        levels = c("lineal", "no_lineal"),
        labels = c("Lineal", "No lineal")
      ),
      n_f = factor(
        n,
        levels = c(100, 500, 1000, 5000),
        labels = c("n = 100", "n = 500", "n = 1000", "n = 5000")
      ),
      censura_f = factor(
        censura_objetivo,
        levels = c(0.20, 0.50, 0.75),
        labels = c(
          "Censura = 20%",
          "Censura = 50%",
          "Censura = 75%"
        )
      )
    )
}


# Base para graficar IBS.
# Aquí sí se exige IBS válido.

resultados_modelos_super_plot_ibs <- resultados_modelos_super %>%
  dplyr::filter(
    !is.na(ibs),
    is.finite(ibs)
  ) %>%
  preparar_factores_grafico()


# Base para graficar C-index.
# Aquí NO se exige IBS válido, porque el C-index puede existir aunque pec() falle.

resultados_modelos_super_plot_cindex <- resultados_modelos_super %>%
  dplyr::filter(
    !is.na(c_index),
    is.finite(c_index)
  ) %>%
  preparar_factores_grafico()


resultados_super_lineal_ibs <- resultados_modelos_super_plot_ibs %>%
  dplyr::filter(tipo_relacion == "lineal")

resultados_super_no_lineal_ibs <- resultados_modelos_super_plot_ibs %>%
  dplyr::filter(tipo_relacion == "no_lineal")


resultados_super_lineal_cindex <- resultados_modelos_super_plot_cindex %>%
  dplyr::filter(tipo_relacion == "lineal")

resultados_super_no_lineal_cindex <- resultados_modelos_super_plot_cindex %>%
  dplyr::filter(tipo_relacion == "no_lineal")


graficar_boxplot_super <- function(data,
                                   metrica,
                                   titulo_grafico,
                                   escalas = "fixed") {
  
  ggplot2::ggplot(
    data = data,
    ggplot2::aes(
      x = modelo,
      y = .data[[metrica]],
      fill = modelo
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.70,
      outlier.alpha = 0.35,
      outlier.size = 1.2
    ) +
    ggplot2::facet_grid(
      censura_f ~ n_f,
      scales = escalas,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = titulo_grafico,
      subtitle = "Validación externa en superpoblación independiente n = 10000",
      x = "Modelo",
      y = dplyr::case_when(
        metrica == "ibs" ~ "Integrated Brier Score (IBS)",
        metrica == "c_index" ~ "C-index",
        TRUE ~ metrica
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}

### 9.7.1 Boxplot IBS

plot_super_ibs_lineal <- graficar_boxplot_super(
  data = resultados_super_lineal_ibs,
  metrica = "ibs",
  titulo_grafico = "Comparación del IBS entre modelos en escenarios lineales"
) +
  ggplot2::coord_cartesian(ylim = c(0, 0.5))

plot_super_ibs_no_lineal <- graficar_boxplot_super(
  data = resultados_super_no_lineal_ibs,
  metrica = "ibs",
  titulo_grafico = "Comparación del IBS entre modelos en escenarios no lineales"
) +
  ggplot2::coord_cartesian(ylim = c(0, 0.5))

plot_super_ibs_lineal
plot_super_ibs_no_lineal


### 9.7.2 Boxplot C index

plot_super_cindex_lineal <- graficar_boxplot_super(
  data = resultados_super_lineal_cindex,
  metrica = "c_index",
  titulo_grafico = "Comparación del C-index entre modelos en escenarios lineales"
) +
  ggplot2::coord_cartesian(ylim = c(0.3, 1.0))

plot_super_cindex_no_lineal <- graficar_boxplot_super(
  data = resultados_super_no_lineal_cindex,
  metrica = "c_index",
  titulo_grafico = "Comparación del C-index entre modelos en escenarios no lineales"
) +
  ggplot2::coord_cartesian(ylim = c(0.3, 1.0))

plot_super_cindex_lineal
plot_super_cindex_no_lineal



### 9.7.3 Guardar figuras

carpeta_figuras_24 <- file.path(
  ruta_super,
  "figuras_24_ibs_pec"
)

dir.create(
  carpeta_figuras_24,
  showWarnings = FALSE,
  recursive = TRUE
)

ggplot2::ggsave(
  file.path(carpeta_figuras_24, "superpob_ibs_lineal.png"),
  plot_super_ibs_lineal,
  width = 13,
  height = 10,
  dpi = 300
)

ggplot2::ggsave(
  file.path(carpeta_figuras_24, "superpob_ibs_no_lineal.png"),
  plot_super_ibs_no_lineal,
  width = 13,
  height = 10,
  dpi = 300
)

ggplot2::ggsave(
  file.path(carpeta_figuras_24, "superpob_cindex_lineal.png"),
  plot_super_cindex_lineal,
  width = 13,
  height = 10,
  dpi = 300
)

ggplot2::ggsave(
  file.path(carpeta_figuras_24, "superpob_cindex_no_lineal.png"),
  plot_super_cindex_no_lineal,
  width = 13,
  height = 10,
  dpi = 300
)



## 9.8 Comparación estadistica entre modelos mediante Kruskal Wallis


### 9.8.1 Verificación mínima

stopifnot(exists("resultados_modelos_super"))


resultados_kw_base <- resultados_modelos_super %>%
  dplyr::filter(
    modelo %in% c("CoxPH", "SRF", "GBS"),
    censura_objetivo %in% c(0.20, 0.50, 0.75)
  ) %>%
  dplyr::mutate(
    modelo = factor(
      modelo,
      levels = c("CoxPH", "SRF", "GBS")
    ),
    tipo_relacion = factor(
      tipo_relacion,
      levels = c("lineal", "no_lineal")
    ),
    n = as.numeric(n),
    censura_objetivo = as.numeric(censura_objetivo)
  )

# Verificación de disponibilidad por escenario y modelo

resultados_kw_base %>%
  dplyr::group_by(
    tipo_relacion,
    n,
    censura_objetivo,
    modelo
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_ibs = sum(!is.na(ibs) & is.finite(ibs)),
    n_c_index = sum(!is.na(c_index) & is.finite(c_index)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    tipo_relacion,
    n,
    censura_objetivo,
    modelo
  ) %>%
  print(n = Inf, width = Inf)



### 9.8.2 Función general para Kruskal-Wallis


calcular_kw_por_escenario <- function(data,
                                      metrica) {
  
  data_metrica <- data %>%
    dplyr::filter(
      !is.na(.data[[metrica]]),
      is.finite(.data[[metrica]])
    )
  
  resultado_kw <- data_metrica %>%
    dplyr::group_by(
      tipo_relacion,
      n,
      censura_objetivo
    ) %>%
    dplyr::group_modify(
      function(.x, .y) {
        
        modelos_disponibles <- .x %>%
          dplyr::filter(
            !is.na(.data[[metrica]]),
            is.finite(.data[[metrica]])
          ) %>%
          dplyr::distinct(modelo) %>%
          dplyr::pull(modelo)
        
        n_modelos <- length(modelos_disponibles)
        
        if (n_modelos < 2) {
          return(
            tibble::tibble(
              metrica = metrica,
              estadistico_kw = NA_real_,
              gl = NA_real_,
              p_valor = NA_real_,
              n_total = nrow(.x),
              n_modelos = n_modelos,
              estado_kw = "menos_de_dos_modelos"
            )
          )
        }
        
        prueba <- tryCatch(
          stats::kruskal.test(
            formula = stats::as.formula(
              paste(metrica, "~ modelo")
            ),
            data = .x
          ),
          error = function(e) {
            NULL
          }
        )
        
        if (is.null(prueba)) {
          return(
            tibble::tibble(
              metrica = metrica,
              estadistico_kw = NA_real_,
              gl = NA_real_,
              p_valor = NA_real_,
              n_total = nrow(.x),
              n_modelos = n_modelos,
              estado_kw = "fallo_kw"
            )
          )
        }
        
        tibble::tibble(
          metrica = metrica,
          estadistico_kw = as.numeric(prueba$statistic),
          gl = as.numeric(prueba$parameter),
          p_valor = as.numeric(prueba$p.value),
          n_total = nrow(.x),
          n_modelos = n_modelos,
          estado_kw = "ok"
        )
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_valor_ajustado_bh = p.adjust(
        p_valor,
        method = "BH"
      ),
      significancia = dplyr::case_when(
        is.na(p_valor) ~ "No evaluable",
        p_valor < 0.001 ~ "***",
        p_valor < 0.01 ~ "**",
        p_valor < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      etiqueta_p = dplyr::case_when(
        is.na(p_valor) ~ "KW: no evaluable",
        p_valor < 0.001 ~ "KW: p < 0.001",
        TRUE ~ paste0("KW: p = ", format(round(p_valor, 3), nsmall = 3))
      )
    ) %>%
    dplyr::arrange(
      tipo_relacion,
      n,
      censura_objetivo
    )
  
  resultado_kw
}



### 9.8.3 Kruskal-Wallis para IBS y C-index


kw_ibs <- calcular_kw_por_escenario(
  data = resultados_kw_base,
  metrica = "ibs"
)

kw_c_index <- calcular_kw_por_escenario(
  data = resultados_kw_base,
  metrica = "c_index"
)

kw_ibs %>%
  print(n = Inf, width = Inf)

kw_c_index %>%
  print(n = Inf, width = Inf)

# Guardar resultados

carpeta_kw <- file.path(
  ruta_super,
  "comparacion_kruskal_wallis_24"
)

dir.create(
  carpeta_kw,
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  kw_ibs,
  file.path(
    carpeta_kw,
    "kruskal_wallis_ibs.csv"
  ),
  row.names = FALSE
)

write.csv(
  kw_c_index,
  file.path(
    carpeta_kw,
    "kruskal_wallis_c_index.csv"
  ),
  row.names = FALSE
)



### 9.8.4 Preparar datos para boxplots con p-valor

preparar_factores_grafico_kw <- function(data) {
  
  data %>%
    dplyr::mutate(
      modelo = factor(
        modelo,
        levels = c("CoxPH", "SRF", "GBS")
      ),
      tipo_relacion_f = factor(
        tipo_relacion,
        levels = c("lineal", "no_lineal"),
        labels = c("Lineal", "No lineal")
      ),
      n_f = factor(
        n,
        levels = c(100, 500, 1000, 5000),
        labels = c("n = 100", "n = 500", "n = 1000", "n = 5000")
      ),
      censura_f = factor(
        censura_objetivo,
        levels = c(0.20, 0.50, 0.75),
        labels = c(
          "Censura = 20%",
          "Censura = 50%",
          "Censura = 75%"
        )
      )
    )
}


preparar_p_kw_grafico <- function(tabla_kw,
                                  y_pos) {
  
  tabla_kw %>%
    dplyr::mutate(
      tipo_relacion_f = factor(
        tipo_relacion,
        levels = c("lineal", "no_lineal"),
        labels = c("Lineal", "No lineal")
      ),
      n_f = factor(
        n,
        levels = c(100, 500, 1000, 5000),
        labels = c("n = 100", "n = 500", "n = 1000", "n = 5000")
      ),
      censura_f = factor(
        censura_objetivo,
        levels = c(0.20, 0.50, 0.75),
        labels = c(
          "Censura = 20%",
          "Censura = 50%",
          "Censura = 75%"
        )
      ),
      x_pos = 2,
      y_pos = y_pos
    )
}


# Bases separadas para IBS y C-index

datos_plot_ibs_kw <- resultados_kw_base %>%
  dplyr::filter(
    !is.na(ibs),
    is.finite(ibs)
  ) %>%
  preparar_factores_grafico_kw()

datos_plot_cindex_kw <- resultados_kw_base %>%
  dplyr::filter(
    !is.na(c_index),
    is.finite(c_index)
  ) %>%
  preparar_factores_grafico_kw()


# Tablas de p-valor para anotar en los gráficos

p_ibs_graf <- preparar_p_kw_grafico(
  tabla_kw = kw_ibs,
  y_pos = 0.48
)

p_cindex_graf <- preparar_p_kw_grafico(
  tabla_kw = kw_c_index,
  y_pos = 0.97
)



### 9.8.5 Función para boxplot con p-valor

graficar_boxplot_kw <- function(data,
                                tabla_p,
                                metrica,
                                titulo_grafico,
                                y_label,
                                ylim_vals) {
  
  ggplot2::ggplot(
    data = data,
    ggplot2::aes(
      x = modelo,
      y = .data[[metrica]],
      fill = modelo
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.70,
      outlier.alpha = 0.35,
      outlier.size = 1.2
    ) +
    ggplot2::geom_text(
      data = tabla_p,
      ggplot2::aes(
        x = x_pos,
        y = y_pos,
        label = etiqueta_p
      ),
      inherit.aes = FALSE,
      size = 3.8
    ) +
    ggplot2::facet_grid(
      censura_f ~ n_f,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = ylim_vals
    ) +
    ggplot2::labs(
      title = titulo_grafico,
      subtitle = "Prueba global de Kruskal-Wallis entre CoxPH, SRF y GBS",
      x = "Modelo",
      y = y_label
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 13
      ),
      axis.title = ggplot2::element_text(
        size = 13
      ),
      axis.text = ggplot2::element_text(
        size = 11
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      legend.position = "none",
      panel.spacing = grid::unit(0.6, "lines")
    )
}


### 9.8.6 Boxplots IBS con p-valor de Kruskal-Wallis

plot_ibs_kw_lineal <- graficar_boxplot_kw(
  data = datos_plot_ibs_kw %>%
    dplyr::filter(tipo_relacion == "lineal"),
  tabla_p = p_ibs_graf %>%
    dplyr::filter(tipo_relacion == "lineal"),
  metrica = "ibs",
  titulo_grafico = "Comparación del IBS entre modelos en escenarios lineales",
  y_label = "Integrated Brier Score (IBS)",
  ylim_vals = c(0, 0.5)
)

plot_ibs_kw_no_lineal <- graficar_boxplot_kw(
  data = datos_plot_ibs_kw %>%
    dplyr::filter(tipo_relacion == "no_lineal"),
  tabla_p = p_ibs_graf %>%
    dplyr::filter(tipo_relacion == "no_lineal"),
  metrica = "ibs",
  titulo_grafico = "Comparación del IBS entre modelos en escenarios no lineales",
  y_label = "Integrated Brier Score (IBS)",
  ylim_vals = c(0, 0.5)
)

plot_ibs_kw_lineal
plot_ibs_kw_no_lineal




### 9.8.7 Boxplots C-index con p-valor de Kruskal-Wallis

plot_cindex_kw_lineal <- graficar_boxplot_kw(
  data = datos_plot_cindex_kw %>%
    dplyr::filter(tipo_relacion == "lineal"),
  tabla_p = p_cindex_graf %>%
    dplyr::filter(tipo_relacion == "lineal"),
  metrica = "c_index",
  titulo_grafico = "Comparación del C-index entre modelos en escenarios lineales",
  y_label = "C-index",
  ylim_vals = c(0.3, 1.0)
)

plot_cindex_kw_no_lineal <- graficar_boxplot_kw(
  data = datos_plot_cindex_kw %>%
    dplyr::filter(tipo_relacion == "no_lineal"),
  tabla_p = p_cindex_graf %>%
    dplyr::filter(tipo_relacion == "no_lineal"),
  metrica = "c_index",
  titulo_grafico = "Comparación del C-index entre modelos en escenarios no lineales",
  y_label = "C-index",
  ylim_vals = c(0.3, 1.0)
)

plot_cindex_kw_lineal
plot_cindex_kw_no_lineal


### 9.8.8 Guardar boxplots con p-valor

carpeta_figuras_kw <- file.path(
  ruta_super,
  "figuras_kruskal_wallis_24"
)

dir.create(
  carpeta_figuras_kw,
  showWarnings = FALSE,
  recursive = TRUE
)

ggplot2::ggsave(
  file.path(carpeta_figuras_kw, "boxplot_ibs_kw_lineal.png"),
  plot_ibs_kw_lineal,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(carpeta_figuras_kw, "boxplot_ibs_kw_no_lineal.png"),
  plot_ibs_kw_no_lineal,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(carpeta_figuras_kw, "boxplot_cindex_kw_lineal.png"),
  plot_cindex_kw_lineal,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(carpeta_figuras_kw, "boxplot_cindex_kw_no_lineal.png"),
  plot_cindex_kw_no_lineal,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

