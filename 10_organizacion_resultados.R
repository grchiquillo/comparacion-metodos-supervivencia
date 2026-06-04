##########################################################################

# 10. organización de resultados

## 10.1 Verificación inicial


resultados_modelos_super %>%
  dplyr::summarise(
    n_filas = dplyr::n(),
    n_escenarios = dplyr::n_distinct(escenario),
    n_modelos = dplyr::n_distinct(modelo),
    n_repeticiones = dplyr::n_distinct(rep),
    n_ibs = sum(!is.na(ibs)),
    n_c_index = sum(!is.na(c_index)),
    incluye_censura_95 = any(censura_objetivo == 0.95, na.rm = TRUE),
    existe_error_estimacion = "error_estimacion" %in% names(resultados_modelos_super)
  )

stopifnot(!("error_estimacion" %in% names(resultados_modelos_super)))
stopifnot(!any(resultados_modelos_super$censura_objetivo == 0.95, na.rm = TRUE))


## 10.2 Funciones auxiliares


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

mediana_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

q_na <- function(x, prob) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs = prob, na.rm = TRUE))
}

colapsar_unicos <- function(x) {
  x <- unique(na.omit(as.character(x)))
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

primer_valor_no_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[1]
}


## 10.3 Tabla resumen por escenario y modelo


resumen_resultados_por_escenario_modelo <- resultados_modelos_super %>%
  dplyr::filter(
    modelo %in% c("CoxPH", "SRF", "GBS"),
    censura_objetivo %in% c(0.20, 0.50, 0.75)
  ) %>%
  dplyr::group_by(
    escenario,
    n,
    censura_objetivo,
    tipo_relacion,
    modelo
  ) %>%
  dplyr::summarise(
    n_repeticiones = dplyr::n(),
    n_repeticiones_ok = sum(estado_ajuste == "ok", na.rm = TRUE),
    n_ibs_validos = sum(!is.na(ibs) & is.finite(ibs)),
    n_c_index_validos = sum(!is.na(c_index) & is.finite(c_index)),
    
    ibs_media = media_na(ibs),
    ibs_sd = sd_na(ibs),
    ibs_mediana = mediana_na(ibs),
    ibs_p25 = q_na(ibs, 0.25),
    ibs_p75 = q_na(ibs, 0.75),
    ibs_min = ifelse(all(is.na(ibs)), NA_real_, min(ibs, na.rm = TRUE)),
    ibs_max = ifelse(all(is.na(ibs)), NA_real_, max(ibs, na.rm = TRUE)),
    
    c_index_media = media_na(c_index),
    c_index_sd = sd_na(c_index),
    c_index_mediana = mediana_na(c_index),
    c_index_p25 = q_na(c_index, 0.25),
    c_index_p75 = q_na(c_index, 0.75),
    c_index_min = ifelse(all(is.na(c_index)), NA_real_, min(c_index, na.rm = TRUE)),
    c_index_max = ifelse(all(is.na(c_index)), NA_real_, max(c_index, na.rm = TRUE)),
    
    tiempo_total_segundos = sum(tiempo_ajuste_validacion_segundos, na.rm = TRUE),
    tiempo_medio_segundos = media_na(tiempo_ajuste_validacion_segundos),
    
    metodo_ibs = colapsar_unicos(metodo_ibs),
    estados_ajuste = colapsar_unicos(estado_ajuste),
    mensajes_ibs = colapsar_unicos(mensaje_ibs),
    
    mtry = primer_valor_no_na(mtry),
    nodesize = primer_valor_no_na(nodesize),
    ntree = primer_valor_no_na(ntree),
    splitrule = primer_valor_no_na(splitrule),
    
    eta = primer_valor_no_na(eta),
    max_depth = primer_valor_no_na(max_depth),
    nrounds = primer_valor_no_na(nrounds),
    subsample = primer_valor_no_na(subsample),
    colsample_bytree = primer_valor_no_na(colsample_bytree),
    best_nrounds = primer_valor_no_na(best_nrounds),
    nrounds_usado = primer_valor_no_na(nrounds_usado),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    modelo = factor(
      modelo,
      levels = c("CoxPH", "SRF", "GBS")
    )
  ) %>%
  dplyr::arrange(
    tipo_relacion,
    censura_objetivo,
    n,
    modelo
  )

resumen_resultados_por_escenario_modelo %>%
  print(n = Inf, width = Inf)


## 10.4 Ranking de modelos por escenario


resumen_resultados_por_escenario_modelo <- resumen_resultados_por_escenario_modelo %>%
  dplyr::group_by(
    escenario,
    n,
    censura_objetivo,
    tipo_relacion
  ) %>%
  dplyr::mutate(
    ranking_ibs = dplyr::min_rank(ibs_media),
    ranking_c_index = dplyr::min_rank(dplyr::desc(c_index_media)),
    mejor_ibs = ranking_ibs == 1,
    mejor_c_index = ranking_c_index == 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    tipo_relacion,
    censura_objetivo,
    n,
    modelo
  )


## 10.5 Tabla fila por escenario


resumen_resultados_por_escenario_wide <- resumen_resultados_por_escenario_modelo %>%
  dplyr::select(
    escenario,
    n,
    censura_objetivo,
    tipo_relacion,
    modelo,
    n_repeticiones_ok,
    n_ibs_validos,
    n_c_index_validos,
    ibs_media,
    ibs_sd,
    ibs_mediana,
    c_index_media,
    c_index_sd,
    c_index_mediana,
    ranking_ibs,
    ranking_c_index,
    mejor_ibs,
    mejor_c_index
  ) %>%
  tidyr::pivot_wider(
    names_from = modelo,
    values_from = c(
      n_repeticiones_ok,
      n_ibs_validos,
      n_c_index_validos,
      ibs_media,
      ibs_sd,
      ibs_mediana,
      c_index_media,
      c_index_sd,
      c_index_mediana,
      ranking_ibs,
      ranking_c_index,
      mejor_ibs,
      mejor_c_index
    ),
    names_glue = "{.value}_{modelo}"
  ) %>%
  dplyr::arrange(
    tipo_relacion,
    censura_objetivo,
    n
  )

resumen_resultados_por_escenario_wide %>%
  print(n = Inf, width = Inf)


## 10.6 Identificar mejor modelo por escenario

mejor_modelo_por_escenario <- resumen_resultados_por_escenario_modelo %>%
  dplyr::group_by(
    escenario,
    n,
    censura_objetivo,
    tipo_relacion
  ) %>%
  dplyr::summarise(
    mejor_modelo_ibs = paste(
      as.character(modelo[ranking_ibs == 1]),
      collapse = " | "
    ),
    mejor_ibs_media = min(ibs_media, na.rm = TRUE),
    
    mejor_modelo_c_index = paste(
      as.character(modelo[ranking_c_index == 1]),
      collapse = " | "
    ),
    mejor_c_index_media = max(c_index_media, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    tipo_relacion,
    censura_objetivo,
    n
  )

mejor_modelo_por_escenario %>%
  print(n = Inf, width = Inf)



## 10.7 Guardar tablas resumen


carpeta_resumen_escenarios <- file.path(
  ruta_super,
  "resumen_resultados_24_ibs_pec"
)

dir.create(
  carpeta_resumen_escenarios,
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  resumen_resultados_por_escenario_modelo,
  file.path(
    carpeta_resumen_escenarios,
    "resumen_resultados_por_escenario_modelo.csv"
  ),
  row.names = FALSE
)

write.csv(
  resumen_resultados_por_escenario_wide,
  file.path(
    carpeta_resumen_escenarios,
    "resumen_resultados_por_escenario_wide.csv"
  ),
  row.names = FALSE
)

write.csv(
  mejor_modelo_por_escenario,
  file.path(
    carpeta_resumen_escenarios,
    "mejor_modelo_por_escenario.csv"
  ),
  row.names = FALSE
)

saveRDS(
  resumen_resultados_por_escenario_modelo,
  file.path(
    carpeta_resumen_escenarios,
    "resumen_resultados_por_escenario_modelo.rds"
  )
)

saveRDS(
  resumen_resultados_por_escenario_wide,
  file.path(
    carpeta_resumen_escenarios,
    "resumen_resultados_por_escenario_wide.rds"
  )
)

saveRDS(
  mejor_modelo_por_escenario,
  file.path(
    carpeta_resumen_escenarios,
    "mejor_modelo_por_escenario.rds"
  )
)

