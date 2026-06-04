##############################################################################

# 7. Simulación de los 24 escenarios

escenarios <- expand.grid(
  n = c(100, 500, 1000, 5000),
  censura_objetivo = c(0.20, 0.50, 0.75),
  tipo_relacion = c("lineal", "no_lineal"),
  rep = 1:200,
  stringsAsFactors = FALSE
)

beta0_usar <- beta0_calibrado

lista_escenarios <- vector("list", nrow(escenarios))

for (i in seq_len(nrow(escenarios))) {
  
  esc <- escenarios[i, ]
  
  sim_i <- simular_dataset_supervivencia(
    n = esc$n,
    seed = 1000 + i,
    tipo_relacion = esc$tipo_relacion,
    censura_objetivo = esc$censura_objetivo,
    gamma_evento = 1,
    beta0 = beta0_usar,
    gamma_cens = 1,
    admin_time = 120
  )
  
  attr(sim_i, "escenario_id") <- i
  attr(sim_i, "rep") <- esc$rep
  
  lista_escenarios[[i]] <- sim_i
  
  if (i %% 100 == 0) {
    cat("Simulados:", i, "de", nrow(escenarios), "\n")
  }
}

cat("Total de bases simuladas:", length(lista_escenarios), "\n")



## 7.1 Guardar los modelos simulados
#Se guardan los 4800 escenarios simulados (24*200 repeticiones)

dir.create("escenarios_simulados", showWarnings = FALSE)

for (i in seq_along(lista_escenarios)) {
  
  esc <- escenarios[i, ]
  
  nombre <- paste0(
    "escenarios_simulados/",
    "esc_n", esc$n,
    "_c", sprintf("%03d", as.integer(esc$censura_objetivo * 100)),
    "_", esc$tipo_relacion,
    "_rep", sprintf("%03d", esc$rep),
    ".rds"
  )
  
  saveRDS(lista_escenarios[[i]], file = nombre)
  
  if (i %% 500 == 0) {
    cat("Guardados:", i, "de", length(lista_escenarios), "\n")
  }
}

saveRDS(escenarios, "escenarios_simulados/tabla_escenarios.rds")

cat("Archivos guardados en carpeta escenarios_simulados\n")


### 7.1.1 Reagrupación de archivos

ruta_escenarios <- "escenarios_simulados"

archivos_rds <- list.files(
  path = ruta_escenarios,
  pattern = "\\.rds$",
  full.names = TRUE
)

# Excluir archivos auxiliares
archivos_rds <- archivos_rds[
  !grepl("tabla_escenarios|resumen_verificacion|escenarios_24_listas", archivos_rds)
]



### 7.1.2 Función para extarer la información del nombre del archivo:

extraer_info_archivo <- function(archivo) {
  
  nombre <- basename(archivo)
  
  n_val <- str_extract(nombre, "(?<=esc_n)\\d+") |> as.numeric()
  c_txt <- str_extract(nombre, "(?<=_c)\\d+")
  rep_val <- str_extract(nombre, "(?<=_rep)\\d+") |> as.integer()
  
  tipo_val <- ifelse(
    str_detect(nombre, "_no_lineal_"),
    "no_lineal",
    "lineal"
  )
  
  censura_val <- case_when(
    c_txt %in% c("020", "02") ~ 0.20,
    c_txt %in% c("050", "05", "0500") ~ 0.50,
    c_txt == "075" ~ 0.75,
    TRUE ~ NA_real_
  )
  
  tibble(
    archivo = archivo,
    nombre_archivo = nombre,
    n = n_val,
    censura_objetivo = censura_val,
    tipo_relacion = tipo_val,
    rep = rep_val
  )
}


### 7.1.3 Tabla con la metadata de los 4800 archivos

tabla_archivos <- bind_rows(lapply(archivos_rds, extraer_info_archivo)) %>%
  filter(
    !is.na(censura_objetivo),
    censura_objetivo %in% c(0.20, 0.50, 0.75)
  ) %>%
  arrange(n, censura_objetivo, tipo_relacion, rep)

dim(tabla_archivos)
head(tabla_archivos)

tabla_archivos %>%
  count(censura_objetivo)


### 7.1.4 Crear la tabla de los 24 escenarios únicos

tabla_24_escenarios <- expand.grid(
  n = c(100, 500, 1000, 5000),
  censura_objetivo = c(0.20, 0.50, 0.75),
  tipo_relacion = c("lineal", "no_lineal"),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    escenario_id = row_number(),
    escenario_nombre = paste0(
      "esc_n", n,
      "_c", sprintf("%03d", as.integer(censura_objetivo * 100)),
      "_", tipo_relacion
    )
  )

tabla_24_escenarios
nrow(tabla_24_escenarios)



### 7.1.5 Agrupar las bases simuladas en una lista por escenario

escenarios_24_listas <- vector("list", length = nrow(tabla_24_escenarios))
names(escenarios_24_listas) <- tabla_24_escenarios$escenario_nombre

for (i in seq_len(nrow(tabla_24_escenarios))) {
  
  esc_i <- tabla_24_escenarios[i, ]
  
  archivos_i <- tabla_archivos %>%
    filter(
      n == esc_i$n,
      censura_objetivo == esc_i$censura_objetivo,
      tipo_relacion == esc_i$tipo_relacion
    ) %>%
    arrange(rep)
  
  lista_i <- lapply(archivos_i$archivo, readRDS)
  
  names(lista_i) <- paste0("rep_", archivos_i$rep)
  
  attr(lista_i, "n") <- esc_i$n
  attr(lista_i, "censura_objetivo") <- esc_i$censura_objetivo
  attr(lista_i, "tipo_relacion") <- esc_i$tipo_relacion
  attr(lista_i, "escenario_id") <- esc_i$escenario_id
  
  escenarios_24_listas[[i]] <- lista_i
}


### 7.1.6 Guardar escenarios agrupados

dir.create("escenarios_agrupados", showWarnings = FALSE)

for (i in seq_along(escenarios_24_listas)) {
  
  saveRDS(
    escenarios_24_listas[[i]],
    file = file.path(
      "escenarios_agrupados",
      paste0(names(escenarios_24_listas)[i], ".rds")
    )
  )
}

saveRDS(
  escenarios_24_listas,
  file = file.path("escenarios_agrupados", "escenarios_24_listas.rds")
)

write.csv(
  tabla_24_escenarios,
  file = file.path("escenarios_agrupados", "tabla_24_escenarios.csv"),
  row.names = FALSE
)


### 7.1.7 Ejemplo de la primera repetición del primer escenario que
#es n=100, cencura=0.2 y lineal

#Esto nos sirve para ver si los escenarios simulados verdaderamente 
#muestran esas caracteristicas y si los valores de las covariables se 
#basan en la base original.


# Primer escenario
escenario_1 <- escenarios_24_listas[[1]]

# Primera repetición del primer escenario
train_data_ejm <- escenario_1[[1]]

dim(train_data_ejm)
names(train_data_ejm)
head(train_data_ejm)

summary(train_data_ejm )



