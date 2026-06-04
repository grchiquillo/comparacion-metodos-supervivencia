#################################################################

# 2. Descripción general

skim(datos1) %>% print(n = Inf)

## 2.1 Revisión valores faltantes

### 2.1.1 Valores faltantes por variable
na_por_columna <- colSums(is.na(datos1))
porcentaje_na <- round((na_por_columna / nrow(datos1)) * 100, 2)
resultado_na <- data.frame(
  Variable = names(na_por_columna),
  Conteo_NA = na_por_columna,
  Porcentaje_NA = porcentaje_na,
  stringsAsFactors = FALSE
)

resultado_na

### 2.1.2 Filas con valores faltantes

sum(is.na(datos1))

### 2.1.3 Imputación de datos

datos1 <- kNN(datos1, k = 5, imp_var = FALSE)

sum(is.na(datos1))


## 2.2 Descripción de cada variable

#Se hace una revisión visual de cada una de las variables para conocer su distribución,
#identificar si hay valores atípicos, propoción de valores por categoría, etc.


### 2.2.1 Variables continuas

#age_at_diagnosis
summary(datos1$age_at_diagnosis)

ggplot(datos1, aes(x = age_at_diagnosis)) +
  geom_histogram(binwidth = 5, fill = "#FFA07A", color = "black", alpha=0.7) +
  labs(title = "Distribución de la edad de los pacientes",
       x = "Edad (años)", y = "Frecuencia") +
  theme_minimal()


#tumor_size
summary(datos1$tumor_size)

ggplot(datos1, aes(x = tumor_size)) +
  geom_histogram(binwidth = 5, fill = "#FFA07A", color = "black", alpha=0.7) +
  labs(title = "Distribución del tamaño del tumor",
       x = "Tamaño tumor", y = "Frecuencia") +
  theme_minimal()


#lymph_nodes_examined_positive
summary(datos1$lymph_nodes_examined_positive)

ggplot(datos1, aes(x = lymph_nodes_examined_positive)) +
  geom_histogram(binwidth = 5, fill = "#FFA07A", color = "black", alpha=0.7) +
  labs(title = "Distribución del numero nodos examinados positivos",
       x = "Numero nodos examinados positivos", y = "Frecuencia") +
  theme_minimal()


#overall_survival_months
summary(datos1$overall_survival_months)

ggplot(datos1, aes(x = overall_survival_months)) +
  geom_histogram(binwidth = 5, fill = "#FFA07A", color = "black", alpha=0.7) +
  labs(title = "Distribución del tiempo de supervivencia",
       x = "Tiempo supervivencia", y = "Frecuencia") +
  theme_minimal()


#Se evalúa la normalidad en las variables continuas

# Test de Shapiro-Wilk variable edad
shapiro.test(datos1$age_at_diagnosis)

# Test de Shapiro-Wilk variable tamaño tumor
shapiro.test(datos1$tumor_size)

# Test de Shapiro-Wilk variable tamaño tumor
shapiro.test(datos1$lymph_nodes_examined_positive)


### 2.2.2 Varaibles categóricas

# type_of_breast_surgery
tablatipocirugia<-table(datos1$type_of_breast_surgery)
tablatipocirugia

bp2<-barplot(tablatipocirugia,                       
             main="Distribución por tipo de cirugía",               
             xlab="Tipo cirugía",                   
             ylab="Numero pacientes",                   
             col=c("#DDA0DD","#87CEEB"))

text(x = bp2, 
     y = tablatipocirugia, 
     label = tablatipocirugia, 
     pos = 1,  
     cex = 0.8, 
     col = "black")



#her2_status
tablaher2<-table(datos1$her2_status)
tablaher2

bp2<-barplot(tablaher2,                       
             main="Distribución por her2 status",               
             xlab="Estado Her2",                   
             ylab="Numero pacientes",                   
             col=c("#DDA0DD","#87CEEB"))

text(x = bp2, 
     y = tablaher2, 
     label = tablaher2, 
     pos = 1,  
     cex = 0.8, 
     col = "black")


#chemotherapy
tablaquim<-table(datos1$chemotherapy)
tablaquim

bp2<-barplot(tablaquim,                       
             main="Distribución por quimioterapia",               
             xlab="Quimioterapia",                   
             ylab="Numero pacientes",                   
             col=c("#DDA0DD","#87CEEB"))

text(x = bp2, 
     y = tablaquim, 
     label = tablaquim, 
     pos = 1,  
     cex = 0.8, 
     col = "black")


#neoplasm_histologic_grade
tablagradohis<-table(datos1$neoplasm_histologic_grade)
tablagradohis

bp2<-barplot(tablagradohis,                       
             main="Distribución por grado histologico",               
             xlab="Grado histológic",                   
             ylab="Numero pacientes",                   
             col=c("#DDA0DD","#87CEEB", "#66CDAA"))

text(x = bp2, 
     y = tablagradohis, 
     label = tablagradohis, 
     pos = 1,  
     cex = 0.8, 
     col = "black")



#overall_survival
tablasuperv<-table(datos1$overall_survival)
tablasuperv

bp2<-barplot(tablasuperv,                       
             main="Distribución por supervivencia",               
             xlab="Supervivencia",                   
             ylab="Numero pacientes",                   
             col=c("#DDA0DD","#87CEEB"))

text(x = bp2, 
     y = tablasuperv, 
     label = tablasuperv, 
     pos = 1,  
     cex = 0.8, 
     col = "black")



## 2.3 Verificación y conversión de tipos de datos

#Para facilitar la simulación de la base de datos a partir de estos datos reales
#se procede a volver factor las variables categóricas, a renombrar sus categorias
#y a verificar que las variables continuas si estén como numericas.


### 2.3.1 Volver factor variables categóricas y renombrar categorias

datos1$her2_status<- factor(datos1$her2_status, 
                            labels = c("Negativo", "Positivo"))

datos1$neoplasm_histologic_grade<- factor(datos1$neoplasm_histologic_grade, 
                                          labels = c("Grado 1", "Grado 2", "Grado 3"))

datos1$chemotherapy<- factor(datos1$chemotherapy, 
                             labels = c("No", "Si"))


datos1$type_of_breast_surgery<- factor(datos1$type_of_breast_surgery,
                                       labels = c("Conserva_la_mama", "Mastectomia"))




### 2.3.2 Verificar que las variables numericas esten como numericas

datos1$age_at_diagnosis <- as.numeric(datos1$age_at_diagnosis)
datos1$lymph_nodes_examined_positive  <- as.numeric(datos1$lymph_nodes_examined_positive)
datos1$tumor_size  <- as.numeric(datos1$tumor_size)
datos1$overall_survival_months  <- as.numeric(datos1$overall_survival_months)


## 2.4 Tabla resumen estadisticos

#Con el fin de justificar los parámetros empleados en la simulación, 
#se presenta una tabla resumen de las covariables seleccionadas en la base real. 
#Para las variables continuas se reportan media, desviación estándar, 
#mínimo y máximo; para las variables categóricas se presentan 
#frecuencias absolutas y relativas.


### 2.4.1 Definir variables según su naturaleza

vars_numericas <- c(
  "age_at_diagnosis",
  "lymph_nodes_examined_positive",
  "tumor_size",
  "overall_survival_months"
)

vars_categoricas <- c(
  "type_of_breast_surgery",
  "chemotherapy",
  "neoplasm_histologic_grade",
  "her2_status",
  "overall_survival"
)


### 2.4.2 Tabla variables numericas
tabla_numericas <- datos1 %>%
  dplyr::select(all_of(vars_numericas)) %>%
  summarise(
    across(
      everything(),
      list(
        n = ~sum(!is.na(.)),
        media = ~round(mean(., na.rm = TRUE), 2),
        sd = ~round(sd(., na.rm = TRUE), 2),
        mediana = ~round(median(., na.rm = TRUE), 2),
        q1 = ~round(quantile(., 0.25, na.rm = TRUE), 2),
        q3 = ~round(quantile(., 0.75, na.rm = TRUE), 2),
        min = ~round(min(., na.rm = TRUE), 2),
        max = ~round(max(., na.rm = TRUE), 2)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", ".value"),
    names_sep = "_(?=[^_]+$)"
  )

kable(
  tabla_numericas,
  caption = "Tabla resumen de variables numéricas"
)


### 2.4.3 Tabla variables categóricas

tabla_categoricas <- datos1 %>%
  dplyr::select(all_of(vars_categoricas)) %>%
  mutate(
    chemotherapy = as.factor(chemotherapy),
    her2_status = as.factor(her2_status),
    neoplasm_histologic_grade = as.factor(neoplasm_histologic_grade),
    overall_survival = as.factor(overall_survival),
    type_of_breast_surgery = as.factor(type_of_breast_surgery)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Categoria"
  ) %>%
  filter(!is.na(Categoria)) %>%
  count(Variable, Categoria, name = "Frecuencia") %>%
  group_by(Variable) %>%
  mutate(
    Porcentaje = round(100 * Frecuencia / sum(Frecuencia), 2)
  ) %>%
  ungroup()

kable(
  tabla_categoricas,
  caption = "Tabla resumen de variables categóricas"
)


