# Comparación de modelos de supervivencia: CoxPH, SRF y GBS

Este repositorio contiene el código en R desarrollado para el trabajo de grado de la Maestría en Bioestadística. El objetivo del estudio es comparar el desempeño predictivo de tres modelos de análisis de supervivencia bajo diferentes escenarios simulados: el modelo de riesgos proporcionales de Cox (CoxPH), los bosques aleatorios de supervivencia (Survival Random Forest, SRF) y el modelo de potenciación por gradiente para supervivencia (Gradient Boosting Survival, GBS).

La estrategia metodológica consiste en generar datos simulados a partir de una estructura derivada de datos reales y evaluar el desempeño de los modelos en escenarios controlados. La comparación se realiza considerando dos dimensiones principales: la calibración, medida mediante el Integrated Brier Score (IBS), y la discriminación, evaluada a través del C-index.

Los datos reales provienen de la base METABRIC y se utilizan únicamente como referencia para estimar la estructura estadística de las covariables y definir los parámetros de la simulación. A partir de esta información se generan las bases simuladas empleadas en el análisis comparativo. De esta forma, el estudio permite evaluar el comportamiento relativo de los modelos bajo condiciones previamente definidas de tamaño de muestra, porcentaje de censura y tipo de relación entre covariables y tiempo de supervivencia (lineal o no lineal).

## Datos

El código parte del archivo `METABRIC_RNA_Mutation.csv`, del cual se seleccionan
nueve variables:

- `age_at_diagnosis`
- `type_of_breast_surgery`
- `chemotherapy`
- `neoplasm_histologic_grade`
- `lymph_nodes_examined_positive`
- `her2_status`
- `tumor_size`
- `overall_survival_months`
- `overall_survival`

El archivo de datos no se incluye en el repositorio. Debe ubicarse en el
directorio de trabajo antes de ejecutar el primer script, ya que la carga se
hace con una ruta relativa.

## Estructura del repositorio

El código original se encuentra subdividido en diez scripts de R, uno por cada
sección del análisis, con el fin de facilitar su reproducibilidad. El orden de
los archivos corresponde al orden de ejecución.

```
.
├── 01_seleccion_variables.R
├── 02_descripcion_general.R
├── 03_simulacion_covariables.R
├── 04_simulacion_tiempos_supervivencia.R
├── 05_simulacion_tiempos_censura.R
├── 06_simulacion_completa.R
├── 07_simulacion_24_escenarios.R
├── 08_validacion_superpoblacion.R
├── 09_ajuste_validacion_modelos.R
├── 10_organizacion_resultados.R
├── METABRIC_RNA_Mutation.csv        (no incluido; debe agregarse)
├── escenarios_simulados/            (generada al ejecutar)
├── escenarios_agrupados/            (generada al ejecutar)
└── prueba_superpoblacion_n10000/    (generada al ejecutar)
```

## Descripción de los scripts

A continuación se describe el contenido de cada script, conservando la
numeración de las secciones del código original.

1. **Selección de variables.** Se cargan las librerías y la base METABRIC y se
   seleccionan las nueve covariables de interés.

2. **Descripción general.** Se revisan los valores faltantes y se imputan, se
   examina la distribución de cada variable, se evalúa la normalidad de las
   continuas mediante el test de Shapiro-Wilk, se convierten los tipos de datos
   y se construye la tabla resumen de las covariables en la base real.

3. **Simulación de covariables.** Se simulan las variables en cadena a partir de
   modelos ajustados sobre los datos reales: normal truncada para la edad,
   regresión binomial negativa para el número de nodos positivos, regresión
   logística para `her2_status`, `chemotherapy` y `type_of_breast_surgery`,
   regresión lineal para el tamaño del tumor y modelo multinomial para el grado
   histológico. Se incluye la validación de cada modelo y las funciones
   auxiliares que generan las variables.

4. **Simulación de los tiempos de supervivencia.** Se estandarizan las
   covariables continuas y se simula el tiempo al evento bajo una distribución
   Weibull. La escala basal se calibra con la curva Kaplan-Meier de la base
   METABRIC y el riesgo se resume en un predictor lineal que admite relaciones
   lineales o no lineales con las covariables.

5. **Simulación de los tiempos de censura.** Se incorpora la censura mediante una
   variable Bernoulli y el tiempo observado se construye como el mínimo entre el
   tiempo al evento y el tiempo de censura. Se incluye una función que calibra la
   escala de censura para alcanzar las proporciones objetivo.

6. **Simulación completa.** Se integran covariables, tiempo al evento y censura
   en una sola rutina y se realiza una prueba de la simulación.

7. **Simulación de los 24 escenarios.** Los escenarios resultan de combinar
   cuatro tamaños de muestra (100, 500, 1000 y 5000), tres niveles de censura
   objetivo (20%, 50% y 75%) y dos tipos de relación (lineal y no lineal), lo que
   da 24 escenarios. Con 200 repeticiones por escenario se obtienen 4800 bases
   simuladas, que se guardan en disco.

8. **Validación en superpoblación de n=10.000.** Se genera una superpoblación de
   diez mil observaciones por escenario, se definen las funciones de ajuste de
   CoxPH, SRF y GBS y se optimizan los hiperparámetros de SRF y GBS mediante
   validación cruzada interna. Los hiperparámetros seleccionados se almacenan en
   disco.

9. **Ajuste y validación de modelos.** Cada modelo se ajusta sobre las
   repeticiones simuladas y se valida contra la superpoblación externa del
   escenario correspondiente. Como métricas de desempeño se emplean el índice de
   concordancia (C-index) y el Integrated Brier Score (IBS) calculado con `pec`.
   Se consolidan los resultados, se construye la tabla resumen por escenario, se
   generan los diagramas de caja y se comparan los modelos mediante la prueba de
   Kruskal-Wallis.

10. **Organización de resultados.** Se verifican y resumen los resultados por
    escenario y modelo, se construye el ranking de modelos, se identifica el
    mejor modelo en cada escenario y se guardan las tablas finales.

## Requisitos

El código se ejecuta en R. Se requieren los siguientes paquetes, que deben
instalarse previamente:

`skimr`, `tidyverse`, `dplyr`, `purrr`, `ggcorrplot`, `nnet`, `truncnorm`,
`MASS`, `VIM`, `knitr`, `tibble`, `AER`, `pROC`, `ResourceSelection`, `car`,
`survival`, `stringr`, `pec`, `prodlim`, `randomForestSRC`, `xgboost`,
`ggplot2`, `forcats` y `readr`.

La instalación puede hacerse con:

```r
install.packages(c(
  "skimr", "tidyverse", "ggcorrplot", "nnet", "truncnorm", "MASS", "VIM",
  "knitr", "AER", "pROC", "ResourceSelection", "car", "survival", "pec",
  "prodlim", "randomForestSRC", "xgboost"
))
```

## Ejecución

Los diez scripts forman un único flujo secuencial y comparten objetos en memoria
(rutas, tablas de escenarios y funciones de ajuste), además de los resultados
intermedios que se guardan en disco. Por lo tanto deben ejecutarse en orden, del
`01` al `10`, dentro de una misma sesión de R; no es posible ejecutar un script
de forma aislada sin haber corrido los anteriores.

Antes de comenzar conviene tener en cuenta tres puntos:

- El directorio de trabajo debe ser la carpeta que contiene
  `METABRIC_RNA_Mutation.csv`. Se recomienda trabajar dentro de un proyecto de
  RStudio (`.Rproj`).
- Los archivos están codificados en UTF-8. Si se ejecutan con `source()`, debe
  indicarse `encoding = "UTF-8"`.
- Los scripts `08` y `09` son los de mayor costo computacional, ya que incluyen
  la optimización de hiperparámetros y el ajuste de los tres modelos sobre todos
  los escenarios. Su primera ejecución puede tardar de forma considerable.

La ejecución manual consiste en abrir cada script y correrlo completo de forma
ordenada. De manera equivalente, los scripts pueden ejecutarse en cadena con un
script maestro:

```r
setwd("ruta/al/proyecto")   # carpeta que contiene el CSV

archivos <- c(
  "01_seleccion_variables.R",
  "02_descripcion_general.R",
  "03_simulacion_covariables.R",
  "04_simulacion_tiempos_supervivencia.R",
  "05_simulacion_tiempos_censura.R",
  "06_simulacion_completa.R",
  "07_simulacion_24_escenarios.R",
  "08_validacion_superpoblacion.R",
  "09_ajuste_validacion_modelos.R",
  "10_organizacion_resultados.R"
)

for (f in archivos) {
  message(">>> Ejecutando ", f)
  source(f, encoding = "UTF-8")
}
```

## Resultados generados

Durante la ejecución se crean varias carpetas con los productos del análisis. En
`escenarios_simulados/` y `escenarios_agrupados/` se almacenan las bases
simuladas y la tabla de escenarios. En `prueba_superpoblacion_n10000/` se
guardan las superpoblaciones, los hiperparámetros seleccionados para SRF y GBS,
los resultados de la validación externa, las tablas consolidadas y las figuras.
Los resultados se conservan en formato `.rds` y `.csv`, de modo que las
ejecuciones posteriores reutilizan lo ya calculado.
