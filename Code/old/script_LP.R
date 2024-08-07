library(lpSolve)
library(readxl)
library(writexl)
library(ggplot2)
library(viridis)
library(scales)
library(corrplot)
library(sf)
library(dplyr)

choose.files()

{

parques <- read_xlsx("C:\\Users\\gtarr\\OneDrive\\Mestrado\\GIS_Data\\Data\\Apt_comite\\parques_ex_comite.xlsx")

# Par�metros do problema
orcamento_total <- 100000000
num_min <- 20

# N�mero de parques
num_parques <- nrow(parques)

# Matriz de custos
custos <- parques$custo

# Coeficientes da fun��o objetiva
obj_coef <- parques$nor_lst # maximizar rela��o oferta do CRS

# Dire��o do problema (maximizar ou minimizar)
direction <- "max"

# Restri��o de custo total
A <- matrix(custos, nrow = 1)
dir <- "<="
rhs <- orcamento_total

# Restri��o n� minimo
B <- matrix(rep(1, num_parques), nrow = 1)
dir_b <- ">="
rhs_b <- num_min

# restri��o de selecionar pelo menos um parque
A_eq <- matrix(rep(1, num_parques), nrow = 1)
dir_eq <- ">="
rhs_eq <- 1

# Adicionando vari�vel bin�ria para cada parque (0 quando n�o selecionado e 1 quando selecionado)
A_bin <- diag(num_parques)

# Combinando todas as matrizes
A_final <- rbind(A,B,A_eq, A_bin)
dir_final <- c(rep(dir, 1), dir_b, dir_eq,rep("<=", num_parques))
rhs_final <- c(rhs,rhs_b, rhs_eq, rep(1, num_parques))

# Resolva o problema de programa��o linear
sol <- lp(direction, objective.in = obj_coef, const.mat = A_final, const.dir = dir_final, const.rhs = rhs_final)
print(sol$status)
best_sol <- round(sol$solution, digits = 3)

# Parques selecionados
parques_selecionados <- parques[best_sol == 1, ]

# Confirmando se restri��es foram respeitadas
sum(parques_selecionados$custo) < orcamento_total

# Parques n�o selecionados

parques_nao_selecionados <- parques[best_sol < 1, ]

}
