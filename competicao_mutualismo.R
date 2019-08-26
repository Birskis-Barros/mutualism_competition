getwd()

#necessary packages  

library(bipartite)
library("ggplot2")
library("GGally")
library("sna")
library(igraph)
library("network")


#Parameters 
sp_p = 5 #number of plants
sp_a = 4 #number of animals 
N = sp_a + sp_p #total richness
alpha = -2 

#Initial Trait 
z_a = runif(sp_a, 0, 10)# for the animals 
z_p = runif(sp_p, 0, 10)# for the plants 

#data <- matrix(NA, ncol = 2, nrow = t_max)

t_max = 1000
#for(i in 1:t_max){
  
### Calculating trait matching among plants and animals (alfa must be low)
  z_dif = sapply(z_a, "-", z_p) # difference between traits 
  A <- matrix(1, nrow=sp_p, ncol=sp_a) #total connected matrix 
  
  # In A, if the value of the difference between traits is more than 3 (arbitrary), 
  #we remove the interaction (=0)
  A[abs(z_dif)> 3] = 0 
  
  # Calculating competition among animals
  mat.nij <- matrix(NA, ncol=sp_a, nrow=sp_a) #matrix with only animals
  
  #for how much overleap exist (related to the above part in eq 1)
  #(quanto de sobreposicao nas interacoes existem (parte superior da eq 1))
  mat.int <- t(A) %*% A  
  
  tot <- diag(mat.int)
  
  for (i in 1:nrow(mat.nij))
  {
    for(j in 1:ncol(mat.nij)){
     mat.nij[i,j] <- mat.int[i,j]/(tot[i]+tot[j]-mat.int[i,j]) 
    }
  }
  
  #quanto maior o valor de z_dif_animal, maior a similaridade entre 
  #as especies e, portando, maior a competicao
  
  z_dif_animal = sapply(z_a, "-", z_a)^2 #(equacao 2 - sem a exponencial)
  
  matriz.c <-   mat.nij*exp(alpha*z_dif_animal) #(equacao 3)
  
  diag(matriz.c) <- NA
  
  1- apply(matriz.c,2,mean, na.rm=TRUE)
  #o primeiro valor daqui vai multiplicar com todos os elementos da primeira coluna
  # o segundo valor vai multiplicar com todos os elementos da segunda coluna.. e assim sucessivamente 
  
  #depois de calcular isso, precisamos criar uma matriz igual a A, mas dadas as probabilidades calculadas 
  
  

 