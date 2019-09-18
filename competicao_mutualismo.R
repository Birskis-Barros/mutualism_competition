#necessary packages  

library("bipartite")
library("ggplot2")
library("GGally")
library("sna")
library("igraph")
library("network")
require("plyr")
require("ggnetwork")
require("dplyr")

#Parameters 
sp_p = 10 #number of plants
sp_a = 20 #number of animals 
N = sp_a + sp_p #total richness
alpha = -2 

#Initial Trait 
z_a = runif(sp_a, 0, 10)# for the animals 
z_p = runif(sp_p, 0, 10)# for the plants 

#data <- matrix(NA, ncol = 2, nrow = t_max)

t_max = 200
redes = list()

### Calculating trait matching among plants and animals (alfa must be low)
z_dif = sapply(z_a, "-", z_p) # difference between traits 
A <- matrix(1, nrow=sp_p, ncol=sp_a) #total connected matrix 

# In A, if the value of the difference between traits is more than 3 (arbitrary), 
#we remove the interaction (=0)
A[abs(z_dif)> 3] = 0 
A.ini <- A

for(n in 1:t_max){
  
  # Calculating competition among animals
  mat.nij <- matrix(NA, ncol=sp_a, nrow=sp_a) #matrix with only animals
  
  #for how much overleap exist (related to the above part in eq 1)
  #(quanto de sobreposicao nas interacoes existem (parte superior da eq 1))
  mat.int <- t(A) %*% A  
  
  tot <- diag(mat.int) #o numero total de interações que a espécie da coluna tem 
  
  for (i in 1:nrow(mat.nij))
  {
    for(j in 1:ncol(mat.nij)){
     mat.nij[i,j] <- mat.int[i,j]/(tot[i]+tot[j]-mat.int[i,j]) 
    }
  }
  
  #mat.int[i,j]  = numero de interacoes que a especie da coluna compartilha com a da linha 
  
  #quanto maior o valor de z_dif_animal, maior a similaridade entre 
  #as especies e, portando, maior a competicao
  
  z_dif_animal = sapply(z_a, "-", z_a)^2 #(equacao 2 - sem a exponencial)
  
  matriz.c <-   mat.nij*exp(alpha*z_dif_animal) #(equacao 3)
  diag(matriz.c) <- NA #removendo a competicao intraespecifica 
 
  #Calculating the probability to maintain the interaction 

  manter.interacao <- 1- apply(matriz.c,2,mean, na.rm=TRUE)  #média de competicao pra cada especie de coluna com a especie da linha
  
  manter.interacao <- matrix(manter.interacao, ncol=sp_a, nrow=sp_p, byrow=TRUE)
  
  teste <- matrix(runif(20,0,1), ncol=sp_a, nrow=sp_p, byrow=TRUE)
  
  A <- ifelse(teste>manter.interacao,0,A)
  
 redes[[n]] = A
  
  }

###Analyzing 

ggnet2(network(redes[[10]]), label=TRUE)

plot(x=1:50, y=plyr::laply(redes, function(x) sum(x)))
ggnet2(network(redes[[2]]), label=TRUE)


length(which(apply(redes[[10]], 2, sum)==0))

llply(redes, function(x) length(which(apply(x, 2, sum)==0)))
llply(redes, function(x) length(which(apply(x, 1, sum)==0)))

llply(redes, function(x) networklevel(x, index="NODF" ))

llply(redes, function(x) networklevel(x, index="connectance" ))

redes_analise = data.frame()
for ( i in 1:t_max){
  redes_analise[i,1] = llply(redes, function(x) networklevel(x, index="connectance"))[[i]]
  redes_analise[i,2] = llply(redes, function(x) networklevel(x, index="NODF"))[[i]]
}

colnames(redes_analise) = c("connectance", "aninhamento", "t_max")
redes_analise$t_max = 1:t_max
plot(redes_analise$connectance~redes_analise$t_max)
plot(redes_analise$aninhamento~redes_analise$t_max)

