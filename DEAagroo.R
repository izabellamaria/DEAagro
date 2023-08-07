#### Remover dados existentes

rm(list=ls())

## Carregar Pacotes (é necessário instalar os pacotes previamente)

library(sfsmisc)	
library(AER)
library(FEAR)
library(boot) 

## A rotina utiliza os quatro pacotes. Favor atentar para possíveis licenças para utilização dos pacotes. 

#### Ler a base de dados

agro<-read.table("C:\Users\agro.txt", h=T)

## xi(input), yi(output), zi(background)

## y1= renda bruta total
## x1= custos com mão de obra
## x2= custos com energia e agua
## x3= custos com irrigação
## z1= idade
## z2= escolaridade
## z3= capital total empregado
## z4= Nº total de empregados
## z5= Índice de introdução de inovação
## z6= Índice de inovações realizadas
## z7= Quantidade de tecnologia agrícola
## z8= Treinamento
## z9= Índice de fonte de informação
## z10= cooperativa



###############################################
## Teste para detecção de outlier - leverage ##
###############################################



## Criar novas matrizes input e output

x=matrix(nrow=3,ncol=85)
x[1,]=agro$x1
x[2,]=agro$x2
x[3,]=agro$x3


y=matrix(nrow=2,ncol=85)
y[1,]=agro$y1
y[2,]=1

#######################################################################################################################


## Estimador Sampaio-Stosic #### retornos variáveis ####

## Como são poucos dados, podemos construir o leverage sem utilizar o artifício das 
## Bubbles de Sampaio-Stosic.

## Estimando os coeficientes de eficiência para as K dmu`s

deat <- dea(XOBS=x,YOBS=y, ORIENTATION=1)


## Primeiro desenvolve-se uma função para calcular a DEA original, com todas as informações

first <- function(x, y, rts, orientation){ 
  
  j <- seq(1:ncol(x)) 
  
  rts <- seq(1:3) 
  
  orientation <- c(1,2) 
  
  saida1 <- NULL
  
  dea_t <- dea(XOBS=x,YOBS = y, RTS = rts[1], ORIENTATION = orientation[1])
  
  saida1 <- matrix(c(saida1,dea_t), nrow=length(j), ncol=1)
  
}

dea_t <- first(x,y, rts=1, orientation=2)

dea_t

## A segunda parte é retirar as dmu`s uma a uma, formulando para tanto outra função.
## A função dea.out serve para identificar as observações outliers, retirando as dmu's uma a uma e estimando o dea. Essa
## função requer os mesmos parâmetros que a função dea, porém o primeiro parâmetro é o número de observações neste caso = 85 

dea.out <- function(j, x, y, rts, orientation){
  
  j <- seq(1:ncol(x)) 
  
  jn <- length(j)
  
  rts <- seq(1:3) 
  
  orientation <- c(1,2) 
  
  saida2 <- NULL
  
  xi <- 0
  yi <- 0
  
  for (i in 1:jn){
    xi <- x[,-i]
    yi <- y[,-i]
    dea <- dea(XOBS=xi,YOBS=yi, RTS = rts[1], ORIENTATION = orientation[1]) 
    
    saida2 <- c(saida2, dea)
  }
  
  saida2
  
}

dea <- dea.out(85,x,y, rts=1, orientation=1)

## Lembrar de transformar dea em uma matrix com 'jn-1' linhas e 'jn' colunas:

dea_j <- matrix(dea, nrow =84, ncol=85)

## Agora é preciso colocar a fórmula presente em Sampaio de Sousa e Stosic (2005, pp.162)
## Função construída Jackstrap (jacknife+bootstrap):
## Primeiro inseri-se os valores "0" nas diagonais [i,i] no lugar de cada observação a fim de calcular o leverage:

orientation <- c(1,2)
rts <- seq(1:3)
comp1 <- 0
comp2 <- 0
comp3 <- 0
j <- seq(1:ncol(x)) 
jn <- length(j)
jnm <- length(j)-1
completar<- matrix(0, nrow=jn, ncol=jn)
saida3 <- NULL


## Fazemos a primeira e a última coluna da matriz de eficiência com os 0's enxertados.

comp1 <- c(completar[1,1], dea_j[1:jn-1,1])
comp3 <- c(dea_j[1:jn-1,jn], completar[jn,jn])

dea_j1 <- 0
dea_j2 <- 0
dea_j3 <- 0

## o miolo da matrix enxertada:

comp2 <- for (i in 2:jnm){
  
  dea_j1 <- dea_j[1:(i-1),i]
  dea_j2 <- completar[i,i]
  dea_j3 <- dea_j[i:jnm,i]
  vetor <- as.vector(c(dea_j1, dea_j2, dea_j3))
  
  saida3 <- c(saida3,vetor)
}
comp2 <- saida3

dea_J <- matrix(c(comp1, comp2, comp3),nrow=jn, ncol=jn)


dea_J

## Sai a matriz com os índices de eficiência após a retirada de uma a uma das observações
## linhas = n observações, colunas = n observações.  


## Agora a função da fórmula leverage propriamente, sem os bubbles.

jackknife.lev <- function (DEA_t,DEA_j)
  
{
  
  K <- length(j)
  
  dea_t <- dea(XOBS=x,YOBS = y, RTS = rts[1], ORIENTATION = orientation[1])
  
  theta <- c(dea_t)
  
  saida4 <- NULL
  
  lev <- for (i in 1:jn){
    theta_e <- dea_J[,i]
    square <- (theta_e-theta)^2
    soma <- sum(square[-i])
    leverage <- sqrt(soma/(K-1))
    saida4 <- c(saida4,leverage)
  }
  
  nleverage <- round(saida4, digits = 5)
  
  histograma <- hist(nleverage)
  kernel_lev <- density(nleverage)
  grafo <- lines(kernel_lev$x, kernel_lev$y, col="blue", lwd=1.5)
  
  
  ## Quais são os outliers? Usamos a regra de Bolso em Sampaio-Stosic lev >= 0.02
  
  outliers <- nleverage >= 0.02
  posicao <- matrix(c(nleverage,j), ncol=2, byrow=FALSE) 
  
  ## Pedindo para listar apenas os outliers
  
  nlev <- nleverage[outliers]
  
  outliers
  
  c(nlev,posicao[,2][posicao[,1]>=0.02])
  
}

jackknife.lev(dea_t, dea_J)

## Caso sejam identificados outliers eles devem ser extraídos da amostra.


#######################################################################################################################


############################################################
## Primeiro Estágio (DEA com viés corrigido por Bootstrap ##
############################################################

#### Remover dados existentes

rm(list=ls())

## Carregar Pacotes (é necessário instalar os pacotes previamente)

library(sfsmisc)	
library(AER)
library(FEAR)
library(boot) 

## A rotina utiliza os quatro pacotes. Favor atentar para possíveis licenças para utilização dos pacotes. 

#### Ler a nova base de dados sem outliers

agro<-read.table("C:\Users\Izabella Maria\Desktop\agro.txt", h=T)

## Criar matrizes input e output

x=matrix(nrow=1,ncol=85)
x[1,]=agro$x1
x[2,]=agro$x2
x[3,]=agro$x3

y=matrix(nrow=2,ncol=85)
y[1,]=agro$y1
y[2,]=1


## Estimando a eficiência (Shepard) Retornos constantes e Orientação output

dhat=dea(XOBS=x,YOBS=y,RTS=1,ORIENTATION=1)

## Estimando e corrigindo o viés , a variância e o intervalo de confiança do índice de eficiência por bootstrap homogêneo

tmp=boot.sw98(XOBS=x,YOBS=y,DHAT=dhat, RTS=1, ORIENTATION=1,NREP=2000)

## Formando tabela com os resultados

n=ncol(x)      #número de DMUs
table.in=matrix(nrow=n,ncol=7)
table.in[,1]=c(1:n)
table.in[,2]=dhat
table.in[,3]=dhat-tmp$bias   #estiamtiva de correção de viés
table.in[,4]=tmp$bias
table.in[,5]=tmp$var
table.in[,6:7]=tmp$conf.int
table.in[1:9,1]=paste(" ",table.in[1:9,1],sep="")

## Formatando a tabela

table.in[,2]=ifelse(nchar(table.in[,2])==1,
                    paste(table.in[,2],".",sep=""), table.in[,2])
table.in[,2:7]=paste(table.in[,2:7],"000000",sep="")
table.in[,c(2:3,5:7)]=substr(table.in[,c(2:3,5:7)],1,6)
table.in[,4]=substr(table.in[,4],1,7)
table.in=paste(table.in[,1]," & ",
               table.in[,2]," & ",
               table.in[,3]," & ",
               table.in[,4]," & ",
               table.in[,5]," & ",
               table.in[,6]," & ",
               table.in[,7]," \\",sep="")

## Coluna 1 = DMU?s (1-85)
## Coluna 2 = índice de eficiência Shepard orientação produto
## Coluna 3 = índice de eficiência corrigido, extrai-se o viés
## Coluna 4 = viés do índice
## Coluna 5 = variância
## Coluna 6-7 = intervalo de confiança.

## Exibir tabela

table.in

## Gerar estimadores dea com correção de viés

dea_corr=(1/(dhat-tmp$bias))  # invertido para uso no segundo estágio como variável explicada

DEA=matrix(nrow=85,ncol=1)   # formar matriz para inserir dentro do banco de dados
DEA[,1]=dea_corr
agro[,24]=DEA              # inverso do estimador dea corrigido inserido dentro do banco de dados como V24

################################################################################################
## Segundo Estágio (DEA regredido em relação às variáveis ambientais corrigido por Bootstrap) ##
################################################################################################


############################  Formulação Modelo  #########################

## Definição das variáveis

DEA<-agro[,24]  ## Atentar para definir as variáveis (agro[,24]=V24=dea_corr)

# V24 = dea_corr : inverso dos parâmetros de eficiência com viés corrigido por bootstrap

Z1<-c("z1","z2","z10")                             ## modelo 1
Z2<-c("z3","z4","z8")		                  ## modelo 2
Z3<-c("z5","z6","z7","z9")           	         ## modelo 3


## As variáveis já estão, na base de dados disponibilizadas,

alfa<-c(0.005,0.025,0.050,0.5,0.950,0.975,0.995) ## Níveis de significância

L=2000		## Número de repetições do bootstrap

dea_o<-DEA

z1<-agro[,Z1]
z2<-agro[,Z2]
z3<-agro[,Z3]



############################        Modelo 1       #########################



##* ef<-which(dea_o==1)      ##* comandos utilizados para retirar dmu's eficientes

##* Y_dea<-dea_o[-ef]

##* zi<-as.matrix(z[-ef,])

ni<-nrow(z1)

cons=matrix(rep(1,ni,nrow=ni,ncol=1))

## Estimar o modelo TOBIT utilizando o pacote AER

fit1 = tobit(DEA ~ z1 + z2 + z10,left=1, data = agro)                        ## modelo 1 


## Extrair Parâmetros modelo 1

coef1<-coef(fit1)		                        # Coeficientes

desvpad1<-sqrt(diag(vcov(fit1)[-(length(coef1)+1),]))	# Desvio Padrão dos Coefientes

er1<-residuals(fit1)	                                # Resíduos

sdpb1<-sd(er1) 		                                # Desvio padrão do resíduo


#####  Passo iii  #####

## gerar planilha de coeficientes gerados

coefsboot<-matrix(,nrow=ncol(z1)+1,ncol=L)

## Calcular os valores da truncagem

Bzb1=as.matrix(cbind(1,z1))%*%coef1

medb=cons - Bzb1

rept=1

repeat{
  
  rept=1+rept
  
  #### Passo iii.1 Gerar erro aleatório truncado
  
  erb=matrix(,ncol=1,nrow=ni)
  
  for (n in 1:ni){
    
    erb[n,]=rnorm.trunc(1, t.left = medb[n], t.right = Inf, mean = 0, sigma = sdpb1)
  }
  
  
  #### Passo iii.2  gerar DEA - Corrigido a partir do erro gerado no passo anterior
  
  dea_ob2=Bzb1+erb
  
  #### Passo iii.3 Fazer novamente uma regressão truncada com dea gerado no passo anterior
  
  fit2b = tobit(dea_ob2 ~ z1 + z2 + z10,left = 1, data = agro)
  
  ## Extrair Coeficientes
  
  coefb2<-coef(fit2b)
  
  ## coletar resultado
  
  coefsboot[,(rept-1)]<-coefb2
  
  if(rept==(L+1)) break
  
}

#### Resultados modelo 1 ####

##  calcular coeficientes e intervalos de confiança

coef_mean1<-rowMeans(coefsboot, na.rm = TRUE)	        # Coeficientes

sd_mean1<-sd(t(coefsboot))				# Desvio Padrão

ci1 = t(apply(coefsboot, 1, quantile, probs=alfa,
              type = 8, na.rm = TRUE))		# Quantis

## Salvar resultado

Coeficientes1<-data.frame(coef1,desvpad1,
                          coef_mean1,sd_mean1,ci1)



############ Salvar resultados modelo 1 ############

#### Gerar planilha em csv com os resultados do modelo 1

write.table(Coeficientes1, file = "MODELO1.CSV", sep = ";",
            row.names = TRUE,col.names = TRUE)



############################        Modelo 2       #########################

##* ef<-which(dea_o==1)      ##* comandos utilizados para retirar dmu's eficientes

##* Y_dea<-dea_o[-ef]

##* zi<-as.matrix(z[-ef,])

ni<-nrow(z2)

cons=matrix(rep(1,ni,nrow=ni,ncol=1))

## Estimar o modelo TOBIT utilizando o pacote AER

fit2 = tobit(V24 ~ z3 + z4 + z8,left=1, data = agro) ## modelo 2


## Extrair Parâmetros modelo 2

coef2<-coef(fit2)		                        # Coeficientes

desvpad2<-sqrt(diag(vcov(fit2)[-(length(coef2)+1),]))	# Desvio Padrão dos Coefientes

er2<-residuals(fit2)	                                # Resíduos

sdpb2<-sd(er2) 		                                # Desvio padrão do resíduo


#####  Passo iii  #####

## gerar planilha de coeficientes gerados

coefsboot<-matrix(,nrow=ncol(z2)+1,ncol=L)

## Calcular os valores da truncagem

Bzb2=as.matrix(cbind(1,z2))%*%coef2

medb=cons - Bzb2

rept=1

repeat{
  
  rept=1+rept
  
  #### Passo iii.1 Gerar erro aleatório truncado
  
  erb=matrix(,ncol=1,nrow=ni)
  
  for (n in 1:ni){
    
    erb[n,]=rnorm.trunc(1, t.left = medb[n], t.right = Inf, mean = 0, sigma = sdpb2)
  }
  
  
  #### Passo iii.2  gerar DEA - Corrigido a partir do erro gerado no passo anterior
  
  dea_ob2=Bzb2+erb
  
  #### Passo iii.3 Fazer novamente uma regressão truncada com dea gerado no passo anterior
  
  fit2b = tobit(dea_ob2 ~ z3 + z4 + z8,left = 1, data = agro)
  
  ## Extrair Coeficientes
  
  coefb2<-coef(fit2b)
  
  ## coletar resultado
  
  coefsboot[,(rept-1)]<-coefb2
  
  if(rept==(L+1)) break
  
}

#### Resultados modelo 2 ####

##  calcular coeficientes e intervalos de confiança

coef_mean2<-rowMeans(coefsboot, na.rm = TRUE)	        # Coeficientes

sd_mean2<-sd(t(coefsboot))				# Desvio Padrão

ci2 = t(apply(coefsboot, 1, quantile, probs=alfa,
              type = 8, na.rm = TRUE))		# Quantis

## Salvar resultado

Coeficientes2<-data.frame(coef2,desvpad2,
                          coef_mean2,sd_mean2,ci2)



############ Salvar resultados modelo 2 ############

#### Gerar planilha em csv com os resultados do modelo 1

write.table(Coeficientes2, file = "MODELO2.CSV", sep = ";",
            row.names = TRUE,col.names = TRUE)



############################        Modelo 3       #########################



##* ef<-which(dea_o==1)      ##* comandos utilizados para retirar dmu's eficientes

##* Y_dea<-dea_o[-ef]

##* zi<-as.matrix(z[-ef,])

ni<-nrow(z3)

cons=matrix(rep(1,ni,nrow=ni,ncol=1))

## Estimar o modelo TOBIT utilizando o pacote AER

fit3 = tobit(V24 ~ z5 + z6 + z7 + z9,left=1, data = agro)                    ## modelo 3


## Extrair Parâmetros modelo 1

coef3<-coef(fit3)		                        # Coeficientes

desvpad3<-sqrt(diag(vcov(fit3)[-(length(coef3)+1),]))	# Desvio Padrão dos Coefientes

er3<-residuals(fit3)	                                # Resíduos

sdpb3<-sd(er3) 		                                # Desvio padrão do resíduo


#####  Passo iii  #####

## gerar planilha de coeficientes gerados

coefsboot<-matrix(,nrow=ncol(z3)+1,ncol=L)

## Calcular os valores da truncagem

Bzb3=as.matrix(cbind(1,z3))%*%coef3

medb=cons - Bzb3

rept=1

repeat{
  
  rept=1+rept
  
  #### Passo iii.1 Gerar erro aleatório truncado
  
  erb=matrix(,ncol=1,nrow=ni)
  
  for (n in 1:ni){
    
    erb[n,]=rnorm.trunc(1, t.left = medb[n], t.right = Inf, mean = 0, sigma = sdpb3)
  }
  
  
  #### Passo iii.2  gerar DEA - Corrigido a partir do erro gerado no passo anterior
  
  dea_ob2=Bzb3+erb
  
  #### Passo iii.3 Fazer novamente uma regressão truncada com dea gerado no passo anterior
  
  fit2b = tobit(dea_ob2 ~ z5 + z6 + z7 + z9,left = 1, data = agro)
  
  ## Extrair Coeficientes
  
  coefb2<-coef(fit2b)
  
  ## coletar resultado
  
  coefsboot[,(rept-1)]<-coefb2
  
  if(rept==(L+1)) break
  
}

#### Resultados modelo 3 ####

##  calcular coeficientes e intervalos de confiança

coef_mean3<-rowMeans(coefsboot, na.rm = TRUE)	        # Coeficientes

sd_mean3<-sd(t(coefsboot))				# Desvio Padrão

ci3 = t(apply(coefsboot, 1, quantile, probs=alfa,
              type = 8, na.rm = TRUE))		# Quantis

## Salvar resultado

Coeficientes3<-data.frame(coef3,desvpad3,
                          coef_mean3,sd_mean3,ci3)



############ Salvar resultados modelo 3 ############

#### Gerar planilha em csv com os resultados do modelo 1

write.table(Coeficientes3, file = "MODELO3.CSV", sep = ";",
            row.names = TRUE,col.names = TRUE)




