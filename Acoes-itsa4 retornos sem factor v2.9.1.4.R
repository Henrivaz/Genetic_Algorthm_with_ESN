#A partir daqui AG para escolher a distribuição
#Treinando e validando com retornos sem factor e testando com retornos sem factor

#USANDO O MAE

#getwd()
setwd("C:\\Users\\henri\\OneDrive\\Área de Trabalho\\BolsaFURG\\ArqItsa4")

#BIBLIOTECAS
library(ggplot2)               #Plota os gráficos
library(PerformanceAnalytics)  #Assimetria e curtose
library(GA)                    #Algoritmo genético
library(pracma)                #Tic e Toc (demarcadores de tempo de execução)
library(rugarch)
library(fGarch)
library(skewt)                 #T assimetrica


#SÉRIES TEMPORAIS
#data = as.matrix(read.csv2('Retornos ITSA4_close sem factor_2001-30-06-2023.txt',header=F))       #5571 sem factor
#data = as.numeric(data)

#data
#length(data)
#View(data)                    # 5571 observações
#describe(data)
#summary(data)

#data_simulate = as.matrix(read.csv2('Retornos Log simulados Arma(7,0) 02-01-2001-28-12-2018 v1.3.txt',header=F))
#data_simulate = as.numeric(data_simulate)


#data_simulate = as.matrix(read.csv2('Retornos Log simulados v1.4 Arma(7,0) 02-01-2001-28-12-2018.txt',header=F))
#data_simulate = as.numeric(data_simulate)

#data_simulate = as.matrix(read.csv2('Retornos Log simulados v1.5 Arma(7,0) 02-01-2001-28-12-2018.txt',header=F))
#data_simulate = as.numeric(data_simulate)

data_simulate = as.matrix(read.csv2('Retornos Log simulados v1.6 Arma(7,0) 02-01-2001-28-12-2018.txt',header=F))
data_simulate = as.numeric(data_simulate)



#TAMANHOS TREINO, VALIDAÇÃO E TESTE
#treino   = 4453                #80% (02/01/2001 a 28/12/2018)  
treino   = 3335                #60% (02/01/2001 a 23/06/2014) 
valida   = 1118                #20% (24/06/2014 a 28/12/2018)   
#teste    = 1118                #20% (02/01/2019 a 30/6/2023) 

inSize = outSize = 1


#ANÁLISE DOS DADOS
#summary(data)                                 #min=, median=, mean=, max=
#skewness(data,method="moment")                #assimetria = (|assimétrica|>0.5)
#kurtosis(data,method="excess")                #curtose    = (platicúrtica<0)

par(mfrow=c(1,1))                               #plotar interaçoes
#hist(data)
#hist(data, ylim=c(0,1000), breaks = (c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)),xlab = "Closing price of ITSA4 stock (raw data)" ,main = "Histogram of daily closing price without factor \n of ITSA4 stock between 2001 and 30-06-2023")
#boxplot(data,ylim=c(0,15),ylab = "Closing price of ITSA4 stock (R$)",main="(a)") #com outliers

#par(mfrow=c(1,3))
#data_plot_train = data[1:treino]
#qplot(x = 1:treino , y = data_plot_train, geom = 'line') + geom_line(color = 'darkblue') + 
  #labs(x = 'Training days' , y = 'Price (R$)' , title = "Daily closing price of ITSA4 stock between 02-01-2001 and 23-06-2014") + geom_hline(yintercept = mean(data_plot_train) , color = 'red')

#data_plot_valid = data[3336:4453]
#qplot(x = 3336:4453 , y = data_plot_valid, geom = 'line') + geom_line(color = 'darkblue') + 
  #labs(x = 'Test days' , y = 'Price (R$)' , title = "Daily closing price of ITSA4 stock between 24-06-2014 and 28-12-2018") + geom_hline(yintercept = mean(data_plot_valid) , color = 'red')

#data_plot_test = data[4454:5571]
#qplot(x = 4454:5571 , y = data_plot_test, geom = 'line') + geom_line(color = 'darkblue') + 
  #labs(x = 'Test days' , y = 'Price (R$)' , title = "Daily closing price of ITSA4 stock between 02-01-2019 and 30-06-2023") + geom_hline(yintercept = mean(data_plot_test) , color = 'red')


#ALGORITMO GENÉTICO
#Função de monitoramento
monitora <- function(obj){
  plot(obj)
  #ggplot(y = c(obj@summary[,1],obj@summary[,2],obj@summary[,6]), geom = 'line') 
  #plot(type="l",cbind(obj@summary[,1],obj@summary[,2],obj@summary[,6]))
  #conta = conta + 1L
  if (old_obj < obj@summary[[obj@iter]]){
     #linha_fitness = cbind(obj@iter,obj@summary[[obj@iter]],old_obj)
     #write.table(linha_fitness,file = "Dados ITSA4 melhores_fitness ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
     
     a_monit   = sum(obj@bestSol[[obj@iter]][1:17] * 2^(rev(seq(along=obj@bestSol[[obj@iter]][1:17])) - 1))/131071
     #a_monit   = binary2decimal(obj@bestSol[[obj@iter]][1:17])/131071
     if (a_monit == 0) {a_monit=7.62939453125e-6}

     sr_monit  = sum(obj@bestSol[[obj@iter]][18:34] * 2^(rev(seq(along=obj@bestSol[[obj@iter]][18:34])) - 1))/131071
     #sr_monit  = binary2decimal(obj@bestSol[[obj@iter]][18:34])/131071
     if (sr_monit == 0) {sr_monit=7.62939453125e-6}
     
     iL_monit  = sum(obj@bestSol[[obj@iter]][35:41] * 2^(rev(seq(along=obj@bestSol[[obj@iter]][35:41])) - 1))
     #iL_monit  = binary2decimal(obj@bestSol[[obj@iter]][35:41])
     iL_monit  = iL_monit+2
     
     tr_monit  = sum(obj@bestSol[[obj@iter]][42:46] * 2^(rev(seq(along=obj@bestSol[[obj@iter]][42:46])) - 1))
     #tr_monit  = binary2decimal(alg_gen@bestSol[[i]][42:46])
     tr_monit  = tr_monit+2
     
     reg_monit = sum(obj@bestSol[[obj@iter]][47:55] * 2^(rev(seq(along=obj@bestSol[[obj@iter]][47:55])) - 1))/511
     #reg_monit = binary2decimal(alg_gen@bestSol[[obj@iter]][47:55])/511
     reg_monit = (reg_monit + 1e-4)*(1e-4 - 1e-6)
     
     linha_bestSol = cbind(obj@iter,a_monit,sr_monit,iL_monit,tr_monit,reg_monit,obj@summary[[obj@iter]])
     write.table(linha_bestSol,file = "Dados ITSA4 bestSol melhores_fitness mae_otim40x60 retornos sem factor 10000_3.csv", append = T, row.names = F, col.names = F, sep = "\t", eol = "\n")     
     
     old_obj  <<- obj@summary[[obj@iter]]
    }
  #sumario[[obj@iter]] <<- obj@summary[[obj@iter]]   #Primeiro valor do sumário (melhor)
  #sumario[[obj@iter]] <<- obj@summary[[obj@iter]][]
  #sumario[[obj@iter]] <<- c(obj@summary[obj@iter,]) #Todos valores do sumário
  sumario[[obj@iter]] <<- c(obj@summary[[obj@iter,1]],obj@summary[[obj@iter,2]],obj@summary[[obj@iter,6]])   #Todos específicos do sumário
}

#Dados para treino e validação
treino_valida = data_simulate[1:(treino+valida)]                                       #Sem factor
#head(treino_valida[1:treino])
#tail(treino_valida[1:treino])
#treino_valida[3118]
#treino_valida
length(treino_valida)         #4453

#Para cálculo do RRSE
train             = data_simulate[1:treino]                                               #Sem factor
par(mfrow=c(1,1))
#plot(train,   type="l")
train_media       = mean(train)

validate          = data_simulate[(treino+1):(treino+valida)]                             #Sem factor
#plot(validate,type="l")
valid_media       = mean(validate)

#Função de Fitness do AG
#f <- function(i, j, k){
f <- function(cromossoma=c()){

   a_GA             = binary2decimal(cromossoma[1:17])                 #real
   a_GA             = a_GA/131071                                      #taxa de vazão
   #a_GA             = sum(cromossoma[1:17] * 2^(rev(seq(along=cromossoma[1:17])) - 1))/131071
   if (a_GA == 0)    {a_GA=7.62939453125e-6}
   
   sr_GA            = binary2decimal(cromossoma[18:34])                #real
   sr_GA            = sr_GA/131071                                     #raio espectral
   #sr_GA            = sum(cromossoma[18:34] * 2^(rev(seq(along=cromossoma[18:34])) - 1))/131071
   if (sr_GA == 0)   {sr_GA=7.62939453125e-6}
   
   initLen_GA        = binary2decimal(cromossoma[35:41])               #int
   #initLen_GA       = sum(cromossoma[35:41] * 2^(rev(seq(along=cromossoma[35:41])) - 1))
   initLen_GA       = initLen_GA+2                                     #valores iniciais desconsiderados
   
   tam_reservoir_GA = binary2decimal(cromossoma[42:46])                   #int
   #tam_reservoir_GA = sum(cromossoma[42:46] * 2^(rev(seq(along=cromossoma[42:46])) - 1))
   tam_reservoir_GA = tam_reservoir_GA+2                               #tamanho do reservatório
   
   #reg_GA           = binary2decimal(cromossoma[24:24])/100000000      #real
   reg_GA           = binary2decimal(cromossoma[47:55])/511            #coeficiente de regularização da regressão
   #reg_GA           = sum(cromossoma[47:55] * 2^(rev(seq(along=cromossoma[47:55])) - 1))/511
   reg_GA           = (reg_GA + 1e-4)*(1e-4 - 1e-6)
   #if (reg_GA == 0) {reg_GA=1e-9}      
   
   distr_Win = "Uniforme, clássica"
   #Win_minus = -sum(cromossoma[56:67] * 2^(rev(seq(along=cromossoma[56:67])) - 1))/4095
   #Win_plus  =  sum(cromossoma[68:79] * 2^(rev(seq(along=cromossoma[68:79])) - 1))/4095
   #Win_GA    = matrix(runif(tam_reservoir_GA*(1+inSize),Win_minus,Win_plus),tam_reservoir_GA)
   #Win_par   = sum(cromossoma[56:67] * 2^(rev(seq(along=cromossoma[56:67])) - 1))/4095
   #Win_GA    = matrix(runif(tam_reservoir_GA*(1+inSize),-Win_par,Win_par),tam_reservoir_GA)
   Win_GA    = matrix(runif(tam_reservoir_GA*(1+inSize),-1.0,1.0),tam_reservoir_GA)
   
   distr_W   = "Uniforme, clássica"
   #W_minus   = -sum(cromossoma[80:91] * 2^(rev(seq(along=cromossoma[80:91])) - 1))/4095
   #W_plus    =  sum(cromossoma[92:103] * 2^(rev(seq(along=cromossoma[92:103])) - 1))/4095 
   #W_GA      = matrix(runif(tam_reservoir_GA*tam_reservoir_GA,W_minus,W_plus),tam_reservoir_GA)
   #W_par     = sum(cromossoma[68:79] * 2^(rev(seq(along=cromossoma[68:79])) - 1))/4095
   #W_GA      = matrix(runif(tam_reservoir_GA*tam_reservoir_GA,-W_par,W_par),tam_reservoir_GA)
   W_GA      = matrix(runif(tam_reservoir_GA*tam_reservoir_GA,-1.0,1.0),tam_reservoir_GA)
   
   #distr_Win = "Normal_clássica"
   #Win_GA    = matrix(rnorm(tam_reservoir_GA*(1+inSize),mean=-0.08,sd=1),tam_reservoir_GA)
   
   #distr_W    = "Normal_clássica"
   #W_GA       = matrix(rnorm(tam_reservoir_GA*tam_reservoir_GA,mean=-0.08,sd=1),tam_reservoir_GA)
   
   #fitdist(treino_valida,distribution="std")
   #Win_GA    = matrix(rt(tam_reservoir_GA*(1+inSize),df=4.6135422773),tam_reservoir_GA) #Cria matriz
   #distr_Win = "t_de_Student_rugarch"
   
   #fitdist(treino_valida,distribution="std") 
   #distr_W  = "t_de_Student_rugarch"
   #W_GA      = matrix(rt(tam_reservoir_GA*tam_reservoir_GA,df=4.6135422773),tam_reservoir_GA)
   
   #fitdist(distribution = "ged", treino_valida) 
   #Win_GA = matrix(rged(tam_reservoir_GA*(1+inSize), mean= 5.573398e-06, sd= 2.202071e-02, nu= 1.164987e+00),tam_reservoir_GA)
   #distr_Win= "Ged_distr"
   
   #fitdist(distribution = "ged", treino_valida)
   #W_GA = matrix(rged(tam_reservoir_GA*tam_reservoir_GA,mean= 5.573398e-06, sd= 2.202071e-02, nu= 1.164987e+00), tam_reservoir_GA)
   #distr_W= "Ged_distr"
   
   #fitdist(treino_valida,distribution="sstd")
   #distr_Win = "Assimétrica_t_as_rugarch"
   #Win_GA    = matrix(rskt(tam_reservoir_GA*(1+inSize), df=4.6135422773,gamma= 1.0475870739) ,tam_reservoir_GA) #Cria matriz
   
   #fitdist(treino_valida,distribution="sstd")
   #distr_W   = "Assimétrica_t_as_rugarch"
   #W_GA      = matrix(rskt(tam_reservoir_GA*tam_reservoir_GA, df=4.6135422773,gamma= 1.0475870739),tam_reservoir_GA) #Cria matriz
   
   
     
     rhoW_GA = abs(eigen(W_GA,only.values=TRUE)$values[1])                  #Autovalores da matriz W
   W_GA = sr_GA * W_GA / rhoW_GA
   X = matrix(0,1+inSize+tam_reservoir_GA,treino-initLen_GA)
   Yt = matrix(treino_valida[(initLen_GA+2):(treino+1)],1)
   x = rep(0,tam_reservoir_GA)                                            #Preenche o vetor x (do tamanho do reservatório) com 0

#Treinando   
   for (t in 1:treino){
       u = treino_valida[t]
       x = (1-a_GA)*x + a_GA*tanh( Win_GA %*% rbind(1,u) + W_GA %*% x )         #tanh é a função de ativação
       if (t > initLen_GA)
          X[,t-initLen_GA] = rbind(1,u,x)
      }
   X_T = t(X)
   Wout_GA = Yt %*% X_T %*% solve(X %*% X_T + reg_GA*diag(1+inSize+tam_reservoir_GA))
   
#Prevendo os dados de teste
   Y = matrix(0,outSize,(valida-1))
   u = treino_valida[treino+1]

   for (t in 1:(valida-1)){
       x = (1-a_GA)*x + a_GA*tanh( Win_GA %*% rbind(1,u) + W_GA %*% x )
       y = Wout_GA %*% rbind(1,u,x)
       Y[,t] = y
       u = treino_valida[treino+t+1]
   }

#Prevendo os dados de treino
   Ytr = matrix(0,outSize,(treino-1))
   u   = treino_valida[1]
   
   for (j in 1:(treino-1)) {
       x = (1-a_GA)*x + a_GA*tanh( Win_GA %*% rbind(1,u) + W_GA %*% x )
       y = Wout_GA %*% rbind(1,u,x)
       Ytr[,j] = y
       u = treino_valida[j+1]
      }
   
#Cálculo dos erros para otimização
   #mse_treino_GA  = (sum((        treino_valida[1:treino] - Ytr[outSize,1:treino])^2 )/treino)
   #mse_valida_GA  = (sum((        treino_valida[(treino+1):(treino+valida)]-Y[outSize,1:valida])^2)/valida)
   rmse_treino_GA = sqrt(sum((    treino_valida[2:treino] - Ytr[outSize,1:(treino-1)])^2)/(treino-1))
   rmse_valida_GA = sqrt(sum((    treino_valida[(treino+2):(treino+valida)]-Y[outSize,1:(valida-1)])^2)/(valida-1))
   rrse_treino_GA = sqrt(sum((    treino_valida[2:treino] - Ytr[outSize,1:(treino-1)])^2)/sum((treino_valida[2:treino] - train_media)^2))
   rrse_valida_GA = sqrt(sum((    treino_valida[(treino+2):(treino+valida)]-Y[outSize,1:(valida-1)])^2)/sum((treino_valida[(treino+2):(treino+valida)]-valid_media)^2))
   mae_treino_GA  =      sum(abs( treino_valida[2:treino] - Ytr[outSize,1:(treino-1)])) /(treino-1)
   mae_valida_GA  =      sum(abs( treino_valida[(treino+2):(treino+valida)]-Y[outSize,1:(valida-1)]))/(valida-1)
   #mdae_treino_GA =      sum(abs((treino_valida[2:treino] - Ytr[outSize,1:(treino-1)])  /treino_valida[2:treino]))/(treino-1)
   #mdae_valida_GA =      sum(abs((treino_valida[(treino+2):(treino+valida)]-Y[outSize,1:(valida-1)])/treino_valida[(treino+2):(treino+valida)]))/(valida-1)
   mdae_treino_GA =      median(abs( treino_valida[2:treino] - Ytr[outSize,1:(treino-1)]))
   mdae_valida_GA =      median(abs( treino_valida[(treino+2):(treino+valida)]-Y[outSize,1:(valida-1)]))
   
#Fitness = erro de treinamento * 0.4 e erro de validação * 0.6
   #otimiza            = (-rmse_treino_GA*0.4 - rmse_valida_GA*0.6)
   otimiza            = (-mae_treino_GA*0.4 - mae_valida_GA*0.6)
   contar           <<- contar + 1L
   old_otim_backup  <<- old_otim
   old_rmse_backup  <<- old_rmse
   
#Gravação dos resultados parciais   
   Win_GA_t         <<- t(Win_GA)
   W_GA_t           <<- t(W_GA)
   Wout_GA_t        <<- t(Wout_GA)
   if (old_otim < otimiza){
      #linha_entrada      = cbind(contar,Win_GA_t, distr_Win_GA)
      #linha_entrada      = cbind(contar,Win_GA_t,"Uniforme",Win_minus,Win_plus)
      #linha_entrada      = cbind(contar,Win_GA_t,"Uniforme",-Win_par,Win_par)
      linha_entrada      = cbind(contar,Win_GA_t,"Uniforme",-1.0,1.0)
      #linha_entrada      = cbind(contar,Win_GA_t,"Normal_classica",distr_Win,0,1)
      #linha_entrada      = cbind(contar,Win_GA_t,distr_Win,4.6135422773)  #T de student 
      #linha_entrada      = cbind(contar,Win_GA_t,distr_Win,5.573398e-06, 2.202071e-02, 1.164987e+00) #Ged
      #linha_entrada      = cbind(contar,Win_GA_t,distr_Win,4.6135422773,1.0475870739) #T de Student Assimetrica
      
      write.table(linha_entrada,     file = "Dados ITSA4 Win ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
      #linha_reservatorio = cbind(contar,W_GA_t, distr_W_GA)
      #linha_reservatorio = cbind(contar,W_GA_t, "Uniforme",W_minus,W_plus)
      #linha_reservatorio = cbind(contar,W_GA_t, "Uniforme",-W_par,W_par)
      linha_reservatorio = cbind(contar,W_GA_t, "Uniforme",-1.0,1.0)
      #linha_reservatorio = cbind(contar,W_GA_t,"Normal_classica",  distr_W,0,1)
      #linha_reservatorio = cbind(contar,W_GA_t,distr_W, 5.573398e-06, 2.202071e-02, 1.164987e+00)  #Ged
      #linha_reservatorio = cbind(contar,W_GA_t,"T de Student",4.6135422773)
      #linha_reservatorio = cbind(contar,W_GA_t,"T de Student Assimetrica",4.6135422773,1.0475870739)
      
      write.table(linha_reservatorio,file = "Dados ITSA4 W reservatório ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
      linha_saida        = cbind(contar,Wout_GA_t)
      write.table(linha_saida,       file = "Dados ITSA4 Wout ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
      old_otim         <<- otimiza
     }
   linha_ESN   = cbind(contar,a_GA,sr_GA,initLen_GA,tam_reservoir_GA,reg_GA,mae_treino_GA,mae_valida_GA,otimiza,old_otim_backup,rmse_treino_GA,rmse_valida_GA,rmse_treino_GA*0.4+rmse_valida_GA*0.6,mdae_treino_GA ,mdae_valida_GA,mdae_treino_GA*0.4+mdae_valida_GA*0.6,rrse_treino_GA,rrse_valida_GA,rrse_treino_GA*0.4+rrse_valida_GA*0.6)
   write.table(linha_ESN,         file = "Dados ITSA4 ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
  
   return(otimiza)
}

#Inicialização dos arquivos para gravação
linha_ESN          = cbind("Parâmetro","Contar","a","sr","initLen","tam_reservoir","reg","MAE_treino40%","MAE_valida60%","Otimiza","old_Otimiza","RMSE_treino","RMSE_valida","RMSE_t0.4+RMSE_v0.6","mdae_treino","mdae_valida","mdae_t0.4+mdae_v0.6","RRSE_treino","RRSE_valida","RRSE_t0.4+RRSE_v0.6")
write.table(linha_ESN,file     = "Dados ITSA4 ESN_mae_otim40x60 retornos sem factor 10000_3.csv", append = F, row.names = F, col.names = F, sep = "\t", eol = "\n")
linha_bestSol      = cbind("Época","a","sr","iL","tr","reg","Fitness")
write.table(linha_bestSol,file = "Dados ITSA4 bestSol melhores_fitness mae_otim40x60 retornos sem factor 10000_3.csv", append = F, row.names = F, col.names = F, sep = "\t", eol = "\n")     
linha_entrada      = cbind("Época","Win","Distribuição","Win_Minus","Win_Plus")
write.table(linha_entrada,file = "Dados ITSA4 Win ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = F, col.names = F, sep = "\t", eol = "\n")
linha_reservatorio = cbind("Época","W","Distribuição","W_Minus","W_Plus")
write.table(linha_reservatorio,file = "Dados ITSA4 W reservatório ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = F, col.names = F, sep = "\t", eol = "\n")
linha_saida        = cbind("Época","Wout")
write.table(linha_saida,file = "Dados ITSA4 Wout ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = F, col.names = F, sep = "\t", eol = "\n")
#linha_fitness = cbind("Época","Fitness","old_fitness")
#write.table(linha_fitness,file  = "Dados ITSA4 melhor_fitness ESN_mae_otim40X60 brutosretornos sem factor 10000_3.csv", append = F, col.names = F, sep = "\t", eol = "\n")

#Inicialização dos dados
old_obj  <- -Inf
old_otim <- -Inf
old_rmse <-  Inf
old_mae   =  Inf
conta    <- contar <- 0L
rmse_teste= rmse_treino = mae_teste = mae_treino = mdae_teste = mdae_treino = otimiza = 0
itera     = 10000L
sumario  <- vector("list", length = itera)
#semente   = 17161

#Chamada da função do AG
tic()
alg_gen <- ga("binary", fitness=function(x) f(x), nBits = 67, popSize = 10, pcrossover = 0.8, population = gabin_Population,
              selection = gabin_tourSelection, crossover = gabin_spCrossover, mutation = gabin_raMutation, run=3500, #máximo de iterações sem melhora
              #selection = gabin_nlrSelection, crossover = gabin_uCrossover, mutation = gabin_raMutation, run=2500, #máximo de iterações sem melhora
              #selection = gabin_lrSelection, crossover = gabin_uCrossover, mutation = gabin_raMutation, run=2500, #máximo de iterações sem melhora
              #selection = gabin_rwSelection, crossover = gabin_spCrossover, mutation = gabin_raMutation, run=2500, #máximo de iterações sem melhora
              #pmutation=0.1, elitism = 1, parallel = 8, maxiter= itera,     #n de gerações - aumentar depois
              #pmutation=0.05, elitism = 0.999, parallel = F, maxiter= 1535, #n de gerações - aumentar depois
              pmutation=0.1, elitism = 1, parallel = F, maxiter= itera,      #n de gerações - aumentar depois
              #keepBest = T,  monitor=plot, optim = F, seed = semente)
              #keepBest = T,  monitor=monitora, optim = F, seed = semente)
              keepBest = T,  monitor=monitora, optim = F)                #o keepBest salva os melhores individuos de cada geração
              #keepBest = T, monitor=plot, seed=107, optim = T, optimArgs = list(method = "Brent", poptim = 0.1)) #o keepBest salva os melhores individuos de cada geração
              #keepBest = T, monitor=plot, seed=107, optim = T, optimArgs = list(method = "CG",    poptim = 0.1)) #o keepBest salva os melhores individuos de cada geração
              #keepBest = T, monitor=plot, seed=107, optim = T, optimArgs = list(method = "Nelder-Mead", poptim = 0.1)) #o keepBest salva os melhores individuos de cada geração
              #keepBest = T, monitor=plot, seed=107, optim = T, optimArgs = list(method = "BFGS",  poptim = 0.1)) #o keepBest salva os melhores individuos de cada geração
              #keepBest = T, monitor=plot, seed=107, optim = T, optimArgs = list(method = "SANN",  poptim = 0.1)) #o keepBest salva os melhores individuos de cada geração
toc()
#sumario
#View(sumario)
#old_obj
#par(mfrow=c(1,1))
#View(alg_gen)
#alg_gen@summary
#alg_gen@fitness[1, ]
#plot(alg_gen@solution[1, ])

#Gravação dos resultados finais
linha_ESN          = cbind("Parâmetro","Contar","a","sr","initLen","tam_reservoir","reg","MAE_treino40%","MAE_valida60%","Otimiza","old_Otimiza","RMSE_treino","RMSE_valida","RMSE_t0.4+RMSE_v0.6","mdae_treino","mdae_valida","mdae_t0.4+mdae_v0.6","RRSE_treino","RRSE_valida","RRSE_t0.4+RRSE_v0.6")
write.table(linha_ESN,file     = "Dados ITSA4 ESN_mae_otim40x60 retornos sem factor 10000_3.csv", append = T, row.names = F, col.names = F, sep = "\t", eol = "\n")

linha_bestSol      = cbind("Época","a","sr","iL","tr","reg","Fitness")
write.table(linha_bestSol,file = "Dados ITSA4 bestSol melhores_fitness mae_otim40x60 retornos sem factor 10000_3.csv", append = T, row.names = F, col.names = F, sep = "\t", eol = "\n")     

linha_entrada      = cbind("Época","Win","Distribuição")
write.table(linha_entrada,file = "Dados ITSA4 Win ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
linha_reservatorio = cbind("Época","W","Distribuição")
write.table(linha_reservatorio,file = "Dados ITSA4 W reservatório ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
linha_saida        = cbind("Época","Wout")
write.table(linha_saida,file = "Dados ITSA4 Wout ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")

linha_resumo       = cbind("Melhor_Fitness","Fitness_médio","Fitness_mínimo")
write.table(linha_resumo,file  = "Dados ITSA4 resumo fitness ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = F, col.names = F, sep = "\t", eol = "\n")
linha_resumo       = cbind(sumario)
write.table(linha_resumo,file  = "Dados ITSA4 resumo fitness ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")
linha_resumo       = cbind("Melhor_Fitness","Fitness_médio","Fitness_mínimo")
write.table(linha_resumo,file  = "Dados ITSA4 resumo fitness ESN_mae_otim40X60 retornos sem factor 10000_3.csv", append = T, col.names = F, sep = "\t", eol = "\n")

itera   = 8381

linha_bestSol      = cbind("Época","a","sr","iL","tr","reg")
write.table(linha_bestSol,file = "Dados ITSA4 bestSol por época mae_otim40x60 retornos sem factor 10000_3.csv", append = F, row.names = F, col.names = F, sep = "\t", eol = "\n")
for (i in 1:itera) {
    resulta_a      = binary2decimal(alg_gen@bestSol[[i]][1:17])/131071
    if (resulta_a == 0) {resulta_a=7.62939453125e-6}
    resulta_sr     = binary2decimal(alg_gen@bestSol[[i]][18:34])/131071
    if (resulta_sr == 0) {resulta_sr=7.62939453125e-6}
    resulta_iL     = binary2decimal(alg_gen@bestSol[[i]][35:41])
    resulta_iL     = resulta_iL+2
    #if (resulta_iL == 0) {resulta_iL=2}
    resulta_tr     = binary2decimal(alg_gen@bestSol[[i]][42:46])
    resulta_tr     = resulta_tr+2
    resulta_reg    = binary2decimal(alg_gen@bestSol[[i]][47:55])/511
    resulta_reg    = (resulta_reg + 1e-4)*(1e-4 - 1e-6)

    linha_bestSol  = cbind(i,resulta_a,resulta_sr,resulta_iL,resulta_tr,resulta_reg)
    write.table(linha_bestSol,file = "Dados ITSA4 bestSol por época mae_otim40x60 retornos sem factor 10000_3.csv", append = T, row.names = F, col.names = F, sep = "\t", eol = "\n")
}
linha_bestSol = cbind("Época","a","sr","iL","tr","reg")
write.table(linha_bestSol,file = "Dados ITSA4 bestSol por época mae_otim40x60 retornos sem factor 10000_3.csv", append = T, row.names = F, col.names = F, sep = "\t", eol = "\n")

