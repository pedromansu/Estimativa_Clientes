## Funções para analise

##Medidas Resumo

medidasresumo=function(dados){
  resumo=resultados=matrix(0, nrow = 5, ncol = dim(dados)[2])
  for (i in 1:5){
    for (j in 1:dim(dados)[2]){
      resumo[i,j]=summary(dados[,i])[j]
    }
  }
  return(resumo)
}


##Envelope POISSON

envelopePoisson = function(fit.model,ligacao){
  par(mfrow=c(1,1))
  X = model.matrix(fit.model)
  n = nrow(X)
  p = ncol(X)
  w = fit.model$weights
  W = diag(w)
  H = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h = diag(H)
  td = resid(fit.model,type="deviance")/sqrt((1-h))
  e = matrix(0,n,100)
  #
  for(i in 1:100){
    nresp = rpois(n, fitted(fit.model))
    fit = glm(nresp ~ X, family=poisson(link=ligacao))
    w = fit$weights
    W = diag(w)
    H = solve(t(X)%*%W%*%X)
    H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h = diag(H)
    e[,i] = sort(resid(fit,type="deviance")/sqrt(1-h))}
  #
  e1 = numeric(n)
  e2 = numeric(n)
  #
  for(i in 1:n){
    eo = sort(e[i,])
    e1[i] = (eo[2]+eo[3])/2
    e2[i] = (eo[97]+eo[98])/2}
  #
  med = apply(e,1,mean)
  faixa = range(td,e1,e2)
  
  qqnorm(td,xlab="Percentil da N(0,1)",
         ylab="Resíduo Componente do Desvio", ylim=faixa, pch=16,main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  #
  par(new=T)
  #
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1,main="")
  par(new=T)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1,main="")
  par(new=T)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2,main="")
}


# Diagnostico POISSON

diagnosticoPoisson=function(fit.model){
  X = model.matrix(fit.model)
  n = nrow(X)
  p = ncol(X)
  w = fit.model$weights
  W = diag(w)
  H = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h = diag(H)
  ts = resid(fit.model,type="pearson")/sqrt(1-h)
  td = resid(fit.model,type="deviance")/sqrt(1-h)
  di = (h/(1-h))*(ts^2)
  par(mfrow=c(2,2))
  a = max(td)
  b = min(td)
  plot(td,xlab="Indice", ylab="Residuo Componente do Desvio",
       ylim=c(b-1,a+1), pch=16,cex.axis=1.1,cex.lab=1.1,cex=1.1,cex.axis=1.1,cex.lab=1.1)
  abline(2,0,lty=2)
  abline(-2,0,lty=2)
  abline(0,0,lty=2)
  
  
  fited = fitted(fit.model)
  plot(fited ,td,xlab="valor ajustado (média)", ylab="Residuo Componente do Desvio",ylim=c(b-1,a+1), pch=16,
       main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  abline(2,0,lty=2)
  abline(-2,0,lty=2)
  abline(0,0,lty=2)
  
  #
  hist(td,xlab="Resíduo Componente do Desvio",ylab="densidade",probability=TRUE,main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  #
  eta = predict(fit.model)
  z = eta + resid(fit.model, type="pearson")/sqrt(w)
  plot(predict(fit.model),z,xlab="Preditor Linear",ylab="Variavel z", pch=16,main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  lines(smooth.spline(predict(fit.model), z, df=2))
}


# Envelope Binomial Negativa

envelnbin = function(fit.model){
  X = model.matrix(fit.model)
  n = nrow(X)
  p = ncol(X)
  fi = fit.model$theta
  w = fi*fitted(fit.model)/(fi + fitted(fit.model))
  W = diag(w)
  H = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h = diag(H)
  td = resid(fit.model,type="deviance")/sqrt(1-h)
  fi = fit.model$theta
  e = matrix(0,n,100)
  
  
  for(i in 1:100){
    resp = rnegbin(n, fitted(fit.model),fi)
    fit = glm.nb(resp ~ X,link=ligacaonbin)
    w = fit$weights
    W = diag(w)
    H = solve(t(X)%*%W%*%X)
    H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h = diag(H)
    e[,i] = sort(resid(fit,type="deviance")/sqrt(1-h))}
  #
  e1 = numeric(n)
  e2 = numeric(n)
  #
  for(i in 1:n){
    eo = sort(e[i,])
    e1[i] = (eo[2]+eo[3])/2
    e2[i] = (eo[97]+eo[98])/2}
  #
  med = apply(e,1,mean)
  faixa = range(td,e1,e2)
  
  qqnorm(td, xlab="Percentil da N(0,1)",
         ylab="Resíduo Componente do Desvio", ylim=faixa, pch=16, main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  par(new=T)
  #
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
  par(new=T)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
  par(new=T)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2, main="")
}


##DiagnosticoBinomialNegativa

diagnbin=function(fit.model){
  X = model.matrix(fit.model)
  n = nrow(X)
  p = ncol(X)
  fi = fit.model$theta
  w = fi*fitted(fit.model)/(fi + fitted(fit.model))
  W = diag(w)
  H = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h = diag(H)
  ts = resid(fit.model,type="pearson")/sqrt(1-h)
  td = resid(fit.model,type="deviance")/sqrt(1-h)
  di = (h/(1-h))*(ts^2)
  par(mfrow=c(2,2))
  a = max(td)
  b = min(td)
  plot(td,xlab="Índice", ylab="Resíduo Componente do Desvio",
       ylim=c(b-1,a+1), pch=16,cex=1.1,cex.axis=1.1,cex.lab=1.1)
  abline(2,0,lty=2)
  abline(-2,0,lty=2)
  abline(0,0,lty=2)
  #
  plot(fitted(fit.model),td,xlab="Valor Ajustado", 
       ylab="Resíduo Componente do Desvio", ylim=c(b-1,a+1), pch=16,cex=1.1,cex.axis=1.1,cex.lab=1.1)
  #
  abline(2,0,lty=2)
  abline(-2,0,lty=2)
  abline(0,0,lty=2)
  #
  hist(td,xlab="Resíduo Componente do Desvio",ylab="densidade",probability=TRUE,main="",cex=1.1,cex.axis=1.1,cex.lab=1.1)
  #
  w = fit.model$weights
  eta = predict(fit.model)
  z = eta + resid(fit.model, type="pearson")/sqrt(w)
  plot(predict(fit.model),z,xlab="Preditor Linear", 
       ylab="Variável z", pch=16,cex=1.1,cex.axis=1.1,cex.lab=1.1)
  abline(0,1)
}


##Descritiva

library(ggplot2)
library(grid)
library(gridExtra)
library(xtable)
dados=read.table("C:/Users/Pedro/Downloads/loja.txt")
colnames(dados)=c("Nclientes","Ncasas","Renda","Icasas","DistCONC","DistLOJA")

# (Analise Descritiva)

# Gráficos de Dispersão
d1=ggplot(dados,aes(Ncasas,Nclientes))+geom_point()+theme_minimal()
d2=ggplot(dados,aes(Renda,Nclientes))+geom_point()+theme_minimal()
d3=ggplot(dados,aes(Icasas,Nclientes))+geom_point()+theme_minimal()
d4=ggplot(dados,aes(DistCONC,Nclientes))+geom_point()+theme_minimal()
d5=ggplot(dados,aes(DistLOJA,Nclientes))+geom_point()+theme_minimal()
grid.arrange(d1,d2,d3,d4,d5)

# Medidas Resumo
sumario=as.matrix(0)
sumario=medidasresumo(dados)
rownames(sumario)=c("Ncasas","Renda","Icasas","DistCONC","DistLOJA")
colnames(sumario)=c("Min","1ºQ","2ºQ","Média","3ºQ","Max")


# Ajuste

library(MASS)

# Padronização
Nclientes=dados$Nclientes
Ncasas=scale(dados$Ncasas)
Renda=scale(dados$Renda)
Icasas=scale(dados$Icasas)
DistCONC=scale(dados$DistCONC)
DistLOJA=scale(dados$DistLOJA)

banco=data.frame(Nclientes,Ncasas,Renda,Icasas,DistCONC,DistLOJA)

#Splitando o banco de dados
idx=sample((1:110),size=75,replace=FALSE)
train=banco[idx,]
test=banco[-idx,]


# Ajuste do Modelo

fit.poisson=glm(Nclientes~.,data=train,family = "poisson")
summary(fit.poisson)


# Diagnostico

# Diagnostico do Modelo
diagnosticoPoisson(fit.poisson)
envelopePoisson(fit.poisson,"log")
AIC(fit.poisson)
BIC(fit.poisson)

desvio=deviance(fit.poisson)
n=length(Nclientes)
p=ncol(model.matrix(fit.poisson))
pvdesv=1-pchisq(desvio,df=n-p)


# Seleção de Modelos Stepwise

fit.poisson0=glm(Nclientes~1,family="poisson",data=train)
modelforward=stepAIC(fit.poisson0,scope=list(upper=fit.poisson),direction=c("both"),trace = FALSE)
modelbackward=stepAIC(fit.poisson,scope=list(lower=fit.poisson0),direction=c("both"),trace=FALSE)


# Validação na base de teste

qn=qnorm(0.975)
pred=predict(modelforward,type=c("response"),se.fit = TRUE, newdata = test)
mupred=pred$fit
semupred=pred$se.fit
liIC=mupred-qn*semupred
lsIC=mupred+qn*semupred

plot(test$Nclientes,mupred,pch=20,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="observado",ylab="predito")
abline(0,1,lwd=2)


# Ajuste com a binomial negativa


fit.bn=glm.nb(Nclientes ~ DistLOJA + DistCONC + Renda + Ncasas,data=train,link=log)
summary(fit.bn)
diagnbin(fit.bn)
abline(0,1)
ligacaonbin ="log"
envelnbin(fit.bn)

