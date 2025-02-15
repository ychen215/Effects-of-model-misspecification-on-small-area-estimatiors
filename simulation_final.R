# Point estimators
rm(list = ls(all = TRUE))
source("OBP_ner1.3.R")
library(pps)
library(sae)
library(dplyr)
library(DescTools)
library(psych)
library(OBPSAE) 
require(lme4)
require(rhnerm)

NoSim = 1000
m = 400
OBP_Xi = matrix(NA, NoSim, m)
OBP_FH_kv = matrix(NA, NoSim, m)
# OBP_FH_pol = matrix(NA, NoSim, m)
# OBP_FH_upol = matrix(NA, NoSim, m)
OBP_unit = matrix(NA, NoSim, m)
EBLUP = matrix(NA, NoSim, m)
EBLUP_Xi = matrix(NA, NoSim, m)
Direct = matrix(NA, NoSim, m)
Yibar_mat = matrix(NA, NoSim, m)

ar<-seq(1,m)  
ni=rep(4,m)
Ni=rep(1000,m)
N=sum(Ni)
n=sum(ni)

sigma2e.i<- rgamma(m, 3, 0.5)
X<- matrix(c(rlnorm(N, log(4.5)-0.5, 0.5)), nrow=N, ncol=1) #high
u<- NULL
sigma2u<- 1
u<-rnorm(m, 0, sqrt(sigma2u)) 
u<-rep(u, Ni)
e<-NULL
for (j in 1:m){e<- c(e,rnorm(Ni[j], 0, sqrt(sigma2e.i[j])))}

gr<- rep(1:m, each=Ni[1])
ar<- unique(gr)

y = 10 + u + e
pop.matrix<-cbind(y,X,gr)
pop<-as.data.frame(pop.matrix)
names(pop)<-c("y","x","area")     
Yibar=tapply(pop[,1],pop[,3],mean)

for (i in 1:NoSim) {
  set.seed(i+1000)
  
  Yibar_mat[i, ] = as.numeric(Yibar)
  XMean<-tapply(pop$x,pop$area,mean)
  
  s = sort(stratsrs(pop$area, ni))
  x.s = as.matrix(pop[s,]$x)
  y.s = pop[s,]$y
  regioncode.s = pop[s,]$area
  XsMean = cbind(tapply(x.s, regioncode.s, mean))
  Ysbar = cbind(tapply(y.s, regioncode.s, mean))
  
  Direct[i, ] = Ysbar
  
  mod_unit = obp_ner1.3(
    y = y.s,
    x = cbind(x.s),
    area = regioncode.s,
    ni = ni,
    Ni = Ni,
    Xibar = cbind(XMean),
    xibar = cbind(XsMean))
  OBP_unit[i, ] = mod_unit$OBP

  mod_Xi = obp_ner1.3(
  y = y.s,
  x = cbind(x.s),
  area = regioncode.s,
  ni = ni,
  Ni = Ni,
  Xibar = cbind(XMean),
  xibar = cbind(XMean))
  OBP_Xi[i, ] = mod_Xi$OBP
  #
  # di = tapply(pop[s,]$y,pop[s,]$area,var)/ni
  # modobpXFH = obpFH(Ysbar~XMean-1,errorvar = c(di))
  # OBP_FH_upol[i, ] = modobpXFH$theta.OBP
  #
  # di_pool = var(y.s)/ni
  # modobpXFH_PL = obpFH(Ysbar~XMean-1, errorvar = di_pool)
  # OBP_FH_pol[i, ] = modobpXFH_PL$theta.OBP

  Di_pop = tapply(pop[,1],pop[,3],var)/ni
  modobpXFH_kv = obpFH(Ysbar~XMean-1, errorvar = Di_pop)
  OBP_FH_kv[i, ] = modobpXFH_kv$theta.OBP

  # # EBLUP
  Xmean <- data.frame(ar, XMean)
  PopN  <- data.frame(ar, Ni)
  tmp.dta <- data.frame(y.s, x.s, regioncode.s)
  names(tmp.dta) <- c("ys", "xs", "area.s")
  est_npeblup = eblupBHF(
    formula = ys ~ -1 + xs,
    dom = area.s,
    meanxpop = Xmean,
    popnsize = PopN,
    method = "REML",
    data = tmp.dta
  )
  EBLUP[i, ] = est_npeblup$eblup$eblup
  
  dta.eblup.uc = data.frame(y.s, Xi = rep(XMean, ni), regioncode.s)
  names(dta.eblup.uc) <- c("ys", "xs", "area.s")
  mod_eblup_uc = eblupBHF(    
                  formula = ys ~ -1 + xs,
                  dom = area.s,
                  meanxpop = Xmean,
                  popnsize = PopN,
                  method = "REML",
                  data = dta.eblup.uc)
  EBLUP_Xi[i, ] = mod_eblup_uc$eblup$eblup
  plot(i, ylim = c(0, NoSim), main = paste(i, "/", NoSim))
}




OBP_UC = colMeans((OBP_Xi - Yibar_mat)^2)
OBP_FH = colMeans((OBP_FH_kv - Yibar_mat)^2)
# OBP_FH_PL = colMeans((OBP_FH_pol - Yibar_mat)^2)
# OBP_FH_NPL = colMeans((OBP_FH_upol - Yibar_mat)^2)
EBLUP_mse = colMeans((EBLUP - Yibar_mat)^2)
EBLUP_UC = colMeans((EBLUP_Xi - Yibar_mat)^2)
OBP_unit_mse = colMeans((OBP_unit - Yibar_mat)^2)
DIRECT_mse = colMeans((Direct - Yibar_mat)^2)

mean(OBP_UC)
mean(OBP_FH)
# mean(OBP_FH_PL)
# mean(OBP_FH_NPL)
mean(EBLUP_mse)
mean(EBLUP_UC)
mean(OBP_unit_mse)
mean(DIRECT_mse)

save.image("EBLUP_UC_32.RData")

Method = rep(c("OBP_UC", "OBP_FH", "OBP_FH_PL", "OBP_FH_NPL", "OBP_UNIT", "EBLUP", "DIRECT"), rep(m,7))
Method = rep(c("DIRECT", "OBP_UC", "OBP_FH","OBP_FH_NPL", "OBP_FH_PL", "OBP_UNIT", "EBLUP" ), rep(m,7))
Method = factor(Method, levels = c("DIRECT", "OBP_UC", "OBP_FH", "OBP_FH_NPL", "OBP_FH_PL", "OBP_UNIT", "EBLUP"))
Results1 = c(DIRECT_mse, OBP_UC, OBP_FH, OBP_FH_NPL, OBP_FH_PL, OBP_unit_mse, EBLUP_mse)
data1 = data.frame(Overall_MSPE = Results1, Method = Method)
library(ggplot2)
ggplot(data = data1) + 
  geom_boxplot(aes(y = Overall_MSPE, x = Method)) +
  xlab(NULL)+
  ylab(NULL)+
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x= element_text(face ="bold", size = 14, angle = 15),
        axis.text.y= element_text(size = 14))

# save.image("new.sim1_2.RData")
