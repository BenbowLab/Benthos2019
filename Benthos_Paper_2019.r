#########################################################################################
###                             BENTHOS PAPER 2019                                    ###
#########################################################################################

# Set working directory
setwd("~/PHD Torino 2015/II anno 2016-17/Abroad period USA 2017/Dati")

# Upload data
data<-read.table("Transect2014All_AD_09.10.2017.csv", dec=",", sep=";", header=TRUE)

# Packages and libraries
library(MASS) 
library(vegan)
library(indicspecies)
library(ggplot2)
library(lattice)
library(ggpubr)
source("HighstatLibV9.R")
library(gplots)
install.packages("betapart")
library(betapart)
install.packages("gridExtra")
library(gridExtra)
library(lme4)



### ADIPART -------------------------------------------------------

# Community data
y<-read.table("Comunit?.csv", dec=",", sep=";", header=TRUE)

# spatial scales
x<-read.table("Levels.csv", dec=",", sep=";", header=TRUE)
x$Transect<-as.factor(x$Transect)
x$TGr<-as.factor(x$TGr)
x$Stream<-as.factor(x$Stream)

str(x)

# Taxa richness
adipart(y, x, index="richness", relative=TRUE, weights="unif", nsimul=999)

# graph
Components<-c("Alpha", "Beta1", "Beta2", "Beta3", "Alpha", "Beta1", "Beta2", "Beta3")
OA<-c("Observed","Observed","Observed","Observed","Predicted","Predicted","Predicted","Predicted")
Valori<-c(0.184928, 0.095933, 0.282536, 0.436603,0.25066, 0.11330, 0.26609, 0.36995)
Dati<-data.frame(Components, OA, Valori)
Dati$Valori<-Dati$Valori*100
Colori<-c("white", "grey80", "grey50", "black")

S <- ggplot() + geom_bar(aes(y = Valori, x = OA, fill = forcats::fct_rev(Components)), colour="black",
                         data = Dati, stat="identity") +
  labs(x="", y="%", fill="Components") +
  scale_fill_manual(values = Colori) +
  theme_classic()
  

# EPT richness
y2<-read.table("Comunit?EPT.csv", dec=",", sep=";", header=TRUE)

data<-read.table("Transect2014All_AD_09.10.2017.csv", dec=",", sep=";", header=TRUE)
plot(data$EPTrichness)
identify(data$EPTrichness)
x2<-x[-c(1,217),]

adipart(y2, x2, index="richness", relative=TRUE, weights="unif", nsimul=999)

# graph
Components<-c("Alpha", "Beta1", "Beta2", "Beta3", "Alpha", "Beta1", "Beta2", "Beta3")
OA<-c("Observed","Observed","Observed","Observed","Predicted","Predicted","Predicted","Predicted")
Valori<-c(0.24191, 0.11120, 0.29282, 0.35407, 0.33414, 0.13328, 0.25635, 0.27623)
Dati<-data.frame(Components, OA, Valori)
Dati$Valori<-Dati$Valori*100
Colori<-c("white", "grey80", "grey50", "black")

EPTS <- ggplot() + geom_bar(aes(y = Valori, x = OA, fill = forcats::fct_rev(Components)), colour="black",
                         data = Dati, stat="identity") +
  labs(x="", y="%", fill="Components") +
  scale_fill_manual(values = Colori) +
  theme_classic()


# Shannon
adipart(y, x, index="shannon", relative=TRUE, weights="unif", nsimul=999)

# graph
Components<-c("Alpha", "Beta1", "Beta2", "Beta3", "Alpha", "Beta1", "Beta2", "Beta3")
OA<-c("Observed","Observed","Observed","Observed","Predicted","Predicted","Predicted","Predicted")
Valori<-c(0.58542, 0.13621, 0.18402, 0.09435, 0.705815, 0.128720, 0.131238, 0.034227)
Dati<-data.frame(Components, OA, Valori)
Dati$Valori<-Dati$Valori*100
Colori<-c("white", "grey80", "grey50", "black")

H <- ggplot() + geom_bar(aes(y = Valori, x = OA, fill = forcats::fct_rev(Components)), colour="black",
                         data = Dati, stat="identity") +
  labs(x="", y="%", fill="Components") +
  scale_fill_manual(values = Colori) +
  theme_classic()


# Simpson
adipart(y, x, index="simpson", relative=TRUE, weights="unif", nsimul=999)

# graph
Components<-c("Alpha", "Beta1", "Beta2", "Beta3", "Alpha", "Beta1", "Beta2", "Beta3")
OA<-c("Observed","Observed","Observed","Observed","Predicted","Predicted","Predicted","Predicted")
Valori<-c(0.815526, 0.089035, 0.069619, 0.025820, 0.8886521, 0.0678334, 0.0377477, 0.0057668)
Dati<-data.frame(Components, OA, Valori)
Dati$Valori<-Dati$Valori*100
Colori<-c("white", "grey80", "grey50", "black")

D <- ggplot() + geom_bar(aes(y = Valori, x = OA, fill = forcats::fct_rev(Components)), colour="black",
                         data = Dati, stat="identity") +
  labs(x="", y="%", fill="Components") +
  scale_fill_manual(values = Colori) +
  theme_classic()

tiff("Adipart.tiff",height = 14, width = 17, 
     units = 'cm',res=300)
ggarrange(S, EPTS, H, D, 
          labels = c("a", "b", "c", "d"),
          ncol = 2, nrow = 2)
dev.off()


### NMDS -----------------------------------------------------------------------

simbols<-c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14",
           "T15","T16","T17","T18","T19","T20","T21","T22")

community<-read.table("Comunit?.csv", dec=",", sep=";", header=TRUE)
community<-read.table("Comunit?2.csv", dec=",", sep=";", header=TRUE)

beta_dist <- vegdist(community, index = "bray")

mds <- metaMDS(beta_dist)

mds_data <- as.data.frame(mds$points)

mds_data$SampleID <- rownames(mds_data)
#mds_data<-data.frame(mds_data,data$TGr)
mds_data<-data.frame(mds_data,simbols)
#mds_data$data.TGr<-as.factor(mds_data$data.TGr)
mds_data$simbols<-as.factor(mds_data$simbols)
colnames(mds_data)<-c("NMDS1","NMDS2","SampleID","TG")

tiff("NMDS_TG.tiff",height = 12, width = 17, 
     units = 'cm',res=300)
ggplot(mds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed",color = "black", size=0.5) +
  geom_point(size=3) +
  #ylim(-0.28,0.28)+
  geom_text(data=mds_data, aes(x=NMDS1,y=NMDS2,label=SampleID),size=4,vjust=-0.5)+
  theme_classic() +
  theme(panel.border = element_rect(linetype = 1, fill = NA))
dev.off()



### Beta-diversity (Baselga et al. 2013) ----------------------------------------------------

community<-read.table("Comunit?2.csv", dec=",", sep=";", header=TRUE)
community<-y

A<-bray.part(community)

BRAY <- as.matrix(A$bray)
BRAY.BAL<-as.matrix(A$bray.bal)
BRAY.GRA<-as.matrix(A$bray.gra)

mean(BRAY)
mean(BRAY.BAL)
mean(BRAY.GRA)

range(BRAY)
range(BRAY.BAL)
range(BRAY.GRA)

sd(BRAY)
sd(BRAY.BAL)
sd(BRAY.GRA)

d <- data.frame(x = unlist(A), 
                grp = rep(letters[1:length(A)],times = sapply(A,length)))
ggplot(d,aes(x = grp, y = x)) + geom_boxplot() + theme_classic() +
  labs(title="", y="Bray-Curtis dissimilarity", x="Components of beta-diversity")+
  scale_x_discrete(labels=c("BC_bal","BC_gra","BC"))

Components = c("BC.bal", "BC.gra", "Bray")


### GLMMs e Semivariograms ---------------------------------------------------------------

data<-read.table("Transect2014All_AD_09.10.2017.csv", dec=",", sep=";", header=TRUE)

y<-read.table("Comunit?.csv", dec=",", sep=";", header=TRUE)

data$H<-diversity(y, index="shannon") # add shannon-wiener index 
data$D<-diversity(y, index="simpson") # add simpson's index

data<-na.omit(data)

# Outliers in Y variables?
dotchart(data$Richness, main="Richness") # No outliers
dotchart(data$Ntot, main="Total abundance") # No outliers
dotchart(data$EPTrichness, main="ETP richness") # No outliers
dotchart(data$EPTN, main="ETP abundance") # Similar to Ntot but no outliers 
dotchart(data$H, main="H") # 4 zero, No outliers 
dotchart(data$D, main="D")c# 4 zero

par(mfrow=c(1,1))

max(data$Richness)
min(data$Richness)
mean(data$Richness)
sd(data$Richness)
max(data$Ntot)
min(data$Ntot)
mean(data$Ntot)
sd(data$Ntot)

# Outliers in X variables?
par(mfrow=c(2,2))
dotchart(data$Depth, main="Depth") # No outliers)
table(data$Depth_Bin) # quite balanced despite 0-0.5 has half number observations
dotchart(data$Velocity, main="Velocity") # No outliers
table(data$VelocityBin) # quite balanced with excepiton for 0-0,25 and 1-1,25
dotchart(data$Froude, main="Froude Number") # No outliers
dotchart(data$St_Larvae, main="Number of Sturgeon larvae") # a lot of zero values
table(data$YesNo) # very unbalanced 200 No, 20 Yes !
table(data$Cover_Class) # highly heterogeneous ---> further clustering ?
dotchart(data$NCoverCL)
table(data$TGr)
dotchart(data$Hsub) # no outliers
dotchart(data$SubInd) # no outliers
dotchart(data$Avg_Rock_Size, main="Avg_Rock_Size") # 3 outliers
dotchart(data$Med_Rock_Size, main="Med_Rock_Size") # 2 outliers
dotchart(data$med_RS, main="median_Rock_Size_New") # 1 outlier
dotchart(data$Coeff_Var, main="Coeff of Variation") # no outliers

max(data$Froude)
min(data$Froude)
mean(data$Froude)
max(data2$Med_Rock_Size)
min(data2$Med_Rock_Size)
mean(data2$Med_Rock_Size)

max(data2$Froude)
min(data2$Froude)
mean(data2$Froude)
max(data$med_RS)
min(data$med_RS)
mean(data$med_RS)

plot(data$med_RS)
identify(data$med_RS)
data<-data[-c(90),]
# eliminati gli oultiers

# Collinearity X variables (covariates)

tiff("Collinarity Matrix.tiff",height = 12, width = 17, 
     units = 'cm',res=300)
par(mfrow=c(1,1))
MyVar <- c("Depth","Velocity","Froude","med_SS", "Coeff_Var","Fine")
Mypairs(data[,MyVar])
dev.off()

library(GGally)

Var<-data.frame(data$Depth, data$Velocity, data$Froude, data$med_RS, data$Coeff_Var,
                data$Fine)
colnames(Var) <- c("Depth (m)","Velocity (m/s)","Froude","Substrate size (cm)", "Coeff. Variation","Fine (%)")

tiff("Collinarity Matrix.tiff",height = 13, width = 17, 
     units = 'cm',res=300)
ggpairs(Var)+theme_bw()+theme(text=element_text(size=7))
dev.off()


## Semivariograms with gstat and sp --------------------------

install.packages("gstat")
library(gstat)  
install.packages("sp")
library(sp)


data2<-data.frame(data$TransectX, data$Random, data$Richness, data$Ntot, 
                  data$EPTrichness, data$EPTN, data$H, data$D)
colnames(data2)<-c("x","y","Richness","Ntot","EPTrichness","EPTN","H","D")
data2$LOGS<-log(data2$Richness+1)
data2$LOGN<-log(data2$Ntot+1)
data2$LOGEPTS<-log(data2$EPTrichness+1)
data2$LOGEPTN<-log(data2$EPTN+1)
data2$LOGH<-log(data2$H+1)
data2$LOGD<-log(data2$D+1)

data3<-data.frame(data$x, data$Random, data$Richness, data$Ntot, 
                  data$EPTrichness, data$EPTN, data$H, data$D)
colnames(data3)<-c("x","y","Richness","Ntot","EPTrichness","EPTN", "H","D")
data3$LOGS<-log(data3$Richness+1)
data3$LOGN<-log(data3$Ntot+1)
data3$LOGEPTS<-log(data3$EPTrichness+1)
data3$LOGEPTN<-log(data3$EPTN+1)
data3$LOGH<-log(data3$H+1)
data3$LOGD<-log(data3$D+1)


coordinates(data2)= ~ x+y # usare coordinate relative per plottare
coordinates(data3)= ~ x+y # usare distanza per semivariograms

SB<-bubble(data2, zcol="Richness", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")
NB<-bubble(data2, zcol="Ntot", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")
EPTSB<-bubble(data2, zcol="EPTrichness", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")
EPTNB<-bubble(data2, zcol="EPTN", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")
HB<-bubble(data2, zcol="H", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")
DB<-bubble(data2, zcol="D", fill=TRUE, do.sqrt=FALSE, maxsize=1, col="black")


tiff("Bubbles.tiff",height = 17, width = 17, 
     units = 'cm',res=300)
grid.arrange(SB, NB, EPTSB, EPTNB,H, D, nrow=6, ncol=1)
dev.off()


S_var<-variogram(LOGS~1, data=data3)
S<-plot(S_var, col="black", pch=19, xlab="distance (m)", main="Total richness", ylim=c(0,0.25))

N_var<-variogram(LOGN~1, data=data3)
N<-plot(N_var, col="black", pch=19, xlab="distance (m)", main="Total abundance", ylim=c(0,0.7))

EPTS_var<-variogram(LOGEPTS~1, data=data3)
EPTS<-plot(EPTS_var, col="black", pch=19, xlab="distance (m)", main="EPT richness",ylim=c(0,0.25))

EPTN_var<-variogram(LOGEPTN~1, data=data3)
EPTN<-plot(EPTN_var, col="black", pch=19, xlab="distance (m)", main="EPT abundance",ylim=c(0,0.8))


H_var<-variogram(LOGH~1, data=data3)
H<-plot(H_var, col="black", pch=19, xlab="distance (m)", main="Shannon-Wiener index", ylim=c(0,0.06))

D_var<-variogram(LOGD~1, data=data3)
D<-plot(D_var, col="black", pch=19, xlab="distance (m)", main="Simpson's index", ylim=c(0,0.015))

TheVariogramModel <- vgm(psill=0.6, model="Exp", nugget=0, range=200)
plot(EPTN_var, model=TheVariogramModel)

tiff("Semivariograms.tiff",height = 15, width = 17, 
     units = 'cm',res=300)
grid.arrange(S, N, EPTS, EPTN,H, D, nrow=3, ncol=2)
dev.off()


### Models (GLMMs) ---------------------------------------

data$TransectX<-as.factor(data$TransectX)

# Taxa richness #
M1<-glmer(Richness ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), family = poisson, data = data)
summary(M1)

E1 <- resid(M1, type = "pearson")
N  <- nrow(data2)
p  <- length(fixef(M1)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion
# No overdispersion

drop1(M1)
M1<-glmer(Richness ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), family = poisson, data = data)
summary(M1)


# Froude
M1<-glm(Richness ~ Froude, family=poisson, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

TS=ggplot() + 
  xlab("Froude number") +
  ylab("Taxa richness")+
  geom_point(data = data, aes(y = Richness, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
TS


# Rock size
M1<-glm(Richness ~ med_RS, family=poisson, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

TS2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("Taxa richness")+
  geom_point(data = data, aes(y = Richness, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
TS2


# Total abundance #
M1<-glmer.nb(Ntot ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

E1 <- resid(M1, type = "pearson")
N  <- nrow(data2)
p  <- length(fixef(M1)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion
# No overdispersion

drop1(M1)
M1<-glmer.nb(Ntot ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

# Froude
M1<-glm(Ntot ~ Froude, family=poisson, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

NT=ggplot() + 
  xlab("Froude number") +
  ylab("Total abundance")+
  geom_point(data = data, aes(y = Ntot, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
NT


# Rock size
M1<-glm(Ntot ~ med_RS, family=poisson, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

NT2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("Total abundance")+
  geom_point(data = data, aes(y = Ntot, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
NT2


# EPT richness #
M1<-glmer(EPTrichness ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), family = poisson, data = data)
summary(M1)

E1 <- resid(M1, type = "pearson")
N  <- nrow(data2)
p  <- length(fixef(M1)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion
# No overdispersion

drop1(M1)
M1<-glmer(EPTrichness ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), family = poisson, data = data)
summary(M1)

# Froude
M1<-glm(EPTrichness ~ Froude, family=poisson, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

EPTS=ggplot() + 
  xlab("Froude number") +
  ylab("EPT richness")+
  geom_point(data = data, aes(y = EPTrichness, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
EPTS


# Rock size
M1<-glm(EPTrichness ~ med_RS, family=poisson, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

EPTS2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("EPT richness")+
  geom_point(data = data, aes(y = EPTrichness, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
EPTS2


# EPT abundance #
M1<-glmer.nb(EPTN ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

E1 <- resid(M1, type = "pearson")
N  <- nrow(data2)
p  <- length(fixef(M1)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion
# No overdispersion

drop1(M1)
M1<-glmer.nb(EPTN ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

# Froude
M1<-glm(EPTN ~ Froude, family=poisson, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

EPTN=ggplot() + 
  xlab("Froude number") +
  ylab("EPT abundance")+
  geom_point(data = data, aes(y = EPTN, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
EPTN


# Rock size
M1<-glm(EPTN ~ med_RS, family=poisson, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

EPTN2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("EPT abundance")+
  geom_point(data = data, aes(y = EPTN, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
EPTN2


# Shannon-Wiener #

hist(data$H)

install.packages("lmerTest")
library(lmerTest)

M1<-lmer(H ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

drop1(M1)
M1<-lmer(H ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)


# Froude
M1<-lm(H ~ Froude, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

H=ggplot() + 
  xlab("Froude number") +
  ylab("Shannon-Wiener")+
  geom_point(data = data, aes(y = H, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
H


# Rock size
M1<-lm(H ~ med_RS, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

H2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("Shannon-Wiener")+
  geom_point(data = data, aes(y = H, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
H2


# Simpson #

plot(data$D)
identify(data$D)
data<-data[-c(1,2,146,162,178),]
# oultiers removed

hist(data$D)

M1<-glmer(D ~ Froude*med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)

drop1(M1)
M1<-glmer(D ~ Froude + med_RS + (1|TransectX) + (1|TransectGroup), data = data)
summary(M1)


# Froude
M1<-lm(D ~ Froude, data=data)
MyData1 <- data.frame(Froude =
                        seq(from = min(data$Froude),
                            to = max(data$Froude),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

D=ggplot() + 
  xlab("Froude number") +
  ylab("Simpson")+
  geom_point(data = data, aes(y = D, x = Froude),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = Froude, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = Froude, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = Froude, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
D


# Rock size
M1<-lm(D ~ med_RS, data=data)
MyData1 <- data.frame(med_RS =
                        seq(from = min(data$med_RS),
                            to = max(data$med_RS),
                            length=10))

P1 <- predict(M1, newdata = MyData1, type = "response", se=T)

MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

D2=ggplot() + 
  xlab("Substrate size (cm)") +
  ylab("Simpson")+
  geom_point(data = data, aes(y = D, x = med_RS),shape = 16, size = 1, col="grey50") +
  geom_line(data = MyData1,aes(x = med_RS, y = mu), size=1, col="black") +
  geom_line(data = MyData1,aes(x = med_RS, y = selow), size=0.5, col="black", linetype="dashed")+
  geom_line(data = MyData1,aes(x = med_RS, y = seup), size=0.5, col="black", linetype="dashed")+
  theme_classic()+theme(text=element_text(size=7))
D2



tiff("Models.tiff",height = 14, width = 17, 
     units = 'cm',res=300)
ggarrange(TS, TS2, NT, NT2,
          EPTS, EPTS2,EPTN, EPTN2,
          H, H2, D, D2,
          labels = c("a", "b", "c", "d",
                     "e", "f", "g", "h",
                     "i", "j", "k", "l"),
          ncol = 4, nrow = 3)
dev.off()


tiff("Models_S_e_N.tiff",height = 14, width = 17, 
     units = 'cm',res=300)
ggarrange(TS, TS2, NT, NT2,
          labels = c("a", "b", "c", "d"),
          ncol = 2, nrow = 2)
dev.off()


### Taxa accumulation curve --------------------------------------------------------------

y<-read.table("Comunit?.csv", dec=",", sep=";", header=TRUE)

SP<-specaccum(y, method = "random", permutations = 100,
              conditioned =TRUE)
tiff("Taxa accumulation curve.tiff",height = 11, width = 17, 
     units = 'cm',res=300)
plot(SP, ci.type="polygon", ci.col="grey90", xlab="Samples", ylab="Taxon richness", 
     xlim=c(0,220), main="", las=1)
dev.off()


### OMI ------------------------------------------------------------------------------------

Env3<-read.table("Env3.csv", dec=",", sep=";", header=TRUE)

Y<-read.table("Community.csv", dec=",", sep=";", header=TRUE)

library(ade4)

pca3<-dudi.pca(Env3, scale = TRUE, scan=FALSE, nf=2)
nic3 <- niche(pca3, log(Y + 1), scan = FALSE)
plot(nic3)
niche.param(nic3)
rtest(nic3,999)

{
  par(mfrow = c(2, 2))
  s.traject(pca3$li, clab = 0)
  s.traject(nic3$ls, clab = 0)
  s.corcircle(nic3$as)
  s.arrow(nic3$c1)
  par(mfrow = c(8,5))
  s.arrow(nic3$c1)
  s.arrow(nic3$li, clab = 0.7)
  for (i in 1:38) s.distri(nic3$ls, as.data.frame(Y[,i]),
                           csub = 2, sub = names(Y)[i], cstar=0)
}

par(mfrow=c(2,2))
s.corcircle(nic3$as)
s.arrow(nic3$c1, clab=1)
s.arrow(nic3$li, clab = 0.8)
s.distri(nic3$ls,dfdistri = Y, xax=1, yax=2, cellipse=1.5, 
         axesell=F, cstar=0, label = names(Y))

# OMI Parameters
pca3<-dudi.pca(Env3, scale = TRUE, scan=FALSE, nf=6)
nic3 <- niche(pca3, log(Y + 1), scan = FALSE, nf=6)

# Eigenvalues
pca3$eig
sum(pca3$eig)
barplot(pca3$eig)
(kip <- 100 * pca3$eig/sum(pca3$eig))
cumsum(kip)
# Variance expplained by the axes: 
# 50.2558620 19.9104526 13.2653906 10.6282905  5.0797349  0.8602696

# Scores axes
nic3$c1
sum(pca3$cw * pca3$c1$CS1^2)

# Scores species
nic3$l1

#Significant taxa
Y2<-data.frame(Y$Isonychiidae,
               Y$Heptageniidae,
               Y$Perlidae,
               Y$Ephemerellidae,
               Y$Hydropsychidae,
               Y$Lepidostomatidae,
               Y$Helicopsychidae,
               Y$Leptoceridae,
               Y$Elmidae,
               Y$Brachycentridae,
               Y$Glossosomatidae,
               Y$Simuliidae,
               Y$Odontoceridae)

colnames(Y2)<-c("Isonychiidae", "Heptageniidae", "Perlidae", "Ephemerellidae",
                "Hydropsychidae", "Lepidostomatidae", "Helicopsychidae",
                "Leptoceridae", "Elmidae", "Brachycentridae",
                "Glossosomatidae", "Simuliidae", "Odontoceridae")

tiff("OMI3.tiff",height = 12, width = 17, 
     units = 'cm',res=300)
par(mfrow = c(4,4))
s.corcircle(nic3$as, clabel = 1.3)
s.arrow(nic3$c1, clabel = 1.3)
s.distri(nic3$ls,dfdistri = Y2, xax=1, yax=2, cellipse=1.5, 
         axesell=F, cstar=0, label = names(Y2), clabel=0, cpoint=0)
for (i in 1:13) s.distri(nic3$ls, as.data.frame(Y2[,i]), axesell=F,
                         csub = 2, sub = names(Y2)[i], cstar=0, cpoint=0)
dev.off()



