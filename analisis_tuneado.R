# CON BASE TUNEADA

library(Rmisc)
library(plyr)
library(ggplot2)
library(psych)
library(glmmTMB)
library(readxl)
library(DHARMa)
library(emmeans)
library(car)

setwd("C:/Users/Admin/Desktop/GAB/GAB 2024/JJB")
posgrado=read.csv("BD_tuneada.csv", header = T, stringsAsFactors = T,sep = ";")

posgrado$ensayo<-as.factor(posgrado$ensayo)
posgrado$Replica<-as.factor(posgrado$Replica)


## DESCRIPTIVA

summary(posgrado)
str(posgrado)

## Descriptiva "a lo bruto", sin considerar v. aleatorias

# Hago la descriptiva de cant de parasitos, lo otro ya esta en el script de grado

estad <- summarySE(posgrado, measurevar="parasitos", groupvars=c("cepa"))
estad

# grafico 1 #quedo rari
ggplot(estad, aes(x=cepa, y=parasitos, colour=cepa)) + 
  stat_summary(fun.y = "mean", colour = "darkblue", size = 6, geom = "point") +
  geom_errorbar(aes(ymin=parasitos-sd, ymax=parasitos+sd),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) + ylab("Cant Parasitos")

# grafico 2 -> quedo mejor que el anterior

ggplot(posgrado, aes(x = cepa, y = parasitos,colour=cepa)) +  
  geom_boxplot(fill=24) + geom_jitter() + 
  labs(title="Box Plot",
       #caption="Source: mtcars",
       x="Cepa",
       y="Cant Parasitos")+
  theme_light()


######## MODELOS 

# Modelo a): Anova como LM (parasitos)

ma=glmmTMB(parasitos~cepa, posgrado)
summary(ma)
drop1(ma, test= "Chisq")

# comparaciones

# Se observaron diferencias en todos los pares, salvo en:
# RH-TgMr//RH-TgSb//RH-VEG//TgMr-TgSb//TgSb-VEG

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))

emmeans(ma,
        list(pairwise ~cepa ),
        type = "response")

### SUPUESTOS

simulationOutput <- simulateResiduals(fittedModel = ma, refit = FALSE, n=1000)
plot(simulationOutput)  #estupidas pruebas no parametricas
testCategorical(simulationOutput, catPred = posgrado$cepa)

## 

e<-resid(ma) # residuos de pearson
pre<-predict(ma) #predichos
par(mfrow = c(1, 2))
plot(pre, e, xlab="Predichos", ylab="Residuos de pearson",main="Grafico de dispersi?n de RE vs PRED",cex.main=.8 )
abline(0,0) ## feuchi el grafico
library(car)
qqPlot(e)
shapiro.test(e)
par(mfrow = c(1, 1))

######################################

# Modelo b) Promedio de VR

base2 <- summarySE(posgrado, measurevar="parasitos", groupvars=c("ensayo"))
base2=as.data.frame(base2)
base3=base2[,1:3]
cepas <- rep(c("RH", "ME49", "VEG", "TgMr", "TgSb"), each = 3)

# Agrega la columna 'cepa' al dataframe
base3$cepa <- cepas

# Entonces, el modelo b quedaria

mb=glmmTMB(parasitos~cepa, base3)
Anova(mb) 
summary(mb)
drop1(mb, test= "Chisq") #cepaTgSb no presenta diferencias

# comparaciones # solo se ven diferencias entre ME49 - TgMr

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))

emmeans(mb,
        list(pairwise ~cepa ),
        type = "response")

### SUPUESTOS

simulationOutput <- simulateResiduals(fittedModel = mb, refit = FALSE, n=1000)
plot(simulationOutput)  
testCategorical(simulationOutput, catPred = base3$cepa)

# 

e<-resid(mb) # residuos de pearson
pre<-predict(mb) #predichos
par(mfrow = c(1, 2))
plot(pre, e, xlab="Predichos", ylab="Residuos de pearson",main="Grafico de dispersi?n de RE vs PRED",cex.main=.8 )
abline(0,0) #feisimo el grafico

library(car)
qqPlot(e)
shapiro.test(e)
par(mfrow = c(1, 1))
# horribles!

######################################

# Modelo c) Poisson mixta

mc=glmmTMB(parasitos~cepa+(1|ensayo), posgrado, family = poisson)
Anova(mc)

# comparaciones 
# solo se ven diferencias entre: ME49-TgMr // TgMr-TgSb // ME49-RH

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))
emmeans(mc,
        list(pairwise ~cepa ),
        type = "response")

### SUPUESTOS

simulationOutput <- simulateResiduals(fittedModel = mc, refit = FALSE, n=1000)
plot(simulationOutput)  
testCategorical(simulationOutput, catPred = posgrado$cepa)
