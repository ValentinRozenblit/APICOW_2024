# Libraries needed for this analysis
library(Rmisc)
library(plyr)
library(ggplot2)
library(psych)
library(glmmTMB)
library(readxl)
library(DHARMa)
library(emmeans)
library(car)

setwd("")
df=read.csv("database_APICOW_2024.csv", header = T, stringsAsFactors = T,sep = ";")

df$experiment<-as.factor(df$experiment)
df$rep<-as.factor(df$rep)


## DESCRIPTIVE ANALYSIS

summary(df)
str(df)

## First Descriptive approach, withour taking random effects into account 

estad <- summarySE(df, measurevar="n_parasites", groupvars=c("strain"))
estad

# First plot
ggplot(estad, aes(x=strain, y=n_parasites, colour=strain)) + 
  stat_summary(fun.y = "mean", colour = "darkblue", size = 6, geom = "point") +
  geom_errorbar(aes(ymin=n_parasites-sd, ymax=n_parasites+sd),
                width=.2, # Width of the error bars
                position=position_dodge(.9)) + ylab("Number of Parasites")

# Second plot (better than the first one)

ggplot(df, aes(x = strain, y = n_parasites,colour=strain)) +  
  geom_boxplot(fill=24) + geom_jitter() + 
  labs(title="Box Plot",
       #caption="Source: mtcars",
       x="strain",
       y="Number of Parasites")+
  theme_light()


######## STATISTICAL MODELS

# Model a): Anova as a LM (n_parasites)

ma=glmmTMB(n_parasites~strain, df)
summary(ma)
drop1(ma, test= "Chisq")

# Comparisons

# Differences were detected for every pair except for:
# RH-TgMr//RH-TgSb//RH-VEG//TgMr-TgSb//TgSb-VEG

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))

emmeans(ma,
        list(pairwise ~strain ),
        type = "response")

### ASSUMPTIONS 

simulationOutput <- simulateResiduals(fittedModel = ma, refit = FALSE, n=1000)
plot(simulationOutput)  #estupidas pruebas no parametricas
testCategorical(simulationOutput, catPred = df$strain)

## 

e<-resid(ma) # pearson residuals
pre<-predict(ma) # predicted
par(mfrow = c(1, 2))
plot(pre, e, xlab="Predicted", ylab="Pearson residuals",main="Scatterplot RE vs PRED",cex.main=.8 )
abline(0,0) ## feuchi el grafico
qqPlot(e)
shapiro.test(e)
par(mfrow = c(1, 1))

######################################

# Model b) Averaging the response variable

df2 <- summarySE(df, measurevar="n_parasites", groupvars=c("experiment"))
df2=as.data.frame(df2)
base3=df2[,1:3]
strains <- rep(c("RH", "ME49", "VEG", "TgMr", "TgSb"), each = 3)

# Adds 'strain' column to dataframe
base3$strain <- strains

# Hence, model b is stated as follows:
mb=glmmTMB(n_parasites~strain, base3)
Anova(mb) 
summary(mb)
drop1(mb, test= "Chisq") #strainTgSb does not show differences

# Comparisons. Differences were observed for ME49 - TgMr only.

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))

emmeans(mb,
        list(pairwise ~strain ),
        type = "response")

### Assumptions

simulationOutput <- simulateResiduals(fittedModel = mb, refit = FALSE, n=1000)
plot(simulationOutput)  
testCategorical(simulationOutput, catPred = base3$strain)

# 

e<-resid(mb) # Pearson residuals
pre<-predict(mb) # predicted
par(mfrow = c(1, 2))
plot(pre, e, xlab="Predicted", ylab="Pearson residuals",main="Scatterplot RE vs PRED",cex.main=.8 )
abline(0,0) #feisimo el grafico

qqPlot(e)
shapiro.test(e)
par(mfrow = c(1, 1))

# Assumptions were not met. 

######################################

# Model c) Mixed-effects Poisson

mc=glmmTMB(n_parasites~strain+(1|experiment), df, family = poisson)
Anova(mc)

# Comparisons 
# Differnces were observed for: ME49-TgMr // TgMr-TgSb // ME49-RH only.

options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),
                      contrast = list(infer = c(TRUE, TRUE))))
emmeans(mc,
        list(pairwise ~strain ),
        type = "response")

### ASSUMPTIONS

simulationOutput <- simulateResiduals(fittedModel = mc, refit = FALSE, n=1000)
plot(simulationOutput)  
testCategorical(simulationOutput, catPred = df$strain)
