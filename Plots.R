library(ggplot2)
library(grid)
library(scales)
library(boot)
library(dplyr)
library(tidyverse)
library(MASS)
library(lme4)
library(sjPlot)
library(data.table)
library(cli) 
library(spatstat)
library(rstudioapi)
library(plyr)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

#Pixels to degrees
source("pix2deg.R")
winx <- 1920  
winy <- 1080
distCM<- 70
moWidthCM<- 53.5
px2deg <- pix2deg(1, winx, moWidthCM, distCM)

xmargin=(winx-1500)/2
ymargin=(winy-1000)/2

xrange=c(xmargin*px2deg, (winx-xmargin)*px2deg)
yrange=c(ymargin*px2deg, (winy-ymargin)*px2deg)

win <- owin(xrange,yrange) # Observation window

# Load data
Filename1<-sprintf("RAW_mm.RData")
load(Filename1)

Filename2<-sprintf("SAC_mm.RData")
load(Filename2)

df <- RAW %>%
  dplyr::filter(., startEvent==1 & sacDP==1 &
                  x>xrange[1] & x<xrange[2] &
                  y>yrange[1] & y<yrange[2]) %>%
  dplyr::filter(fixdur<1000)%>% # SET MAX FIXDUR
  dplyr::mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days")) %>%
  dplyr::mutate(Viewing= case_when(
    viewnr==1 ~ "First",
    viewnr==2 ~ "Second")) %>%
  dplyr::mutate(Condition=as.factor(Condition))%>%
  dplyr::mutate(Viewing=as.factor(Viewing))

df=subset(df, df$imtype2 != "Single(BaselineSess2)")  
df$VP <- factor(df$VP)
df$Img <- factor(df$Img)

n <- length(df$Img)
cond <- rep(0, n)
df$cond = cond

df <- df %>%
  dplyr::mutate(cond= case_when(
    Condition=="Minutes" & Viewing=="First" ~ "First viewing - Minutes Condition",
    Condition=="Minutes" & Viewing=="Second" ~ "Second viewing - Minutes Condition",
    Condition=="Days" & Viewing=="First" ~ "First viewing - Days Condition",
    Condition=="Days" & Viewing=="Second" ~ "Second viewing - Days Condition"))

df <- df %>%
  dplyr::select(VP,Img,x,y,fixdur,Condition,Viewing,cond) 


# ________________________Fixation Duration_____________________________

# Distribution of fixation durations
distr<-ggplot(df, aes(x=fixdur, color=cond)) +
  geom_freqpoly(aes(color = cond),
                bins = 25, size = 1.5)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#868686FF")) +
  labs(colour="cond") +
  geom_vline(data=mu, aes(xintercept=grp.mean,color=cond),
             linetype = "longdash", size=1)
print(distr + labs(y="Frequency", x = "Fixation Duration (ms)"))

# Barplot
mu <- ddply(df, "cond", summarise, grp.mean=mean(fixdur))
means.sem <- ddply(df, c("Condition","Viewing"), summarise,
                   mean=mean(fixdur), sem=sd(fixdur)/sqrt(length(fixdur)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)


barplot<- ggplot(means.sem, aes(x=Condition,y=mean,fill=Viewing,width=0.5,space=c(0.2,0.2,0.2,0.2))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymax=upper,
                    ymin=lower),
                data=means.sem,stat = "identity",position = "dodge",width=0.5)

print(barplot + labs(y="Mean Fixation Duration", x = "Condition"))

# Error bar plot
bar_err <- ggplot(means.sem,aes(x=Viewing,y=mean,colour=Condition,group=Condition))+
  geom_point()+geom_line()+ ggtitle("Fixation Duration") +
  geom_errorbar(aes(ymax=upper,
                    ymin=lower),
                data=means.sem,colour="grey", width=.1)
print(bar_err + labs(y="Mean Fixation Duration", x = "Viewing"))


# ________________________Saccade Amplitude_____________________________

df2 <- SAC %>%
  dplyr::filter(., nth != 1 & imtype2 != "Single(BaselineSess2)"&
                  blinkSac!= TRUE & blinkFix != TRUE &
                  x>xrange[1] & x<xrange[2] &
                  y>yrange[1] & y<yrange[2]) %>%
  dplyr::select(.,x,y,VP,viewnr,PeakVel,fixdur,Ampl,session,imtype2,Img, trial)%>%
  mutate(DV=Ampl)%>%
  mutate(DV_log=log(DV))

df2<-df2 %>% 
  mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days"
  )) %>%
  mutate(Viewing= case_when(
    viewnr=="Viewing 1" ~ "First",
    viewnr=="Viewing 2" ~ "Second"
  )) %>%
  mutate(VP=as.factor(VP),Img=as.factor(Img), DV_reci= 1/DV) %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" )) %>%
  mutate(Viewing = fct_relevel(Viewing, "First" ,"Second" )) %>%
  mutate(Condition=as.factor(Condition))%>%
  mutate(Viewing=as.factor(Viewing))

n <- length(df2$Img)
Factor <- rep(0, n)
df2$Factor = Factor

df2 <- df2 %>%
  dplyr::mutate(Factor= case_when(
    Condition=="Minutes" & Viewing=="First" ~ "First viewing - Minutes Condition",
    Condition=="Minutes" & Viewing=="Second" ~ "Second viewing - Minutes Condition",
    Condition=="Days" & Viewing=="First" ~ "First viewing - Days Condition",
    Condition=="Days" & Viewing=="Second" ~ "Second viewing - Days Condition")) %>%
  mutate(Factor=as.factor(Factor)) %>% 
  mutate(Factor = fct_relevel(Factor, "First viewing - Minutes Condition" ,"Second viewing - Minutes Condition", "First viewing - Days Condition", "Second viewing - Days Condition" )) %>% 
  mutate(Factor=as.factor(Factor))
head(df2)

# Distribution of saccades
distr_sac<-ggplot(df2, aes(x=Ampl, color=Factor)) +
  geom_freqpoly(aes(color = Factor),
                bins = 25, size = 1.5)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#868686FF")) +
  labs(colour="Condition") +
  geom_vline(data=mu2, aes(xintercept=grp.mean,color=Factor),
             linetype = "longdash", size=1)
print(distr_sac + labs(y="Frequency", x = "Saccade Amplitude"))

# Barplot
mu2 <- ddply(df2, "Factor", summarise, grp.mean=mean(Ampl))
means.sem <- ddply(df2, c("Condition","Viewing"), summarise,
                   mean=mean(Ampl), sem=sd(Ampl)/sqrt(length(Ampl)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)

barplot_sac <- ggplot(means.sem, aes(x=Condition,y=mean,fill=Viewing,width=0.5,space=c(0.2,0.2,0.2,0.2))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymax=upper,
                    ymin=lower),
                data=means.sem,stat = "identity",position = "dodge",width=0.5)

print(barplot_sac + labs(y="Mean Saccade Amplitude", x = "Condition"))

# Error bar plot
errbar_sac <- ggplot(means.sem,aes(x=Viewing,y=mean,colour=Condition,group=Condition))+
  geom_point()+geom_line()+ ggtitle("Saccade Amplitude") +
  geom_errorbar(aes(ymax=upper,
                    ymin=lower),
                data=means.sem,colour="grey", width=.1)
print(errbar_sac + labs(y="Mean Saccade Amplitude", x = "Viewing"))


# Main Sequence
main_seq<-ggplot(df2, aes(Ampl, PeakVel)) + geom_point(size=0.5, color="blue") +
  geom_smooth(method=lm, se=FALSE, color="darkred")+
  xlab('Saccade Amplitude')+ylab('Peak Velocity') + ggtitle("Main Sequence") 
print(main_seq)


# ________________________Fatigue Effect_____________________________

meanfixtrial <- ddply(df2, c("Condition","trial"),summarise,mean=mean(fixdur))
meanfixtrial<-meanfixtrial %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" )) %>% 
  mutate_at(c('trial'), as.numeric)

(fat_fixdur<-ggplot(meanfixtrial, aes(x=trial, y=mean)) +
    geom_point(size=0.1) +
    geom_line() +
    geom_smooth()+
    labs(y="Fixation Duration (ms)", x = "Trial Number") +
    facet_grid(~Condition))

meanfixtrial <- ddply(df2, c("Condition","trial"),summarise,mean=mean(Ampl))
meanfixtrial<-meanfixtrial %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" ))%>% 
  mutate_at(c('trial'), as.numeric)

(fat_sacampl<-ggplot(meanfixtrial, aes(x=trial, y=mean)) +
    #geom_point(size=0.1) +
    geom_line() +
    geom_smooth()+
    labs(y="Saccade Amplitude (°)", x = "Trial Number") +
    facet_grid(~Condition))
