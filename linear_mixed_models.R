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
library(spatstat)
library(rstudioapi)
library(stargazer)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

# ________________________Fixation duration_____________________________

# Download and filter the data
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

# Load files
Filename1<-sprintf("RAW_mm.RData")
load(Filename1)

df_fixdur <- RAW %>%
  dplyr::mutate(startEvent= replace_na(startEvent, 0))%>%
  dplyr::filter(startEvent==1) %>%
  dplyr::filter(sacDP==1)%>%
  dplyr::filter(fixdur<1000)%>% # SET MAX FIXDUR
  dplyr::filter(fixdur>50)%>% # SET MIN FIXDUR
  dplyr::filter(.,imtype2 != "Single(BaselineSess2)" & x>xrange[1] & x<xrange[2] 
                & y>yrange[1] & y<yrange[2]) %>%
  dplyr::select(.,x,y,VP,viewnr,sacDP,fixdur,session,imtype2,Img)%>%
  dplyr::mutate(DV=fixdur)%>%
  dplyr::mutate(DV_log=log(DV))

df_fixdur<-df_fixdur %>% 
  mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days"
  )) %>%
  mutate(Viewing= case_when(
    viewnr==1 ~ "First",
    viewnr==2 ~ "Second"
  )) %>%
  mutate(Condition=as.factor(Condition),viewnr=as.factor(viewnr),VP=as.factor(VP),Img=as.factor(Img), DV_reci= 1/DV) %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" )) %>%
  mutate(Viewing = fct_relevel(Viewing, "First" ,"Second" ))

# set contrast
contrasts(df_fixdur$Viewing) <- contr.sdif(2) 
contrasts(df_fixdur$Condition) <- contr.sdif(2)

#check contrasts
contrasts(df_fixdur$Condition)
contrasts(df_fixdur$Viewing)

# Factors to numeric values
mat_cond <- model.matrix(~ Condition+Viewing+Condition:Viewing, df_fixdur)
df_fixdur<-data.frame(df_fixdur)
df_fixdur[, 15:18]<- mat_cond
names(df_fixdur)[15:18] <- c("Intercept", "n.Condition", "n.Viewing", "n.Interaction")

# Check if transformation is required
boxcox(DV ~ n.Condition + n.Viewing  + n.Interaction, data=df_fixdur) # lambda is close to 0  => log transform is required

# Log transform dependent variable
df_fixdur <- df_fixdur %>% 
  mutate(DV=DV_log)
labelDV<-sprintf("Log Fixation Duration")

# Visualize subject variance
Fcond_m1.lmm.subj <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + (1 | VP), data=df_fixdur, REML = FALSE)
str(Fcond_m1.re.subj <- ranef(Fcond_m1.lmm.subj))
str(Fcond_df_m1.re.subj <- as.data.frame(Fcond_m1.re.subj))

Fcond_df_m1.re.subj %>%  ggplot(aes(y=grp,x=condval)) +
  geom_point(color="blue") + facet_wrap( ~ term,scales="free_x") + ylab("") + xlab(labelDV)+
  geom_vline(xintercept=0) +
  geom_errorbarh(aes(xmin=condval -2*condsd,
                     xmax=condval +2*condsd), height=0) +
  theme_bw(base_size=14)


# Visualize image variance
Fcond_m1.lmm.img <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + (1 | Img), data=df_fixdur, REML = FALSE)
str(Fcond_m1.re.img <- ranef(Fcond_m1.lmm.img))
str(Fcond_df_m1.re.img <- as.data.frame(Fcond_m1.re.img))

Fcond_df_m1.re.img %>%  ggplot(aes(y=grp,x=condval)) +
  geom_point(color="blue") + facet_wrap( ~ term,scales="free_x") + ylab("") + xlab(labelDV)+
  geom_vline(xintercept=0) +
  geom_errorbarh(aes(xmin=condval -2*condsd,
                     xmax=condval +2*condsd), height=0) +
  theme_bw(base_size=14)

######### LMM FINAL #######
lmm.fixdur <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + 
                     (1 + n.Condition + n.Viewing  + n.Interaction || VP) + 
                     (1 + n.Condition + n.Viewing  + n.Interaction|| Img), 
                   data=df_fixdur, REML=FALSE, lmerControl(optimizer = "bobyqa"))

summary(lmm.fixdur)
tab_model(lmm.fixdur)
coef(lmm.fixdur)

# Plot
df_fixdur %>% 
  group_by(Condition,Viewing) %>%
  dplyr::summarize(N=n(), M=mean(DV, na.rm=T), 
                   SD = sd(DV,na.rm=T),         
                   SE= SD/sqrt(N)) %>%  
  ggplot(aes(x=Viewing, y=M, group=Condition, color=Condition)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=M-SE, ymax=M+SE), width=.1) +
  scale_color_manual(values= c("red", "green"))



# ________________________Saccade Amplitude_____________________________

# Download and filter the data
Filename1<-sprintf("SAC_mm.RData")
load(Filename1)

df_sacampl <- SAC %>%
  dplyr::filter(., nth != 1 & imtype2 != "Single(BaselineSess2)"&
                  blinkSac!= TRUE & blinkFix != TRUE)%>%
  dplyr::filter(., Ampl>0.15 & Ampl<20)%>%
  dplyr::select(.,x,y,VP,viewnr,PeakVel,Ampl,session,imtype2,Img)%>%
  mutate(DV=Ampl)%>%
  mutate(DV_log=log(DV))

df_sacampl<-df_sacampl %>% 
  mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days"
  )) %>%
  mutate(Viewing= case_when(
    viewnr=="Viewing 1" ~ "First",
    viewnr=="Viewing 2" ~ "Second"
  )) %>%
  mutate(Condition=as.factor(Condition),viewnr=as.factor(viewnr),VP=as.factor(VP),Img=as.factor(Img), DV_reci= 1/DV) %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" )) %>%
  mutate(Viewing = fct_relevel(Viewing, "First" ,"Second" )) 

# set contrasts
contrasts(df_sacampl$Viewing) <- contr.sdif(2) 
contrasts(df_sacampl$Condition) <- contr.sdif(2)

# check contrasts
contrasts(df_sacampl$Condition)
contrasts(df_sacampl$Viewing)

mat_cond <- model.matrix(~ Condition+Viewing+Condition:Viewing, df_sacampl)
df_sacampl<-data.frame(df_sacampl)
df_sacampl[, 15:18]<- mat_cond
names(df_sacampl)[15:18] <- c("Intercept", "n.Condition", "n.Viewing", "n.Interaction")

#Check if transformation is required
boxcox(DV~ n.Condition + n.Viewing  + n.Interaction, data=df_sacampl) # lambda is close to 0  => log transform is required

# Log transform dependent variable
df_sacampl <- df_sacampl %>% 
  mutate(DV=DV_log)
labelDV<-sprintf("Log Saccade Amplitude")

# Visualize subject variance
Fcond_m1.lmm.subj <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + (1 | VP), data=df_sacampl, REML = FALSE)
str(Fcond_m1.re.subj <- ranef(Fcond_m1.lmm.subj))
str(Fcond_df_m1.re.subj <- as.data.frame(Fcond_m1.re.subj))

Fcond_df_m1.re.subj %>%  ggplot(aes(y=grp,x=condval)) +
  geom_point(color="blue") + facet_wrap( ~ term,scales="free_x") + ylab("") + xlab(labelDV)+
  geom_vline(xintercept=0) +
  geom_errorbarh(aes(xmin=condval -2*condsd,
                     xmax=condval +2*condsd), height=0) +
  theme_bw(base_size=14)

# Visualize image variance
Fcond_m1.lmm.img <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + (1 | Img), data=df_sacampl, REML = FALSE)
str(Fcond_m1.re.img <- ranef(Fcond_m1.lmm.img))
str(Fcond_df_m1.re.img <- as.data.frame(Fcond_m1.re.img))

Fcond_df_m1.re.img %>%  ggplot(aes(y=grp,x=condval)) +
  geom_point(color="blue") + facet_wrap( ~ term,scales="free_x") + ylab("") + xlab(labelDV)+
  geom_vline(xintercept=0) +
  geom_errorbarh(aes(xmin=condval -2*condsd,
                     xmax=condval +2*condsd), height=0) +
  theme_bw(base_size=14)

##### LMM final #######

lmm.sacampl <- lmer(DV ~ 1 + n.Condition + n.Viewing  + n.Interaction + 
                      (1 + n.Condition + n.Viewing + n.Interaction || VP) + 
                      (1 + n.Condition + n.Viewing  + n.Interaction || Img), data=df_sacampl, 
                    REML=FALSE, lmerControl(optimizer = "bobyqa"))

summary(lmm.sacampl)
tab_model(lmm.sacampl)
coef(lmm.sacampl)


# Plot
df_sacampl %>% 
  group_by(Condition,Viewing) %>%
  dplyr::summarize(N=n(), M=mean(DV, na.rm=T), 
                   SD = sd(DV,na.rm=T),         
                   SE= SD/sqrt(N)) %>%  
  ggplot(aes(x=Viewing, y=M, group=Condition, color=Condition)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=M-SE, ymax=M+SE), width=.1) +
  scale_color_manual(values= c("red", "green"))
