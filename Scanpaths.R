library(tidyverse)
library(spatstat)
library(ggplot2)
library(parallel)
library(grid)
library(scales)
library(boot)
library(dplyr)
library(MASS)
library(mgcv)
library(jpeg)
library(rstudioapi)
library(gridExtra)
library(ggpubr)

setwd(dirname(getActiveDocumentContext()$path))
curr_dir <- getwd()

source("pix2deg.R")
winx <- 1920  
winy <- 1080
distCM<- 70
moWidthCM<- 53.5
px2deg <- pix2deg(1, winx, moWidthCM, distCM)
px2deg
xmargin=(winx-1500)/2
ymargin=(winy-1000)/2

xrange=c(xmargin, (winx-xmargin)) # in pixels!
yrange=c(ymargin, (winy-ymargin))

# Load files
Filename1<-sprintf("RAW_mm.RData")
load(Filename1)

df <- RAW %>%
  dplyr::mutate(y=y/px2deg,x=x/px2deg) %>% # to pixels in order to make plots with image background
  dplyr::mutate(startEvent= replace_na(startEvent, 0))%>%
  dplyr::filter(.,imtype2 != "Single(BaselineSess2)" 
                & x>xrange[1] & x<xrange[2] & y>yrange[1] & y<yrange[2]) %>%
  dplyr::mutate_at(c('trial'), as.numeric) %>%
  dplyr::select(.,x,y,VP,trial,sacDP,fixdur,session,imtype2,Img, startEvent, sacAmp)

df$sacDP[is.na(df$sacDP)] <- 0
df$Viewing <- c(0)
df$mintrial <-c(0)

df <- df %>%
  dplyr::group_by(VP,Img)%>%
  dplyr::mutate(mintrial = min(trial)) %>%
  dplyr::ungroup()

df<-df %>% 
  dplyr::mutate(Viewing = case_when(
    trial==mintrial ~ "First",
    trial!=mintrial ~ "Second"
  )) %>% 
  dplyr::mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days"
  ))

df<-df %>% 
  mutate(Condition=as.factor(Condition),Viewing=as.factor(Viewing)) %>% 
  mutate(Condition = fct_relevel(Condition, "Minutes" ,"Days" )) %>%
  mutate(Viewing = fct_relevel(Viewing, "First" ,"Second" )) %>%
  mutate(Condition=as.factor(Condition))%>%
  mutate(Viewing=as.factor(Viewing))

n <- length(df$Img)
cond <- rep(0, n)
df$cond = cond

df <- df %>%
  dplyr::mutate(cond= case_when(
    Condition=="Minutes" & Viewing=="First" ~ "First viewing - Minutes Condition",
    Condition=="Minutes" & Viewing=="Second" ~ "Second viewing - Minutes Condition",
    Condition=="Days" & Viewing=="First" ~ "First viewing - Days Condition",
    Condition=="Days" & Viewing=="Second" ~ "Second viewing - Days Condition"))

df$sacDP[is.na(df$sacDP)] <- 0

# ________________________Minutes Condition_____________________________

# ________________________1st Viewing_____________________________

# Data parameters - 1st viewing
useVP=3
useImg=103
con = 'First viewing - Minutes Condition'

# Filter data
subsetR <- df %>%
  dplyr::filter(VP==useVP, Img == useImg, cond==con)
subsetR

# Read image
directory <- sprintf("%s/%s.jpg", curr_dir, useImg)
image <- readJPEG(directory, native = TRUE)

g <- rasterGrob(image, interpolate=TRUE)
img <- subsetR$Img[1]
img <- as.numeric(img)

# Set Labels
labelVP<-sprintf('VP:%s',subsetR$VP[1])
labeltrial<-sprintf('Trial:%s',subsetR$trial[1])
labelViewing<-sprintf('Viewing:%s',subsetR$Viewing[1])
labelCondition<-sprintf('Condition:%s',subsetR$Condition[1])

#  Subset fixations
dfix1<-subsetR %>%
  dplyr::filter(.,sacDP==1 & startEvent ==1)

labelFix<-sprintf('Fixations: %d', length(dfix1$VP))
nth <- seq(1, length(dfix1$VP), by=1)
dfix1$nth <- nth

# Plot
(p1<-ggplot(subsetR,aes(x=x,y=y)) + 
    annotation_custom(grob=g,xmin=210, xmax=1710, ymin=40, ymax=1040) +
    coord_fixed(ratio=1)+
    geom_path(data=subsetR,aes(x=x,y=y, colour=as.factor(sacDP), group=1 ),size=1) + #eye movements
    geom_point(data= dfix1, aes(size = fixdur),colour='yellow', alpha =.8) + # yellow disks = fixdur
    geom_point(data= dfix1,aes(x=960, y=540), colour="white", size=4, shape=3)+ # fixation cross
    geom_text(data=dfix1,aes(label=ifelse(nth<10,as.numeric(nth),'')),hjust=.5,vjust=.5,colour='black') +#numerate fixations
    scale_colour_manual('',
                        values = c('0'='red','1'='green'),labels=c('Saccade','Fixation'))+
    xlim(0,1920)+
    ylim(0, 1080)+
    annotate('text', x = 0, y = 1070, size = 6, label = labelVP, hjust=0,colour='darkgreen')+
    annotate('text', x = 250, y = 1070, size = 6, label = labeltrial, hjust=0, colour='darkgreen')+
    annotate('text', x = 700, y = 1070, size = 6, label = labelCondition, hjust=0, colour='darkgreen')+
    annotate('text', x = 1400, y = 1070, size = 6, label = labelViewing, hjust=0, colour='darkred')+
    theme(legend.position='none')+
    xlab('X-Coordinates [px]')+ylab('Y-Coordinates [px]'))

# ___________________ 1st Viewing VS 2nd Viewing_________________

useVP2=useVP
useImg2=useImg
con2 = 'Second viewing - Minutes Condition'

# Filter data - RAW
subsetR2<- df %>%
  filter(VP==useVP2, Img==useImg2, cond==con2) 

dfix2<-subsetR2 %>%
  dplyr::filter(.,sacDP==1 & startEvent ==1)

# Labels
labelVP2<-sprintf('VP:%s',subsetR2$VP[1])
labeltrial2<-sprintf('Trial:%s',subsetR2$trial[1])
labelViewing2<-sprintf('Viewing:%s',subsetR2$Viewing[1])
labelCondition2<-sprintf('Condition:%s',subsetR2$Condition[1])

labelFix2<-sprintf('Fixations: %d', length(dfix2$VP))
nth2 <- seq(1, length(dfix2$VP), by=1)
dfix2$nth <- nth2

(p2<-ggplot(subsetR2,aes(x=x,y=y)) + 
    annotation_custom(grob=g,xmin=210, xmax=1710, ymin=40, ymax=1040) +
    coord_fixed(ratio=1)+
    geom_path(data=subsetR2,aes(x=x,y=y, colour=as.factor(sacDP), group=1 ),size=1) + #eye movements
    geom_point(data= dfix2, aes(size = fixdur),colour='yellow', alpha =.8) + # yellow disks = fixdur
    geom_point(data= dfix2,aes(x=960, y=540), colour="white", size=4, shape=3)+ # fixation cross
    geom_text(data=dfix2,aes(label=ifelse(nth2<10,as.numeric(nth2),'')),hjust=.5,vjust=.5,colour='black') +#numerate fixations
    scale_colour_manual('',
                        values = c('0'='red','1'='green'),labels=c('Saccade','Fixation'))+
    xlim(0,1920)+
    ylim(0, 1080)+
    annotate('text', x = 0, y = 1070, size = 6, label = labelVP2, hjust=0,colour='darkgreen')+
    annotate('text', x = 250, y = 1070, size = 6, label = labeltrial2, hjust=0, colour='darkgreen')+
    annotate('text', x = 700, y = 1070, size = 6, label = labelCondition2, hjust=0, colour='darkgreen')+
    annotate('text', x = 1400, y = 1070, size = 6, label = labelViewing2, hjust=0, colour='darkred')+
    theme(legend.position='none')+
    xlab('X-Coordinates [px]')+ylab('Y-Coordinates [px]'))

grid.arrange(p1, p2, nrow = 1)