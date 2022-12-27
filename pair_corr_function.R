library(ggplot2)
library(grid)
library(scales)
library(boot)
library(dplyr)
library(tidyverse)
library(MASS)
library(sjPlot)
library(data.table)
library(spatstat)
library(rstudioapi)
library(afex)
library(Rmisc) 
library(jpeg) 
library(parallel)
library(stargazer)
library(gridExtra)
library(pROC)
library(MESS)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

winx <- 1920  
winy <- 1080
distCM<- 70
moWidthCM<- 53.5
source("pix2deg.R")
px2deg <- pix2deg(1, winx, moWidthCM, distCM)

xmargin=(winx-1500)/2
ymargin=(winy-1000)/2

xrange=c(xmargin, (winx-xmargin))
yrange=c(ymargin, (winy-ymargin))

# Load files
Filename1<-sprintf("RAW_mm.RData")
load(Filename1)

RAW$trial = as.numeric(RAW$trial)

# Filter data
d <- RAW %>%
  dplyr::mutate(y=y/px2deg,x=x/px2deg) %>% 
  dplyr::filter(., startEvent==1 & sacDP==1 &
                  x>xrange[1] & x<xrange[2] &
                  y>yrange[1] & y<yrange[2]) %>%
  dplyr::filter(., imtype2 != "Single(BaselineSess2)") %>%
  dplyr::select(.,VP,viewnr,imtype2,Img,x,y)%>%
  dplyr::mutate(Condition= case_when(
    imtype2=="rep:Minutes" ~ "Minutes",
    imtype2=="rep:Days" ~ "Days")) %>%
  dplyr::mutate(Viewing= case_when(
    viewnr==1 ~ "First",
    viewnr==2 ~ "Second")) %>%
  dplyr::mutate(Condition=as.factor(Condition))%>%
  dplyr::mutate(Viewing=as.factor(Viewing))%>%
  dplyr::mutate(VP=as.factor(VP))%>%
  dplyr::mutate(Img=as.factor(Img))

n <- length(d$Img)
cond <- rep(0, n)
d$cond = cond

d <- d %>%
  dplyr::mutate(cond= case_when(
    Condition=="Minutes" & Viewing=="First" ~ "First viewing - Minutes Condition",
    Condition=="Minutes" & Viewing=="Second" ~ "Second viewing - Minutes Condition",
    Condition=="Days" & Viewing=="First" ~ "First viewing - Days Condition",
    Condition=="Days" & Viewing=="Second" ~ "Second viewing - Days Condition"))

d <- d %>%
  dplyr::select(VP,Img,x,y,Viewing,cond) 

d <- d[, c(1, 2, 5, 6, 3, 4)]

names(d)[names(d) == "Img"] <- "image"
names(d)[names(d) == "VP"] <- "id"


####### Step 1 - Simulation of inhomogeneous and homogeneous control processes ####### 

# compute sigma
sigma <- NULL
densityImages <- NULL
for (img in unique(d$image)){
  for (cond in unique(d$cond)){
    idx <- d$cond==cond & d$image==img
    dImg <- d[idx,]
    ppImg <- ppp(dImg$x,dImg$y,window=owin(xrange,yrange))
    bwImg <- mean(bw.scott(ppImg))
    sigma <- rbind(sigma,data.frame(img,cond,bwImg))
    denImg <- density(ppImg,sigma=bwImg)
    densityImages <- rbind(densityImages,data.frame(img,cond,denImg$v))
    sim <- rpoint(nrow(dImg),denImg)
    d$xposInhom[idx] <- round(sim$x,digits=3)
    d$yposInhom[idx] <- round(sim$y,digits=3)
  }}

# add inhomogeneous and homogeneous simulations in dataframe
sim <- rpoint(nrow(d))
d$xposHom <- round(sim$x*diff(xrange)+xrange[1],digits=3)
d$yposHom <- round(sim$y*diff(yrange)+yrange[1],digits=3)

dx <- gather(d,'pp','x',c(5,7,9))
dy <- gather(d,'pp','y',c(6,8,10))
d <- cbind(dx[,-c(5,6,7)],y=dy$y)
d$pp <- factor(d$pp,levels=c('x','xposInhom','xposHom'),
               labels=c('Experiment','Inhomogeneous','Homogeneous'))
table(d$pp)

# Plot fixations for inhomogeneous and homogeneous point processes
# Set image and VP
im = 60
VP = 2
con = 'First viewing - Days Condition'
con2 = 'Second viewing - Days Condition'

x <- readJPEG("../60.jpg", native = TRUE)
g <- rasterGrob(x, interpolate=TRUE)
labelCondition<-sprintf('Condition:%s',dImg$cond[1])

# subsetting data for one image
dImg1 <- filter(d,image==im & cond==con)
bwimg1 <- filter(sigma,img==im & cond==con)
bw1 <- bwimg1$bwImg

# Fixation density plot
(p1 <- ggplot(data=dImg1,aes(x=x,y=y,col=pp)) +
    stat_density2d(aes(alpha=..density..),show.legend=FALSE,
                   geom="raster",contour=FALSE,h=bw1*c(1,1),fill="black") +
    geom_point(size=.01) +
    xlim(0,1920)+
    ylim(0, 1080)+
    facet_grid(.~pp) +
    labs(x='x-Coordinate [pix]',y='y-Coordinate [pix]',colour='Point Process') +
    coord_fixed(xlim=xrange,ylim=yrange,expand = FALSE) +
    theme_bw(base_size = 16) +
    theme(legend.position="none"))

# Heatmap of fixations
(p2<-ggplot(data=dImg1,aes(x=x,y=y)) +
    annotation_custom(grob=g,xmin=0, xmax=1920, ymin=0, ymax=1080) +
    stat_density2d(data=dImg1, aes(x=x, y=y,
                                   fill = ..level.., alpha = ..level..), size= 10, bins= 50, geom='polygon') + 
    theme_bw() +scale_fill_gradient(low = "blue", high = "red") +
    scale_alpha_continuous(range=c(0.01,0.25) , guide = FALSE) +
    xlim(0,1920)+
    ylim(0, 1080))


# Plot scanpath from experiment and the corresponding simulated processes
dTrial <- filter(d,image==im & cond==con & id==VP)

(p3<-ggplot(data=dTrial,aes(x=x,y=y,col=pp)) +
    stat_density2d(data=dImg,aes(alpha=..density..),show.legend=FALSE,
                   geom="raster",contour=FALSE,fill="black") +
    geom_point(size=2) +  geom_path(size=.5) +
    labs(x='x-Coordinate [pix]',y='y-Coordinate [pix]',colour='Point Proces') +
    coord_fixed(xlim=xrange,ylim=yrange,expand=FALSE) +
    theme_bw(base_size=16) +
    facet_grid(.~pp,drop=TRUE) +
    theme(legend.position="none"))

####### Step 2 - Choose optimal bandwidth for intensity estimation of PCF ####### 

# Get the data back to degrees
xrange=c(xmargin*px2deg, (winx-xmargin)*px2deg)
yrange=c(ymargin*px2deg, (winy-ymargin)*px2deg)

d <- d %>%
  mutate(y=y*px2deg,x=x*px2deg)

# average PCF
mPCF <- function(pcf){
  r <- pcf$r[1,]
  trans <- sapply(data.frame(pcf$pcf),mean)
  mPCF <- data.frame(r,trans,row.names=NULL)
  return(mPCF)
}

# PCF deviation function
devPcfComp <- function(g,rmin,rmax){
  indx <- g$r>=rmin & g$r<=rmax
  dev <- sum((g$trans[indx]-1)^2)*(g$r[2]-g$r[1])
  return(dev)
}

# compute PCF for fixation sequence
fixSeqPCF <- function(fixSeq,den,p){
  
  # start, end, number of trials
  idx <- which(diff(as.numeric(fixSeq$id))!=0)
  i1 <- c(1,idx+1)
  i2 <- c(idx,nrow(fixSeq))
  Nid <- length(i1)
  numA <-  nrow(fixSeq)
  
  alvr <- NULL
  alvt <- NULL
  alvs <- NULL
  hdev <- NULL
  aid  <- NULL
  VPs <- NULL
  
  for (t in 1:Nid ) {
    scanpath <- fixSeq[i1[t]:i2[t],]
    
    if ( nrow(scanpath)>p$maxlen )  {
      scanpath <- scanpath[(1:p$maxlen),]
    } # maxlen
    
    if ( nrow(scanpath)>=p$minlen ) {
      numX <- nrow(scanpath)
      id <- unique(scanpath$id)
      X <- ppp(scanpath$x,scanpath$y,window=owin(xrange,yrange))
      g <- pcfinhom(X,den,kernel="epanechnikov",r=p$reval,divisor="r") # divisor="d" | "r"
      g$trans <-  g$trans*numA/numX
      ginteg <- devPcfComp(g,p$rmin,p$rmax)
      hdev <- c(hdev,ginteg)
      alvt <- rbind(alvt,g$trans)
      alvs <- rbind(alvs,rep(t,length(g$trans)))
      VPs <- rbind(VPs,rep(id,length(g$trans)))
      alvr <- rbind(alvr,g$r)
      aid <- c(aid,id)
    } # minlen
  } # id
  pcfSummary <- list(r=alvr,pcf=alvt,trial=alvs,dev=hdev,minlen=p$minlen,
                     maxlen=p$maxlen,rmin=p$rmin,rmax=p$rmax,nId=Nid,
                     excluded=Nid-length(hdev),id=aid,VP=VPs)
  return(pcfSummary)
} # function

optimalBandwidthPcf <- function(bw,img,con){
  pExp <- filter(d,pp=='Experiment' & image==img & cond==con)
  ppExp <- ppp(x=pExp$x,y=pExp$y,window=owin(xrange,yrange))
  denExp <- density(ppExp,sigma=bw,edge=FALSE)
  denExp$v[denExp$v<epsilon] <- epsilon
  pInhom <- filter(d,pp=='Inhomogeneous' & image==img & cond==con)
  
  pcfInhom <- fixSeqPCF(pInhom,denExp,p)
  meanPcfInhom <- mPCF(pcfInhom)
  dev <- devPcfComp(meanPcfInhom,p$rmin,p$rmax)
  type <- unique(pExp$type)
  
  return(data.frame(bw,img,con,dev))
}

epsilon <- .Machine$double.eps # set smallest possible computable number

p <- NULL
p$minlen <- 2     # minimum length of a scanpath
p$maxlen <- 100    # maximum length of a scanpath, otherwise truncated
p$rmin <- .1       # evaluation window of PCF (minimum)
p$rmax <- 6.5      # evaluation window of PCF (maximum)
p$sigmaMin <- .1   # minimum bandwidth for optimal bandwidth evaluation  (Step 3)
p$sigmaMax <- 10   # maximum bandwidth for optimal bandwidth evaluation (Step 3)
p$sigmaSteps <- .1 # stepsize for optimal bandwidth evaluation (Step 3)
p$reval <- seq(0,7,length.out=513) # compute PCF from to

# evaluate PCFs for all bandwidths on all images and each presentation separately
eval_bw <- seq(p$sigmaMin,p$sigmaMax,by=p$sigmaSteps)
eval_img <- unique(d$image)
eval_cond <- unique(d$cond)
parOptBw <- cbind(data.frame(bw=rep(eval_bw,each=length(eval_img)),img=rep(eval_img,times=length(eval_bw))),con=rep(eval_cond,each=length(eval_bw)*length(eval_img)))

cl <- makeForkCluster(detectCores()-1)
deviation <- parLapply(cl, seq_len(nrow(parOptBw)), function(i) {
  do.call(optimalBandwidthPcf, as.list(parOptBw[i,]))
})
stopCluster(cl)

# deviation of the inhomogeneous point process for each bandwidth, image, and presentation
devInhom <- do.call(rbind,lapply(deviation,as.data.frame))

# save result
save(d,devInhom,p,xrange,yrange,file='devInhom.Rdata')

# load data
Filename1<-sprintf("devInhom.Rdata")
load(Filename1)

# plot bandwidths
ggplot(devInhom,aes(x=bw,y=dev,group=img)) +
  geom_line() +
  facet_grid(.~con) +
  labs(x='Bandwidth Sigma',y='Deviation') +
  coord_cartesian(ylim=c(0,5)) +
  theme_bw(base_size=16)


####### Step 3 - Compute pair correlation function ####### 

computePcfImage <- function(img,con,bw){
  pExp <- filter(d,pp=='Experiment' & image==img & cond==con)
  ppExp <- ppp(x=pExp$x,y=pExp$y,window=owin(xrange,yrange))
  denExp <- density(ppExp,sigma=bw,edge=FALSE)
  
  pInhom <- filter(d,pp=='Inhomogeneous' & image==img & cond==con)
  ppInhom <- ppp(x=pInhom$x,y=pInhom$y,window=owin(xrange,yrange))
  denInhom <- density(ppInhom,sigma=bw,edge=FALSE)
  
  pHom <- filter(d,pp=='Homogeneous' & image==img & cond==con)
  ppHom <- ppp(x=pHom$x,y=pHom$y,window=owin(xrange,yrange))
  denHom <- density(ppHom,sigma=bw,edge=FALSE)
  
  pcfExp <- fixSeqPCF(pExp,denExp,p)
  pcfInhom <- fixSeqPCF(pInhom,denInhom,p)
  pcfHom <- fixSeqPCF(pHom,denHom,p)
  
  pcfImage <- rbind(data.frame(pp='Experiment',id=c(pcfExp$VP),
                               r=c(pcfExp$r),pcf=c(pcfExp$pcf)),
                    data.frame(pp='Inhomogeneous',id=c(pcfInhom$VP),
                               r=c(pcfInhom$r),pcf=c(pcfInhom$pcf)),
                    data.frame(pp='Homogeneous',id=c(pcfHom$VP),
                               r=c(pcfHom$r),pcf=c(pcfHom$pcf)))
  
  type <- unique(pExp$type)
  
  return(data.frame(img,con,pcfImage))
  return(g)
}

which(is.na(devInhom$bw))

devInhom[is.na(devInhom)] <- 5

parPcf <- devInhom %>% 
  dplyr::group_by(.,img,con) %>%
  dplyr::summarize(.,bw=bw[which.min(dev)])

cl <- makeForkCluster(detectCores()-1)
pcfDataSet <- parLapply(cl, seq_len(nrow(parPcf)), function(i) {
  do.call(computePcfImage, as.list(parPcf[i,]))
})
stopCluster(cl)

pcfDataSet <- do.call(rbind,lapply(pcfDataSet,as.data.frame))

pcfDataSet <- pcfDataSet %>%
  dplyr::filter(.,r>0.1)%>%
  dplyr::mutate(con= case_when(
    con=="First viewing - Minutes Condition" ~ "1 viewing - Minutes",
    con=="Second viewing - Minutes Condition" ~ "2 viewing - Minutes",
    con=="First viewing - Days Condition" ~ "1 viewing - Days",
    con=="Second viewing - Days Condition" ~ "2 viewing - Days"))

pcfTrial <- pcfDataSet %>%
  dplyr::filter(.,img==103 & con=='1 viewing - Minutes') %>%
  dplyr::group_by(.,pp,id,r) %>%
  dplyr::summarize(.,pcf=mean(pcf))

meanPcf <- pcfDataSet %>%
  dplyr::group_by(.,pp,con,r) %>%
  dplyr::summarize(.,pcf = mean(pcf, na.rm=TRUE))

meanPcfImg <- pcfDataSet %>%
  dplyr::group_by(.,pp,con,img,r) %>%
  dplyr::summarize(.,pcf = mean(pcf, na.rm=TRUE))

# Plots
# PCF for each trial
ggplot(pcfTrial,aes(x=r,y=pcf,col=pp)) +
  geom_line(aes(group=id),col='gray70') +
  geom_line(data=filter(pcfTrial,id%in%c(2,3,4)), aes(lty=as.factor(id)),col='black') +
  geom_line(stat='summary',fun.y='mean') +
  facet_grid(.~pp) +
  coord_cartesian(ylim=c(0,3)) +
  labs(x='Distance r [°]',y='Pair Correlation Function g(r)',colour='Point Process',lty='Participant') +
  theme_bw(base_size=16)

# PCF for each image
ggplot(meanPcfImg,aes(x=r,y=pcf,col=pp)) +
  geom_line(aes(group=img),col='gray70') +
  geom_line(data=meanPcf) +
  facet_grid(con~pp) +
  coord_cartesian(ylim=c(0,3)) +
  labs(x='Distance r [°]',y='Pair Correlation Function g(r)',colour='Point Process') 

# PCF for 4 conditions - only experiment
meanPCFs <- meanPcfDataSet %>%
  dplyr::filter(.,pp=='Experiment') %>%
  dplyr::group_by(.,con,r) %>%
  dplyr::summarize(.,pcf = mean(pcf, na.rm=TRUE),SE=sd(pcf))

pPcfDataSet <- ggplot(meanPCFs,aes(x=r,y=pcf,col=con)) +
  geom_line() +
  coord_cartesian(ylim=c(0,6), xlim=c(0,6)) +
  labs(x='Distance r [°]',y='Pair Correlation Function g(r)',colour='Condition') +
  theme_bw(base_size=16)+
  theme(legend.position=c(.8,.75))
print(pPcfDataSet)


####### Compute area under the curve #######

meanPCFs <- meanPCFs %>%
  dplyr::filter(.,r>1)

##### Minutes - First
minutes_first <- meanPCFs %>%
  dplyr::filter(.,con=='1 viewing - Minutes')
int_min_1 <- discreteIntegrationSimpsonRule(minutes_first$pcf)

##### Minutes - Second
minutes_second <- meanPCFs %>%
  dplyr::filter(.,con=='2 viewing - Minutes')
int_min_2 <- discreteIntegrationSimpsonRule(minutes_second$pcf)

##### Days - First
days_first <- meanPCFs %>%
  dplyr::filter(.,con=='1 viewing - Days')
int_days_1 <- discreteIntegrationSimpsonRule(days_first$pcf)

##### Days - Second
days_second <- meanPCFs %>%
  dplyr::filter(.,con=='2 viewing - Days')
int_days_2 <- discreteIntegrationSimpsonRule(days_second$pcf)

auc_all <- setNames(data.frame(matrix(ncol = 3, nrow = 4)), c("Viewing", "Condition", "AUC"))
auc_all$Viewing <- c('First', 'Second', 'First', 'Second')
auc_all$Condition <-c('Minutes', 'Minutes', 'Days', 'Days')
auc_all$AUC <-c(int_min_1,int_min_2,int_days_1,int_days_2)


