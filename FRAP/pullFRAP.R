## Load/install required packages
###
getPackages <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      library( i , character.only = TRUE )
    }
  }
}
packageList <- c("ggplot2","magrittr","tidyr", "dplyr")
getPackages(packageList)

## Load FRAP functions
###
setwd("E:/FRAP")
source("FRAP_functions.R")

##Generate decay curves
###
setwd("E:/FRAP/DecayCurves")
files <- list.files( full.names = T, recursive = T, pattern = ".*.csv")
output<-list()
for (f in files){
  tmp <- generateDecayCurve(f)
  output<-c(output,tmp)
}
n<-0
decaycurve<-NULL
for (i in output){
  n<-n+1
  t<-i[[1]][[1]][[2]]
  y<-i[[1]][[1]][[3]]
  tn<-i[[1]][[1]][[1]]
  name<-i[[2]][[1]]
  
  decaycurve<-rbind(decaycurve,cbind.data.frame(tname=paste(name,tn, sep="_"),tn=tn,t=t,y=y,name=name,n=n))
}
decaycurve2 <- decaycurve %>%
  separate(name, c("exp", "construct","treat"), "_") %>%
  unite(exp, construct, tn, col="tempname") %>%
  group_by(tempname) %>%
  mutate(avg_y = mean(y)) %>%
  distinct(tempname, .keep_all = T) %>%
  separate(tempname, c("exp","construct","tn"), "_")

ggplot(decaycurve2,aes(y=avg_y, x=t, color=exp, shape=construct))+geom_smooth(span=.4)+geom_point()

## Read experimental FRAP data
###
setwd("E:/FRAP/Data")
files <- list.files( full.names = T, recursive = T, pattern = ".*.csv")
output<-NULL
n<-0
for (f in files){
  n<-n+1
  bits<-unlist(strsplit(basename(dirname(f)), "_"))
  ex<-bits[1]
  construc<-bits[2]
  trea<-bits[3]
  model <- filter(decaycurve2, exp==ex, construct==construc) %>% loess(avg_y ~ t, data=.,span=.4)
  output <- readFrap(f, model)  %>% bind_rows(output, .)
}

data <- output %>%  separate(name, c("exp", "construct","treat","area"), "_")