set.seed(1)
setwd("E:/")
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
packageList <- c("plyr","tidyverse", "reshape2","RColorBrewer","readxl", "Hmisc")
getPackages(packageList)

p.adj.round <- function(data, dp, method="holm"){
  return(round(p.adjust(data, method=method),dp))
}

combinatorics <- function(data, N=10000, dp=3){
  movementcounts<- list()
  results<-NULL
  for (i in levels(data$Geno)){
    count <- 0
    for (j in levels(data$Treat)){
      count <- count+1
      movementcounts[[count]]<-data[data$Geno==i & data$Treat==j,]$value
    }
    water <- movementcounts[[1]][!is.na(movementcounts[[1]])]
    chitin<- movementcounts[[2]][!is.na(movementcounts[[2]])]
    mediandifference<- median(water, na.rm=T) - median(chitin, na.rm=T)
    df<-NULL
    df$data<-c(water,chitin)
    df<-as.data.frame(df)
    values<-c("water","chitin")
    permdiff<-array()
    for (k in 1:N){
      df$values[sample(1:nrow(df), nrow(df), FALSE)] <- rep(values, c(length(water), length(chitin)))
      permdiff[k]<- median(df[df$values=="water",]$data) - median(df[df$values=="chitin",]$data) 
    }
    above<-sum(abs(permdiff)>=abs(mediandifference))+1;
    perms<-N+1;
    confcombp<-binconf(above, perms, method='exact')
    results<-rbind(results,
                   cbind.data.frame(Geno=i,
                                    combp=confcombp[1],
                                    combplower=confcombp[2],
                                    combpupper=confcombp[3],
                                    mediandifference=round(mediandifference,dp)))
  }
  results[,2]<-p.adj.round(results[,2],dp)
  results[,3]<-p.adj.round(results[,3],dp)
  results[,4]<-p.adj.round(results[,4],dp)
  return(results)
}

combinatorics2 <- function(data, N=10000, dp=3){
  results<-NULL
  water<-data[data$Geno=="Col-0" & !is.na(data$value),]$value
  for (i in levels(data$Geno)[c(1,3:14)]){
    chitin<- data[data$Geno==i & !is.na(data$value),]$value
    mediandifference<- median(water, na.rm=T) - median(chitin, na.rm=T)
    df<-NULL
    df$data<-c(water,chitin)
    df<-as.data.frame(df)
    values<-c("water","chitin")
    permdiff<-array()
    for (k in 1:N){
      df$values[sample(1:nrow(df), nrow(df), FALSE)] <- rep(values, c(length(water), length(chitin)))
      permdiff[k]<- median(df[df$values=="water",]$data) - median(df[df$values=="chitin",]$data) 
    }
    above<-sum(abs(permdiff)>=abs(mediandifference))+1;
    perms<-N+1;
    confcombp<-binconf(above, perms, method='exact')
    results<-rbind(results,
                   cbind.data.frame(Geno=i,
                                    combp=confcombp[1],
                                    combplower=confcombp[2],
                                    combpupper=confcombp[3],
                                    mediandifference=round(mediandifference,dp)))
  }
  results[,2]<-p.adj.round(results[,2],dp)
  results[,3]<-p.adj.round(results[,3],dp)
  results[,4]<-p.adj.round(results[,4],dp)
  return(results)
}


bombCPK_w<-read_excel("data_Cheval2020PNAS.xlsx", sheet = 16, skip = 1)
bombCPK<-melt(bombCPK_w)%>%separate(variable,c("Geno","Treat"), sep="-")
bombRBOHD_w<-read_excel("data_Cheval2020PNAS.xlsx", sheet = 15, skip = 1)
bombRBOHD<-melt(bombRBOHD_w)%>%separate(variable,c("Geno","Treat"), sep="-")
bombLYK_w<-read_excel("data_Cheval2020PNAS.xlsx", sheet = 11, skip = 1)
bombLYK<-melt(bombLYK_w)%>%separate(variable,c("Geno","Treat"), sep="-")

bombLYK$Geno<-fct_recode(bombLYK$Geno, "Col-0"="col0")
bombLYK$Treat<-fct_relevel(bombLYK$Treat, "control","chitin")

listy<- list()
bombLYKnorm<-NULL
for (i in levels(bombLYK$Geno)){
  count <- 0
  for (j in levels(bombLYK$Treat)){
    count <- count+1
    listy[[count]]<-bombLYK[bombLYK$Geno==i & bombLYK$Treat==j,]$value
  }
  water <- listy[[1]][!is.na(listy[[1]])]
  chitin<- listy[[2]][!is.na(listy[[2]])]
  norm<-mean(water)
  water2<-water/norm
  chitin2<-chitin/norm
  bombLYKnorm<-rbind(bombLYKnorm, cbind.data.frame(Geno=i,Treat="control",value=water2),cbind.data.frame(Geno=i,Treat=j,value=chitin2))
}

ddply(bombLYKnorm, .(Geno,Treat), summarise, mean = mean(value))

bombCPK$Geno<-fct_recode(bombCPK$Geno, "Col-0"="col0")
bombCPK$Treat<-fct_relevel(bombCPK$Treat, "control","chitin")

listy<- list()
bombCPKnorm<-NULL
for (i in levels(bombCPK$Geno)){
  count <- 0
  for (j in levels(bombCPK$Treat)){
    count <- count+1
    listy[[count]]<-bombCPK[bombCPK$Geno==i & bombCPK$Treat==j,]$value
  }
  water <- listy[[1]][!is.na(listy[[1]])]
  chitin<- listy[[2]][!is.na(listy[[2]])]
  norm<-mean(water)
  water2<-water/norm
  chitin2<-chitin/norm
  bombCPKnorm<-rbind(bombCPKnorm, cbind.data.frame(Geno=i,Treat="control",value=water2),cbind.data.frame(Geno=i,Treat=j,value=chitin2))
}

ddply(bombCPKnorm, .(Geno,Treat), summarise, mean = mean(value))

bombRBOHD$Geno<-fct_recode(bombRBOHD$Geno, "Col-0"="col0")
bombRBOHD$Treat<-fct_relevel(bombRBOHD$Treat, "control","chitin")

listy<- list()
bombRBOHDnorm<-NULL
for (i in levels(bombRBOHD$Geno)){
  count <- 0
  for (j in levels(bombRBOHD$Treat)){
    count <- count+1
    listy[[count]]<-bombRBOHD[bombRBOHD$Geno==i & bombRBOHD$Treat==j,]$value
  }
  water <- listy[[1]][!is.na(listy[[1]])]
  chitin<- listy[[2]][!is.na(listy[[2]])]
  norm<-mean(water)
  water2<-water/norm
  chitin2<-chitin/norm
  bombRBOHDnorm<-rbind(bombRBOHDnorm, cbind.data.frame(Geno=i,Treat="control",value=water2),cbind.data.frame(Geno=i,Treat=j,value=chitin2))
}

ddply(bombRBOHDnorm, .(Geno,Treat), summarise, mean = mean(value))

combinatorics(bombLYKnorm)
combinatorics(bombCPKnorm)
combinatorics(bombRBOHDnorm)

allraw_control<-rbind(bombCPK[(bombCPK$Geno=="Col-0" & bombCPK$Treat=="control"),],
                      bombCPK[(bombCPK$Geno!="Col-0" & bombCPK$Treat=="control"),],
                      bombLYK[(bombLYK$Geno!="Col-0" & bombLYK$Treat=="control"),],
                      bombRBOHD[(bombRBOHD$Geno!="Col-0" & bombRBOHD$Treat=="control"),])

allraw_treated<-rbind(bombCPK[(bombCPK$Geno=="Col-0" & bombCPK$Treat=="chitin"),],
                      bombCPK[(bombCPK$Geno!="Col-0" & bombCPK$Treat=="chitin"),],
                      bombLYK[(bombLYK$Geno!="Col-0" & bombLYK$Treat=="chitin"),],
                      bombRBOHD[(bombRBOHD$Geno!="Col-0" & bombRBOHD$Treat=="chitin"),])

combinatorics2(allraw_control)
combinatorics2(allraw_treated)
