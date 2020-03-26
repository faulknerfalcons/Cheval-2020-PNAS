generateDecayCurve<- function(f){
  result<-list()
  dat <- read.csv(f, skip=2, header=F)
  num <- dim(dat)[2]-1
  t<-dat[,1]
  for (j in 1:num){
    a<-decayCurve(t,dat[,j+1])
    value<-a[1]
    #print(a[2])
    tmp<-list(value=value,name=basename(dirname(f)))
    result[[j]]<-tmp
  }
  return(result)
}


readFrap<-function(f,preds=0){
  result<-NULL
  dat <- read.csv(f, skip=2, header=F)
  num <- dim(dat)[2]-1
  t<-dat[,1]
  preds <- 100 - predict(preds, newdata = data.frame(t = t), se=F)
  for (j in 1:num){
    y<-dat[,j+1]
    a<-analyseFrap(t,y,preds)
    ggsave(paste0(j,"-",basename(f),".png"),plotFrap(t,y,preds,a), device="png",path = "E:/FRAP/SanityCheck")
    result<-rbind(result, cbind(as.data.frame(a),name=basename(dirname(f))))
  }
  return(result)
}

analyseFrap <- function(t,y, preds){
  df<-standardiseFrap(t,y,preds)
  t<-df[,1]
  y<-df[,2]
  model<- loess(y ~ t,span=.4)
  upr<-predict(model, newdata=data.frame(t=60))
  hlf<-spline(x = model$fitted, y = model$x, xout=upr/2)$y
  return(cbind(upper=as.numeric(upr),halflife=as.numeric(hlf)))
}

decayCurve <- function(t=numeric(0),y=numeric(0), pre=5){
  maxi <- mean(head(y, pre))
  y1 <- y/maxi*100
  values<-as.data.frame(cbind(tn=1:length(t),t=t,y=y1))
  plot<-ggplot(values, aes(x=t,y=y1))+geom_smooth()+geom_point()
  return(list(values,plot))
}

standardiseFrap <- function(t=numeric(0),y=numeric(0), preds, pre=30,bleach=60, trim=TRUE){
  maxi <- mean(head(y, 5))  
  y1<-y/maxi*100
  y2<-y1+preds
  t<-t[!is.na(y2)]
  y2<-y2[!is.na(y2)]  
  maxi <- mean(head(y2, pre))  
  mini <- y2[which.min(y)]
  y3 <- y2-mini
  
  y4 <- y3/maxi*100
  
  if (trim){
    t1 <- tail(t,-which.min(y)+1)
    t2 <- t1-t1[1]
    y5 <- tail(y4,-which.min(y)+1)
    return(cbind(t=t2,y=y5))
  }else{
    return(cbind(t=t,y=y4))
  }
  
}

plotFrap <- function(t1,y1,preds,a){
  df<-standardiseFrap(t1,y1,preds)
  t<-df[,1]
  y<-df[,2]
  model<- loess(y ~ t,span=.4)
  plot<-data.frame(t=model$x,
                   yexp=model$y,
                   ymod=model$fitted
  )
  graph<-ggplot(plot, aes(x=t))+
    geom_point(aes(y=yexp),alpha=0.2)+
    geom_line(aes(y=ymod))+
    geom_hline(yintercept = a[1], colour="grey40", linetype=5)+
    geom_vline(xintercept = a[2], colour="grey40", linetype=5)+
    xlab("Time (s)")+
    ylab("Intensity (%)")+
    theme_bw()
  return(graph)
}