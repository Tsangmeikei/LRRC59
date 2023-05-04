load('TCGA_XXX.RData')
tcga.data=log(dt+1)
load('Clinical.Rdata')
comsamples=intersect(c(paste0(rownames(clinical),'-01'),paste0(rownames(clinical),'-03'),paste0(rownames(clinical),'-06')),colnames(tcga.data))
gs=c('XX','XX')

clini.cut=clinical[substr(comsamples,1,12),]
times=clini.cut$OS.time/365
status=clini.cut$OS

data.fin=data.frame(time=times,status=status,t(tcga.data[gs,comsamples]))

plotKMCox=function(dat,genes,type,mypal){
  colnames(dat)=c('time','status','groups')
  dat=dat[which(dat[,1]!='NA'&dat[,2]!='NA'&dat[,3]!='NA'),]
  gp=c('Low exp','High exp')
  vls=1:length(gp)
  gvls=vls[match(dat[,3],gp)]
  dt=data.frame(data.frame(dat[,1],dat[,2],gvls))
  aa=coxFun(dt)
  fit<-survfit(Surv(time,status) ~ groups,data=dat)
  # cox=coxph(Surv(time,status) ~ groups,data=dat)
  # b=summary(cox)
  hr=round(as.numeric(aa[2]),3)
  lower.hr=round(as.numeric(aa[3]),3)
  upper.hr=round(as.numeric(aa[4]),3)
  p=round(as.numeric(aa[1]),3)
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  plot(sf, mark.time = TRUE,col=mypal,xlab=paste("survival years (",type," )"),ylab = "survival rate",main=genes,lwd=2,cex.axis=1.3,cex.lab=1.5,font=2)
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),' ( N = ',sdf$n,')'), col = mypal,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE,cex = 1.2)
  text(x=max(fit$time)/2,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,3),"(","95%CI,",signif(lower.hr,3),"-",signif(upper.hr,3),")"),
       bty="n",font=2)
  return(p)
}
timeType=='OS'
mypal = pal_nejm(alpha = 0.7)(7)
for (aaa in gs) {
  label=ifelse(as.numeric(data.fin[,aaa])>median(as.numeric(data.fin[,aaa])),'High exp','Low exp')
  aa=plotKMCox(data.frame(data.fin$time,data.fin$status,label),aaa,timeType,mypal)
}