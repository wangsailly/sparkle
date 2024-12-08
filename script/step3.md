###
bioRiskPlot=function(inputFile=null, project=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
	rt=rt[order(rt$riskScore),]   
		

```R
#RISKPLOT
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
lowMax=max(rt$riskScore[riskClass=="low"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
plot(line, type="p", pch=20,
	 xlab="Patients (increasing risk socre)",
	 ylab="Risk score",
	 col=c(rep("blue",lowLength),rep("red",highLength)) )
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
dev.off()
	
#survival plot
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="blue"
pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
plot(rt$futime, pch=19,
	 xlab="Patients (increasing risk socre)",
	 ylab="Survival time (years)",
	 col=color)
legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

###
ann_colors=list()
bioCol=c("blue", "red")
names(bioCol)=c("low", "high")
ann_colors[["Risk"]]=bioCol

###
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(Risk=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
pheatmap(rt1, 
	     annotation=annotation,
	     annotation_colors = ann_colors, 
	     cluster_cols = FALSE,
	     cluster_rows = FALSE,
	     show_colnames = F,
	     scale="row",
	     color = colorRampPalette(c(rep("blue",3.5), "white", rep("red",3.5)))(50),
	     fontsize_col=3,
	     fontsize=7,
	     fontsize_row=8)
dev.off()
```
}

#tarin
bioRiskPlot(inputFile="risk.train.txt", project="train")
#test
bioRiskPlot(inputFile="risk.test.txt", project="test")
#all
bioRiskPlot(inputFile="risk.all.txt", project="all")

####################
bioSurvival=function(inputFile=null, outFile=null){
	#
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#
	diff=survdiff(Surv(futime, fustat) ~ Risk, data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
		

```R
####
surPlot=ggsurvplot(fit, 
	           data=rt,
	           conf.int=F,
	           pval=pValue,
	           pval.size=6,
	           legend.title="Risk",
	           legend.labs=c("High risk", "Low risk"),
	           xlab="Time(years)",
	           ylab="Overall survival",
	           break.time.by = 2,
	           palette=c("red", "blue"),
	           risk.table=TRUE,
	           risk.table.title="",
	           risk.table.height=.25)
###
pdf(file=outFile, width=6, height=5, onefile=FALSE)
print(surPlot)
dev.off()
```
}

###
bioSurvival(inputFile="risk.TCGAtrain.txt", outFile="surv.TCGAtrain.pdf")
bioSurvival(inputFile="risk.TCGAtest.txt", outFile="surv.TCGAtest.pdf")
bioSurvival(inputFile="risk.TCGAall.txt", outFile="surv.TCGAall.pdf")
bioSurvival(inputFile="risk.GEO.txt", outFile="surv.GEO.pdf")

####################################################
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		

```R
#
pdf(file=forestFile, width=6.5, height=4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
	
#
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
	
#
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
abline(v=1, col="black", lty=2, lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
axis(1)
dev.off()
```
}

#
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)   
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      
	

```R
#
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])

#
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 uniTab=rbind(uniTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

#multiCox
uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
```
}
indep(riskFile="risk.TCGAall.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox.txt",
      multiOutFile="multiCox.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")

####################################

#

riskFile="risk.TCGAall.txt"   
cliFile="clinical.txt"          

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#
bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")


ROC_rt=timeROC(T=risk$futime, delta=risk$fustat,
	           marker=risk$riskScore, cause=1,
	           weighting='aalen',
	           times=c(1,3,5), ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=3)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=3)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=3)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[1:3], lwd=3, bty = 'n')
dev.off()


predictTime=1     
aucText=c()
pdf(file="cliROC.pdf", width=5, height=5)

i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

for(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
legend("bottomright", aucText,lwd=3,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

####################################

riskFile="risk.TCGAall.txt"     
cliFile="clinical.txt"          

###
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

###
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

###
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "Risk")], cli)

###
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
              plots = c("density", "boxes"),    
              dencol="#FFFF00FF", boxcol="#CC99FFFF",    
              clickable=F,
              title="",              
              points=TRUE,            
              droplines=TRUE,       
              observation=rt[1,],     
              rank="sd",
              failtime = c(1,3,5),    
              prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")

###
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

###
pdf(file="calibration.pdf", width=5, height=5)
#1
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
	 xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
	   col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()