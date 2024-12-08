expFile="symbol.txt"           
riskFile="risk.TCGAall.txt"     
geneFile="gene.txt"           

#
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
#
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

#
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)
	
#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"Risk")]

#
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"Risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "Risk")
rt1=rt1[,sigGene]

#
rt1=melt(rt1,id.vars=c("Risk"))
colnames(rt1)=c("Risk","Gene","Expression")
	
#
group=levels(factor(rt1$Risk))
rt1$Risk=factor(rt1$Risk, levels=c("low", "high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

#
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="Risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  #outlier.shape = NA,
				  palette = c("#008B45FF","#EE0000FF") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
	
#
pdf(file="checkdiff.pdf")

#####################

#
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#
risk=read.table("risk.TCGAall.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)
	
#
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
#
boxplot=ggviolin(data, x="Risk", y="TMB", fill="Risk",
			      xlab="",
			      ylab="Tumor mutation burden (log2)",
			      legend.title="",
			      palette = c("#0066FF","#FF0000"),
			      add = "boxplot", add.params = list(fill="white"))+ 
	stat_compare_means(comparisons = my_comparisons)
	
#plot
pdf(file="riskTMB.pdf")
print(boxplot)
dev.off()

#####################################

riskFile="risk.TCGAall.txt"      
TMEfile="TMEscores.txt"          

#
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$Risk=factor(Risk$Risk, levels=c("low","high"))

#
score=read.table(TMEfile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
score=score[row.names(Risk),,drop=F]

#
rt=cbind(Risk[,"Risk",drop=F], score)

#
data=melt(rt, id.vars=c("Risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#
p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
	     xlab="",
	     ylab="TME score",
	     legend.title="Risk",
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("#0088FF", "#FF5555"), width=1)
p=p+rotate_x_text(45)
px=p+stat_compare_means(aes(group=Risk),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#
pdf(file="vioplot.pdf",)
print(px)
dev.off()

##########################################################       
immFile="CIBERSORT-Results.txt"     

#
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="low",])
treatNum=nrow(rt[rt$Risk=="high",])
data=rt
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Immune", "Expression")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Risk",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#
pdf(file="immune.diff.pdf", width=8, height=6)
dev.off()

############################################   
tideFile="TIDE.csv"              

#
tide=read.csv(tideFile, header=T, sep=",", check.names=F, row.names=1)
tide=tide[,"TIDE",drop=F]     

#
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "Risk", drop=F]
data=cbind(tide, risk)
	
#
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$Risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#
ggx=ggviolin(data, x="Risk", y="TIDE", fill = "Risk", 
	         xlab="", ylab="TIDE",
	         palette=c("#0066FF","#FF0000"),
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

#plot	
pdf(file="TIDE.pdf", width=6, height=5)
print(ggx)
dev.off()

######################################

modelFile="multiCox.txt"     
expFile="IMvigor.exp.txt"     
cliFile="IMvigor.time.txt"    

#
coef=read.table(modelFile, header=T, sep="\t", row.names=1)
row.names(coef)=gsub("`", "", row.names(coef))
coxGene=row.names(coef)

#
vigor=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
vigor=t(vigor[coxGene,])

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")

#
sameSample=intersect(row.names(vigor), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], vigor[sameSample,,drop=F])

#
trainFinalGeneExp=rt[,row.names(coef)]
actCoef=coef[,1]

#
trainFinalGeneExp=scale(trainFinalGeneExp)
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",row.names(coef))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore))
res.cut=surv_cutpoint(outTab, time = "futime", event = "fustat", variables =c("riskScore"))
cutoff=as.numeric(res.cut$cutpoint[1])
Risk=as.vector(ifelse(trainScore>cutoff, "high", "low"))
outTab=cbind(outTab, Risk)
outTab=cbind(id=rownames(outTab), outTab)
write.table(outTab, file="risk.IMvigor.txt", sep="\t", quote=F, row.names=F)

#
bioSurvival=function(inputFile=null, outFile=null){
	#
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#
	diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
		

```R
#
surPlot=ggsurvplot(fit, 
	           data=rt,
	           conf.int=F,
	           pval=pValue,
	           pval.size=6,
	           legend.title="Risk",
	           legend.labs=c("High risk", "Low risk"),
	           xlab="Time(years)",
	           break.time.by = 2,
	           palette=c("red", "blue"),
	           risk.table=F,
	       	   risk.table.title="",
	           risk.table.col = "strata",
	           risk.table.height=.25)
#
pdf(file=outFile, width=4.5, height=4, onefile = FALSE)
print(surPlot)
dev.off()
```
}

#
bioSurvival(inputFile="risk.IMvigor.txt", outFile="sur.IMvigor.pdf")

riskFile="risk.IMvigor.txt"       
cliFile="IMvigor.Response.txt"    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)

#
clinical=colnames(rt)[2]
data=rt[c("riskScore", clinical)]
colnames(data)=c("riskScore", "clinical")
data=data[(data[,"clinical"]!="unknow"),]
group=levels(factor(data$clinical))
data$clinical=factor(data$clinical, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#
boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
		          xlab="",
		          ylab="Risk score",
		          legend.title=clinical,
		          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

#
pdf(file="boxplot.pdf")