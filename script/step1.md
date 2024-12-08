###

expFile="PANexp.txt" 
###
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

###
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))

###
exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
sigGene=c()
for(i in colnames(exp)[1:(ncol(exp)-1)]){
	if(sd(exp[,i])<0.001){next}
	wilcoxTest=wilcox.test(exp[,i] ~ exp[,"Type"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene2=c(sigGene, "Type")
exp=exp[,sigGene2]
diffGeneExp=t(exp[,sigGene])
diffOut=cbind(id=row.names(diffGeneExp), diffGeneExp)
write.table(diffOut, file="diffGeneExp.txt", sep="\t", row.names=F, quote=F)

###
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

###
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

###
pdf(file="boxplot.pdf", width=13, height=7)
print(p1)
dev.off()

##################################
expFile="survSigExp.txt"      #
workDir=" "     #
setwd(workDir)       #

###
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

###
maxK= 9     
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123,
              plot="png")

###
clusterNum=2      #Needs to change based on the graphical results
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$Cluster))
cluster$Cluster=letter[match(cluster$Cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F)

##################################

clusterFile="Cluster.txt"   
cliFile="time.txt"            
setwd("")      

###
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

###
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

###
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

###
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

###
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Cluster",
		           legend.labs=levels(factor(rt[,"Cluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 2,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.3)

###
print(surPlot)
dev.off()

##################################

expFile="survSigExp.txt"       
clusterFile="Cluster.txt"      
setwd("")     

###
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

###
data.pca=prcomp(data)
pcaPredict=predict(data.pca)

###
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
Cluster=as.vector(cluster[,1])

###
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
clusterCol=bioCol[1:length(levels(factor(Cluster)))]

###
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
pdf(file="PCA.pdf", width=5.5, height=4.25)
px=ggscatter(data=PCA, x="PC1", y="PC2", color="Cluster", shape="Cluster", 
	         ellipse=T, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
	         palette=clusterCol, size=2, main="PCA", legend="right")+
	         theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
print(px)
dev.off()

##################################

expFile="symbol.txt"          
gmtFile="immune.gmt"           
clusterFile="Cluster.txt"     

###
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

###
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=t(avereps(data))

###
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

###
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
###
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
###
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

###
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

####
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

###
data=melt(scoreCluster, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Immune", "Fraction")

###
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="Cluster",
     xlab="",
     ylab="Immune infiltration",
     legend.title="Cluster",
     palette=bioCol)
p=p+rotate_x_text(50)

####
pdf(file="Box.pdf")
p+stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()