expFile="symbol.txt"        
cluFile="Cluster.txt"      
logFCfilter=0.585           
adj.P.Val.Filter=0.05       

###

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=log2(data+0.1)

###
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=t(avereps(data))

###
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

###
sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]

###
geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)

###
allDiffGenes=c()
for(i in 1:ncol(comp)){
	fit=lmFit(data, design)
	contrast=paste0(comp[2,i], "-", comp[1,i])
	#print(contrast)
	cont.matrix=makeContrasts(contrast, levels=design)
	fit2=contrasts.fit(fit, cont.matrix)
	fit2=eBayes(fit2)
	

```R
###
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)

###
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
geneList[[contrast]]=row.names(diffSig)
```
}

###
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

###
interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)

##############################

pvalueFilter=0.05    
qvalueFilter=0.05     

###
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}


rt=read.table("interGene.txt", header=F, sep="\t", check.names=F)     

###
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

###GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

###
pdf(file="barp1.pdf")
bar1=barp(kk, drop=TRUE, showCategory=10, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar1)
dev.off()
		
#KEGG
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

###
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

###
pdf(file="barp2.pdf")
barp2(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

##############################

coxPfilter=0.05       
tcgaFile="TCGA.expTime.txt"    
geoFile="GEO.expTime.txt"      

###
rt=read.table(tcgaFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365
	#############Group the data#############
	inTrain=createDataPartition(y=rt[,2], p=0.5, list=F)
	train=rt[inTrain,]
	test=rt[-inTrain,]
	trainOut=cbind(id=row.names(train),train)
	testOut=cbind(id=row.names(test),test)
	

```R
#train group
outUniTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(train[,3:ncol(train)])){
	#cox
	cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]

	#####
	if(coxP<coxPfilter){
	    sigGenes=c(sigGenes,i)
		outUniTab=rbind(outUniTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			         )
uniSigExp=train[,sigGenes]
uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
if(length(sigGenes)<5){next}

#lasso
x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
if(nrow(geneCoef)<2){next}

#############COX model#############
multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

####
outMultiTab=data.frame()
outMultiTab=cbind(
	               coef=multiCoxSum$coefficients[,"coef"],
	               HR=multiCoxSum$conf.int[,"exp(coef)"],
	               HR.95L=multiCoxSum$conf.int[,"lower .95"],
	               HR.95H=multiCoxSum$conf.int[,"upper .95"],
	               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)

#train group
riskScore=predict(multiCox,type="risk",newdata=train)         
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,Risk=risk))
	
#test group
riskScoreTest=predict(multiCox,type="risk",newdata=test)      
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,Risk=riskTest))

#GEO
GEO=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
GEO$futime=GEO$futime/365
geoScore=predict(multiCox, type="risk", newdata=GEO)
geoRisk=as.vector(ifelse(geoScore>medianTrainRisk, "high", "low"))
GEO=cbind(GEO[,outCol], riskScore=as.vector(geoScore), Risk=geoRisk)
geoRiskOut=cbind(id=rownames(GEO), GEO)

#High and low pvalue
diff=survdiff(Surv(futime, fustat) ~Risk,data = trainRiskOut)
pValue=1-pchisq(diff$chisq, df=1)
diffTest=survdiff(Surv(futime, fustat) ~Risk,data = testRiskOut)
pValueTest=1-pchisq(diffTest$chisq, df=1)
diffGEO=survdiff(Surv(futime, fustat) ~Risk, data=GEO)
pValueGEO=1-pchisq(diffGEO$chisq, df=1)

#ROC
predictTime=1    
roc=timeROC(T=train$futime, delta=train$fustat,
            marker=riskScore, cause=1,
            times=c(predictTime), ROC=TRUE)
rocTest=timeROC(T=test$futime, delta=test$fustat,
            marker=riskScoreTest, cause=1,
            times=c(predictTime), ROC=TRUE)	

if((pValue<0.03) & (roc$AUC[2]>0.6) & (pValueTest<0.049) & (rocTest$AUC[2]>0.6) & (pValueGEO<0.049)){
	###
	write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
	write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
	###
	write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
	write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    #lasso
    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
	pdf("lasso.lambda.pdf")
	plot(fit, xvar = "lambda", label = TRUE)
	dev.off()
	pdf("lasso.cvfit.pdf")
	plot(cvfit)
	abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
	dev.off()
    ####
    outMultiTab=outMultiTab[,1:2]
	write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
	write.table(trainRiskOut,file="risk.TCGAtrain.txt",sep="\t",quote=F,row.names=F)
	write.table(testRiskOut,file="risk.TCGAtest.txt",sep="\t",quote=F,row.names=F)
	write.table(geoRiskOut,file="risk.GEO.txt",sep="\t",quote=F,row.names=F)
	#####
	allRiskOut=rbind(trainRiskOut, testRiskOut)
	write.table(allRiskOut,file="risk.TCGAall.txt",sep="\t",quote=F,row.names=F)
	break
}
```
}

#############################

cluFile="Cluster.txt"            
riskFile="risk.TCGAall.txt"       

###
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

###
sameSample=intersect(row.names(cluster), row.names(risk))
data=cbind(risk[sameSample,,drop=F], cluster[sameSample,,drop=F])
###
data$Cluster=factor(data$Cluster, levels=levels(factor(data$Cluster)))
group=levels(factor(data$Cluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

###
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cluster)))]
	
###
boxplot=ggboxplot(data, x="Cluster", y="riskScore", color="Cluster",
			      xlab="Cluster",
			      ylab="Risk score",
			      legend.title="Cluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)

#plot
pdf(file="clusterRisk.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()

###
rt=data[,c("Cluster", "Risk", "fustat")]
colnames(rt)=c("Cluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

####
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + 
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()