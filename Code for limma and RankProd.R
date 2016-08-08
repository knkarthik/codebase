		##limma starts here##

library(limma)

targets = readTargets(file = "targets",sep = "\t")

#rawObj has raw intensities; we do not normally directly work on them
rawObj = read.maimages(files = targets, source = "agilent",green.only = TRUE)

#background corrected signals are stored in corrObj
corrObj <- backgroundCorrect(rawObj, method="normexp", offset=16)

#normObj has the normalized values
normObj <- normalizeBetweenArrays(corrObj, method="quantile")

#normObj.ave is same as the 'y.ave' we're familiar of.
normObj.ave <- avereps(normObj, ID=normObj$genes$ProbeName)

#Optionally at this point, we can write the normalized values to a file.
#For RankProd we directly work on normObj.ave to save time! 
write.table(normObj.ave, file="normalized.txt", sep="\t", quote=FALSE)

## This part is only for DE analysis using limma ##

f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(normObj.ave, design)
contrast.matrix <- makeContrasts("Control-SRNS", "Control-SSNS", "SRNS-SSNS", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
output.SRNS <- topTable(fit2, adjust="BH", coef="Control-SRNS", genelist=normObj.ave$genes, number=Inf)
output.SSNS <- topTable(fit2, adjust="BH", coef="Control-SSNS", genelist=normObj.ave$genes, number=Inf)
output.SRNS_vs_SSNS <- topTable(fit2, adjust="BH", coef="SRNS-SSNS", genelist=normObj.ave$genes, number=Inf)
write.table(output.SRNS, file="Control_vs_SRNS.txt", sep="\t", quote=FALSE)
write.table(output.SSNS, file="Control_vs_SSNS.txt", sep="\t", quote=FALSE)
write.table(output.SRNS_vs_SSNS, file="SRNS_vs_SSNS.txt", sep="\t", quote=FALSE)


                ####RankProd part starts here###

library(RankProd)
##create data for rankprod
#Column 1,2 of normObj.ave$E are controls, 3,4 are SRNS and 5,6,7,8 are SSNS
#Row 1, 2 of normObj.ave$E are positional probes, so we take everything from row 3 to the last #row
ctrl_vs_SRNS = log2((normObj.ave$E[3:nrow(normObj.ave$E),1:4]))
ctrl_vs_SSNS = log2((normObj.ave$E[3:nrow(normObj.ave$E),c(1,2,5:8)]))
SRNS_vs_SSNS = log2((normObj.ave$E[3:nrow(normObj.ave$E),c(3,4,5:8)]))

## RankProd for ctrl vs SRNS ##
cl.ctrl_vs_SRNS <- rep(c(0,1),c(2,2))
origin.ctrl_vs_SRNS <- rep(1, 4)
RP.out.SRNS <- RP(ctrl_vs_SRNS,cl.ctrl_vs_SRNS, num.perm=100, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)
top_RP.out.SRNS = topGene(RP.out.SRNS,cutoff=0.05,method="pval",logged=TRUE,logbase=2,gene.names=gene.names)
#Here DOWN_Ctrl_vs_SRNS means down-regulated in control (i.e up-regulated in SRNS).
write.table(top_RP.out.SRNS$Table1, file="DOWN_ctrl_vs_SRNS", sep="\t", quote=FALSE)
write.table(top_RP.out.SRNS$Table2, file="UP_ctrl_vs_SRNS", sep="\t", quote=FALSE)


## Now, for ctrl vs SSNS ##

cl.ctrl_vs_SSNS <- rep(c(0,1),c(2,4))
origin.ctrl_vs_SSNS <- rep(1, 6)
RP.out.SSNS <- RP(ctrl_vs_SSNS,cl.ctrl_vs_SSNS, num.perm=100, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)
top_RP.out.SSNS = topGene(RP.out.SSNS,cutoff=0.05,method="pval",logged=TRUE,logbase=2,gene.names=gene.names)
#Here DOWN_Ctrl_vs_SRNS means down-regulated in control (i.e up-regulated in SSNS).
write.table(top_RP.out.SSNS$Table1, file="DOWN_ctrl_vs_SSNS", sep="\t", quote=FALSE)
write.table(top_RP.out.SSNS$Table2, file="UP_ctrl_vs_SSNS", sep="\t", quote=FALSE)

## Now, for SRNS vs SSNS ##

cl.SRNS_vs_SSNS <- rep(c(0,1),c(2,4))
origin.SRNS_vs_SSNS <- rep(1, 6)
RP.out.SRNS_vs_SSNS <- RP(SRNS_vs_SSNS,cl.SRNS_vs_SSNS, num.perm=100, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)
top_RP.out.SRNS_vs_SSNS = topGene(RP.out.SRNS_vs_SSNS,cutoff=0.05,method="pval",logged=TRUE,logbase=2,gene.names=gene.names)
#Here DOWN_Ctrl_vs_SRNS means down-regulated in SRNS (i.e up-regulated in SSNS).
write.table(top_RP.out.SRNS_vs_SSNS$Table1, file="DOWN_SRNS_vs_SSNS", sep="\t", quote=FALSE)
write.table(top_RP.out.SRNS_vs_SSNS$Table2, file="UP_SRNS_vs_SSNS", sep="\t", quote=FALSE)


		##heatmaps with dendrogram (using limma data)##

#First, get the top 100 DE genes (number = 100) for each comparison.
sel_1 = topTable(fit2, adjust="BH", coef="Control-SRNS", genelist=normObj.ave$genes, number=100)
#optionally, we can export this as a file by:
#write.table(normObj.ave$E[rownames(sel_1),], file="one.txt", sep="\t", quote = FALSE)
#Now, get the averaged expression values for these genes from normObj.ave$E object 
SRNS = normObj.ave$E[rownames(sel_1),]
#Repeat the above steps for all comparisons
sel_2 = topTable(fit2, adjust="BH", coef="Control-SSNS", genelist=normObj.ave$genes, number=100)
#write.table(normObj.ave$E[rownames(sel_2),], file="two.txt", sep="\t", quote = FALSE)
SSNS = normObj.ave$E[rownames(sel_2),]
sel_3 = topTable(fit2, adjust="BH", coef="SRNS-SSNS", genelist=normObj.ave$genes, number=100)
#write.table(normObj.ave$E[rownames(sel_3),], file="three.txt", sep="\t", quote = FALSE)
SRNS-vs-SSNS = normObj.ave$E[rownames(sel_3),]

#combine the three lists into one
combined = rbind(SRNS, SSNS, SRNS-vs-SSNS)

#generate heatmap, adding col=topo.colors(100) will give Yellow-Green-Blue
heatmap(combined)

#heatmap.2 from gplots will give the color-key
install.packages("gplots")
library("gplots")
heatmap.2(exprs(esetSel), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

##Create dendrograms separately from heatmap

#compute pair-wise distance between the columns
xdist <- dist(t(combined), diag=TRUE, upper=TRUE, method='manhattan')
round(xdist,2)
#image(as.matrix(xdist))
#u implement a hierarchical clustering (agglomerative)
hc <- hclust(xdist, method='average')
plot(as.dendrogram(hc), horiz=FALSE)

