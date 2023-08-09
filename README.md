# CancerTFs

**CellNet setup file (script_setupCellNet.R)** <br>

**Set your library path so that it points to the correct platform and annotation libraries** <br>

.libPaths("~/myprog/cellnetr/packages")
library("cellnetr")

**human** <br>
library("org.Hs.eg.db")

**hgu133plus2** <br>
library("hgu133plus2.db");
library("hgu133plus2cdf");

**set up path for the CellNet objects containing the classifiers, GRNs** <br>
path_CN_obj<-"~/myprog/cellnetr/training_data/";

**change to reflect your platform** <br>
myPlatform<-"hgu133plus2"

mydir<-"~/myprog/cellnetr/final_cellnet/CellNet-master/"
source( paste(mydir, "CellNet_sourceme.R", sep='') );

**this sources all of R files** <br>
utils_sourceRs(mydir);


**CellNet main script (script_maincellnetr.R)** <br>

source("script_setupCellNet.R")

mydate<-utils_myDate();

cellnetfiledir <- "~/myprog/cellnetr"
fileName <- "final_stall_11102022"

cat("# Reading stAll file ...\n")
stQuery<-expr_readSampTab(paste(fileName,".csv",sep=""));
stQuery<-geo_fixNames(stQuery);
stAll<-utils_loadObject(paste("stQuery_",fileName,".R",sep=""));


**THIS WILL LOAD THE CEL FILES and make raw gene expression measurements** <br>

library(affy);
cat("# Reading all Cel files ...\n")
expAll<-Norm_cleanPropRaw(stAll, "hgu133plus2")

**select samples for GRN reconstruction** <br>

cat("# Creating stGRN ...\n")
stGRN<-sample_profiles_grn(stAll, minNum=84);

cat("# Creating expGRN ...\n")
expR<-expAll[,rownames(stGRN)];
expGRN<-Norm_quantNorm(expR);

#Determine latest TR annotation
hTFs<-find_tfs("Hs");

**Correlations and Zscores for GRN** <br>
corrX<-grn_corr_round(expGRN);
hTFs<-intersect(hTFs, rownames(expGRN));
zscs<-grn_zscores(corrX, hTFs);

zthresh<-6; # EMPIRICALLY DETERMINED by comparison to 3 GOLD STANDARDS.

ctGRNs<-cn_grnDoRock(stGRN, expGRN, zscs, corrX, "general", dLevelGK="description6", zThresh=zthresh);

**Make the complete CellNet object using all data** <br>
cnProc<-cn_make_processor(expAll, stAll, ctGRNs, dLevel="description1", classWeight=TRUE, exprWeight=TRUE);
fname<-paste(“cnProc_",fileName, mydate, ".R", sep='');
save(cnProc, file=fname);


**Rainbow plot (script_cellnet_rainbowPlot.R)** <br>

source("script_setupCellNet.R")

.libPaths("~/myprog/cellnetr/packages")

library("cellnetr")
library("hgu133plus2.db");
library("hgu133plus2cdf");
library("randomForest")
library(gplots)
library(ggplot2)


path_CN_obj<-"~/cellnet_2022/cancer/output_full/";
outputfileName<- outputfileName
targetCT<- targetCT
csvFile<- csvFile
cnproc_data<- cnproc_data

myPlatform<-"hgu133plus2"
cName<-"description1"

stQuery<-expr_readSampTab(paste("normal_cells_query/",csvFile,sep=""));
stQuery<-geo_fixNames(stQuery);
cnObjName<-switch(myPlatform,hgu133plus2 = paste(path_CN_obj,cnproc_data ,sep=''));
expQuery<-Norm_cleanPropRaw(stQuery, myPlatform);
cnProc<-utils_loadObject(cnObjName);
tmpAns<-cn_apply(expQuery, stQuery, cnProc, dLevelQuery=cName);
tfScores<-cn_nis_all(tmpAns, cnProc, targetCT);

pdf(paste(outputfileName,".pdf",sep=""))

cn_hmClass(tmpAns,isBig = TRUE);

**Gene regulatory network status of starting cell type (esc) GRN** <br>
cn_barplot_grnSing(tmpAns, cnProc, targetCT, c(targetCT), bOrder=NULL);

**Network influence score of HSPC GRN transcriptional regulators** <br>
cn_plotnis(tfScores[[targetCT]], limit=15);

  scoresDF<-tfScores[[targetCT]]
  xmeans<-apply(scoresDF, 2, mean);
  worst<-which.min(xmeans);
  tfs<-rownames(scoresDF)[order(scoresDF[,worst], decreasing=F)];
  scoresDF<-scoresDF[tfs,];
  topTF<-rownames(scoresDF);
  write.table(topTF,paste("TFs_",outputfileName,".txt",sep=""),row.names=F,quote=F)

plot(mp_rainbowPlot(cnProc[['expTrain']],cnProc[['stTrain']],topTF[i], dLevel="description1"))
dev.off()


**GRN_ROCS (script_cellnet_grn_rocs.sh)** <br>

source("script_setupCellNet.R")

cellnetfiledir <- "~/myprog/cellnetr"

filename="final_stall_11102022"
mydate<-utils_myDate();

stQuery<-expr_readSampTab(paste(filename,".csv",sep=""));

stQuery<-geo_fixNames(stQuery);
stAll <- utils_loadObject(paste("stQuery",filename,".R",sep=""))

library(affy);

**load expall** <br>
load(paste("expAll_",filename, mydate,".R",sep=""));
expProp <- expAll

**load stGRN** <br>
load(paste(“stGRN_",filename, mydate,".R",sep=""));

**load expGRN** <br>
load(paste("expGRN_",filename, mydate,".R",sep=""));

**load tfs** <br>
load(paste("tfs_",filename, mydate,".R",sep=""));

**load corrx htd zscs** <br>
load(paste("tmpforctGRN_",filename, mydate,".rda",sep=""));

zthresh<-6;

load(paste("ctGRNs_",filename, mydate,".R",sep=""));

**Grn report** <br>
grn_report(ctGRNs);
fname<-paste("GRN_report_",filename,".pdf", sep='');
ggsave(file=fname, width=8.5, height=11);
dev.copy(pdf, file=fname, width=8.5, height=11);
dev.off()

**select classification training and validation data** <br>
stList<-samp_for_class(stAll, prop=0.5, dLevel="description1")
lapply(stList, nrow)
stVal<-stList[['stVal']];
stTrain<-stList[['stTrain']];
expTrain<-expProp[,rownames(stTrain)];
library(randomForest)
system.time(classifiers<-cn_makeRFs(expTrain, stTrain, ctGRNs$ctGRNs$geneLists));
expVal<-expProp[,rownames(stVal)]
ansVal<-cn_classify(classifiers, expVal, ctGRNs[[3]]$geneLists)
assessed<-cn_classAssess(ansVal,stVal, classLevels="description2", resolution=0.01);

plot_class_rocs(assessed)

fname<-paste("Rocs_final_",filename,".pdf",sep='');
ggsave(file=fname, width=6, height=7);
dev.off();

