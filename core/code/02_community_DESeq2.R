# Portable adaptation; statistical definitions unchanged.
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1]
pkg <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
output_path <- Sys.getenv("NAKDONG_OUTPUT", unset=file.path(pkg, "run"))
dir.create(file.path(output_path, "results"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_path, "qc"), recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(vegan);library(permute);library(DESeq2)})
root <- file.path(pkg,"data")
base <- file.path(pkg,"data/derived")
out <- output_path
meta <- read.csv(file.path(base,"optionC_metadata_preceding30d.csv"),check.names=FALSE)
rel <- read.csv(file.path(root,"bracken_relative.csv"),row.names=1,check.names=FALSE)
comm <- t(rel[,meta$sample_id])
stopifnot(identical(rownames(comm),meta$sample_id))
distance <- vegdist(comm,"bray")
rows <- list()
for(term in c("month","season","site")){
 perm <- how(blocks=factor(if(term=="site") meta$month else meta$site),within=Within(type="free"),nperm=999)
 set.seed(42)
 model <- adonis2(distance~group,data=data.frame(group=factor(meta[[term]])),permutations=perm)
 rows[[term]] <- data.frame(term=term,R2=model$R2[1],F=model$F[1],p=model[1,"Pr(>F)"])
}
write.csv(do.call(rbind,rows),file.path(out,"results/independent_permanova.csv"),row.names=FALSE)
print(do.call(rbind,rows))
countpath <- file.path(root,"bracken_counts.csv")
if(!file.exists(countpath))stop("Missing raw count authority: ",countpath)
counts <- read.csv(countpath,row.names=1,check.names=FALSE)
counts <- counts[,meta$sample_id]
keep <- rowSums(counts>0)>=5 & rowSums(counts)>=10
counts <- round(as.matrix(counts[keep,]));storage.mode(counts)<-"integer"
meta$site<-factor(meta$site);meta$season<-factor(meta$season,levels=c("fall","winter","spring","summer"))
meta$water_temp_C <- as.numeric(scale(meta$water_temp_C))
meta$bloom_status <- factor(ifelse(meta$MC_ppb>1,"Bloom","Non-bloom"),levels=c("Non-bloom","Bloom"))
rownames(meta)<-meta$sample_id
design<-~site+season+water_temp_C+bloom_status
mm<-model.matrix(design,meta);stopifnot(qr(mm)$rank==ncol(mm))
dds<-DESeqDataSetFromMatrix(counts,meta,design)
dds<-DESeq(dds,quiet=TRUE)
z<-as.data.frame(results(dds,contrast=c("bloom_status","Bloom","Non-bloom")))
z$taxon<-rownames(z)
z$display<-!is.na(z$padj)&z$padj<0.05&abs(z$log2FoldChange)>=1&!is.na(z$lfcSE)&z$lfcSE<1.5
old<-read.csv(file.path(base,"v10_site_aware_sensitivity/fig3b_deseq2_site_blocked_all.csv"))
m<-merge(z,old,by="taxon",suffixes=c("_new","_old"))
summary<-data.frame(features=nrow(z),display=sum(z$display),positive=sum(z$display&z$log2FoldChange>0),negative=sum(z$display&z$log2FoldChange<0),max_lfc_error=max(abs(m$log2FoldChange_new-m$log2FoldChange_old),na.rm=TRUE),max_q_error=max(abs(m$padj_new-m$padj_old),na.rm=TRUE),na_q_pattern_mismatch=sum(is.na(m$padj_new)!=is.na(m$padj_old)))
write.csv(z,file.path(out,"results/independent_deseq2.csv"),row.names=FALSE)
write.csv(summary,file.path(out,"results/independent_deseq2_validation.csv"),row.names=FALSE)
print(summary)
writeLines(capture.output(sessionInfo()),file.path(out,"qc/R_sessionInfo.txt"))
