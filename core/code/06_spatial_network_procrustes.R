# Portable adaptation; statistical definitions unchanged.
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1]
pkg <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
output_path <- Sys.getenv("NAKDONG_OUTPUT", unset=file.path(pkg, "run"))
dir.create(file.path(output_path, "results"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_path, "qc"), recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(vegan);library(permute);library(igraph)})
root<-file.path(pkg,"data");base<-file.path(pkg,"data/derived");out<-output_path
meta<-read.csv(file.path(base,"optionC_metadata_preceding30d.csv"),check.names=FALSE)
rel<-read.csv(file.path(root,"bracken_relative.csv"),row.names=1,check.names=FALSE)
comm<-t(rel[,meta$sample_id])
universe<-read.csv(file.path(root,"spatial_variable_universe.csv"))
rows<-lapply(seq_len(nrow(universe)),function(i){
 v<-universe$variable[i];vals<-vapply(sort(unique(meta$site)),function(s)mean(as.numeric(meta[meta$site==s,v]),na.rm=TRUE),0.);good<-is.finite(vals)
 test<-if(sum(good)>=5&&sd(vals[good])>0)suppressWarnings(cor.test(which(good),vals[good],method="spearman",exact=FALSE))else NULL
 data.frame(variable=v,category=universe$category[i],n_sites=sum(good),rho=if(is.null(test))NA else unname(test$estimate),p=if(is.null(test))NA else test$p.value)
})
sp<-do.call(rbind,rows);sp$q<-p.adjust(sp$p,"BH");write.csv(sp,file.path(out,"results/independent_spatial_45.csv"),row.names=FALSE)
print(aggregate(as.integer(!is.na(sp$q)&sp$q<.05),list(category=sp$category),sum))
# Fig3a sequential definitions, fitted to the source's complete cases.
meta$Microcystis_RA<-colSums(rel[grep("^Microcystis ",rownames(rel)),meta$sample_id])
meta$season<-factor(meta$season,levels=c("fall","winter","spring","summer"));meta$site<-factor(meta$site)
required<-c("site","season","water_temp_C","solar_rad_avg_MJ","Chl_a_mg_m3","discharge_mean_m3s","Microcystis_RA","MC_ppb")
ok<-complete.cases(meta[,required]);m<-droplevels(meta[ok,]);dist<-vegdist(comm[ok,],"bray")
continuous<-setdiff(required,c("site","season"));m[,continuous]<-lapply(m[,continuous,drop=FALSE],function(x)as.numeric(scale(x)))
forms<-list(no_season=dist~water_temp_C+solar_rad_avg_MJ+Chl_a_mg_m3+discharge_mean_m3s+Microcystis_RA+MC_ppb,
 season_first=dist~season+water_temp_C+solar_rad_avg_MJ+Chl_a_mg_m3+discharge_mean_m3s+Microcystis_RA+MC_ppb)
seqrows<-lapply(names(forms),function(name){set.seed(42);z<-as.data.frame(adonis2(forms[[name]],data=m,by="terms",permutations=how(blocks=m$site,nperm=999)));z$term<-rownames(z);z$model<-name;z})
write.csv(do.call(rbind,seqrows),file.path(out,"results/independent_Fig3a_sequential.csv"),row.names=FALSE)
# Network reconstruction: same global top100 and documented season thresholds.
top<-order(colMeans(comm),decreasing=TRUE)[1:100];hubs<-list();networkstats<-list()
for(season in c("spring","summer","fall","winter")){
 x<-comm[meta$season==season,top,drop=FALSE];x<-x[,apply(x,2,var)>0,drop=FALSE];n<-nrow(x);rho<-cor(x,method="spearman")
 ix<-which(upper.tri(rho),arr.ind=TRUE);rv<-rho[ix];p<-2*pt(abs(rv*sqrt((n-2)/(1-rv^2))),df=n-2,lower.tail=FALSE);q<-p.adjust(p,"BH")
 sig<-abs(rv)>=if(season=="fall").55 else .6;sig<-sig&q<if(season=="fall").05 else .01
 edges<-data.frame(from=colnames(x)[ix[sig,1]],to=colnames(x)[ix[sig,2]],weight=abs(rv[sig]))
 g<-graph_from_data_frame(edges,directed=FALSE,vertices=colnames(x))
 hubs[[season]]<-data.frame(season=season,taxon=V(g)$name,degree=degree(g),betweenness_weighted=betweenness(g,normalized=TRUE),betweenness_unweighted=betweenness(g,weights=NA,normalized=TRUE))
 networkstats[[season]]<-data.frame(season=season,n=n,nodes=vcount(g),edges=ecount(g))
}
nh<-do.call(rbind,hubs);write.csv(nh,file.path(out,"results/independent_network_hubs.csv"),row.names=FALSE)
write.csv(do.call(rbind,networkstats),file.path(out,"results/independent_network_summary.csv"),row.names=FALSE)
old<-read.csv(file.path(root,"network_hubs_reference.csv"));nn<-merge(nh,old,by=c("season","taxon"),suffixes=c("_new","_old"))
print(data.frame(network_nodes_compared=nrow(nn),degree_max_error=max(abs(nn$degree_new-nn$degree_old)),betweenness_max_error=max(abs(nn$betweenness_weighted-nn$betweenness))))
# Procrustes geometry and repeated-site permutation for global raw comparison.
tax<-t(rel);tax<-tax[,order(colMeans(tax),decreasing=TRUE)[1:200]]
arg<-read.csv(file.path(root,"functions/arg_class_abundance.csv"),row.names=1)
ids<-sort(intersect(rownames(tax),rownames(arg)));ids<-ids[rowSums(arg[ids,])>0]
tax<-tax[ids,];arg<-arg[ids,];md<-meta[match(ids,meta$sample_id),]
pc<-function(v)cmdscale(vegdist(decostand(v,"hellinger"),"bray"),k=2,eig=TRUE)$points
tx<-pc(tax);ax<-pc(arg)
set.seed(42);fit<-protest(tx,ax,permutations=how(blocks=md$site,nperm=9999))
write.csv(data.frame(n=length(ids),r=unname(fit$t0),p=fit$signif),file.path(out,"results/independent_Procrustes.csv"),row.names=FALSE)

