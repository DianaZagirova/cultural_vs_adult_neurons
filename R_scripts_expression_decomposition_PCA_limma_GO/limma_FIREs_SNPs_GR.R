BiocManager::install("limma")

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)

tab=read.table('/media/anna/Samsung_T5/culture_vs_postmortem/data_for_limma_from_unified_fires_add_samples.txt',sep='\t',header = T)
colnames(tab)
tab
group=as.factor(c(rep('neu',6),rep('cult',7)))
restr=as.factor(c('dpn','dpn','hind','ar','ar','ar','dpn','dpn','dpn','ar','ar','ar','dpn'))
samples=as.data.frame(group)
samples$restrict=restr
colnames(samples)=c('group','restrict')
sample_names=c("Plet","Rahman","Hu","Heffel_i","Heffel_a","Tian","Wu", "Li","Rajarajan","Ballerino","Zaghi","Rahman_iPSC","Lagar")
rownames(samples)=sample_names
samples
design=model.matrix(~group+restrict,samples)
design
fit <- lmFit(tab, design)
fit2 <- eBayes(fit, trend=F)
summary(decideTests(fit2))
f=topTable(fit2, coef='groupneu', number=Inf)

dim(f[(f$adj.P.Val<0.05)&(f$logFC<0),])
neu_spec=f[(f$adj.P.Val<0.05)&(f$logFC>0),]
cult_spec=f[(f$adj.P.Val<0.05)&(f$logFC<0),]
ns=tab[rownames(neu_spec),]
cs=tab[rownames(cult_spec),]

boxplot(ns)
boxplot(cs)
final_tab=f
rownames(final_tab)=as.integer(rownames(final_tab))-1

write.table(final_tab,'/media/anna/Samsung_T5/culture_vs_postmortem/limma_result_with_additional_samples_fires_with_Heffel_infant.txt',sep='\t',quote=F)


######################################################################

aa=read.table('/media/anna/Samsung_T5/culture_vs_postmortem/cult_spec_fires_genes_with_Hef_inf.txt',header = F)
aa
bb=read.table('/media/anna/Samsung_T5/culture_vs_postmortem/genes_for_FIREs_background_tpm_more1_in_iPSC_derived.txt',header=T)

ego <- enrichGO(gene          = aa$V1,
                OrgDb         = org.Hs.eg.db,
                universe = bb$name,
                keyType       = 'ENSEMBL',
                ont           = 'BP',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

ego
ego=simplify(ego, by = "qvalue",select_fun = min, cutoff = 0.7,measure="Lin",semData=dmGO)

ego@result$geneID
ego@result$Description[ego@result$Count>5]
dotplot(ego, showCategory=10) 

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = gos, mart = mart)
ego
#########################################################################
tab=read.table('/media/anna/Samsung_T5/culture_vs_postmortem/snp_recalling_for_limma_2_cell_types_BD_from_tab28_no_rahman_aggregated.txt',sep='\t',header = T,row.names = 1)
rownames(tab)
cols=colnames(tab)
dim(tab)
cols
colnames(tab)

group=as.factor(c(rep('cult',6),rep('neu',6)))
restr=as.factor(c('ar','dpn','dpn','ar','dpn','dpn','ar','ar','hind','dpn','dpn','ar'))
samples=as.data.frame(group)
samples$restrict=restr
colnames(samples)=c('group','restrict')

sample_names=c("Ballarino","Li","Wu","Zaghi","Rajarajan","Lagarkova",           
               "Heffel_infant", "Heffel_adult","Hu","Pletenev","Rahman","Tian" )
rownames(samples)=sample_names
samples
design=model.matrix(~group+restrict,samples)
design
fit <- lmFit(tab, design)
fit2 <- eBayes(fit, trend=F)
summary(decideTests(fit2))
f=topTable(fit2, coef='groupneu', number=Inf)
f
dim(f[(f$adj.P.Val<0.05),])
f[(f$adj.P.Val<0.05),]
neu_spec=f[(f$adj.P.Val<0.05)&(f$logFC>0),]
cult_spec=f[(f$adj.P.Val<0.05)&(f$logFC<0),]
dim(neu_spec)
dim(cult_spec)
ns=tab[rownames(neu_spec),]
cs=tab[rownames(cult_spec),]
f[(f$P.Value<0.05)&(f$logFC>0),]
#dev.off()
boxplot(c(ns,cs),ylim=c(0,5))
boxplot(cs)
final_tab=f
#rownames(final_tab)=as.integer(rownames(final_tab))-1
rownames(tab)
final_tab
write.table(final_tab,'/media/anna/Samsung_T5/culture_vs_postmortem/limma_result_recalling_20K_BD_tab28_2_cell_types.txt',sep='\t',quote=F)

