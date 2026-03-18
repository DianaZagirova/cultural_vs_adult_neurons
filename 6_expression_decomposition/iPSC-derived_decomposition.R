library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)

############MuSiC#######################
sc_eset = readRDS('/tank/projects/v_sagitova/regulons/Jeffries2025/pfc.clean.rds')
sc_eset <- subset(sc_eset, subset = group=='adult' & (batch==210430 | batch==210209 | batch== 200128) & nCount_RNA > 500 & nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cm3=read.csv("/tank/projects/hic_decomposition/RNA_decompose/music_input/liter_MkKenzie_PA_add3types_enhanced_union_7_types.txt",sep='\t')
cm3=cm3$gene

bulk=read.table('/tank/projects/hic_decomposition/RNA_decompose/music_input/cultures_samples_salmon_counts_gene_name.txt', header=T)


bulk=bulk[!duplicated(bulk$gene.name),]
rownames(bulk)=bulk$gene.name
bulk=bulk[,colnames(bulk)[2:ncol(bulk)]]

sample_info <- data.frame(
  condition = rep("X",dim(bulk)[2]), 
  row.names = colnames(bulk)
)

pdata_annotated <- new("AnnotatedDataFrame", data = sample_info)
eset_bulk <- ExpressionSet(assayData = as.matrix(bulk),phenoData = pdata_annotated)
bulk.mtx = exprs(eset_bulk)

genes_from_sc=rownames(sc_eset[['RNA']]$data)
featrs=intersect(rownames(bulk.mtx),cm3)
featrs=intersect(featrs,genes_from_sc)

sobj <- subset(sc_eset, features=featrs)
counts_mat <- LayerData(sobj, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(list(counts=as.matrix(counts_mat)),
                            colData=DataFrame(pdata_df))

bulk.mtx=bulk.mtx[featrs,]

music_prop = music_prop(bulk.mtx = bulk.mtx, sc.sce = sce, clusters = 'new_clusters3',
                        samples = 'orig.ident', verbose = F)

write.table(music_prop$Est.prop.weighted,'/tank/projects/hic_decomposition/RNA_decompose/music_input/cultures_batch_210430_210209_200128_union_liter_MkKenzie_PA.txt',quote=F,sep='\t')

samp <- RunUMAP(sc_eset, features = featrs)
DimPlot(samp, reduction = "umap",group.by = "new_clusters3")


rr=read.table('/media/anna/Samsung_T5/SCZ_new/rnaseq_counts/music_input/cultures_batch_210430_210209_200128_union_liter_MkKenzie_PA.txt')
Li_iPSCder=colMeans(rr["CULTURE_SRR3339422_Li",])
Lu_iPSCder=colMeans(rr[c("CULTURE_SRR9676400_Lu","CULTURE_SRR9676401_Lu" ,"CULTURE_SRR9676402_Lu" ),])
Ballarino_iPSCder=colMeans(rr[c("CULTURE_SRR17648092_PRJNA798046", "CULTURE_SRR17648093_PRJNA798046", "CULTURE_SRR17648094_PRJNA798046"),])
Zaghi_iPSCder=colMeans(rr[c("CULTURE_SRR21307327_Zaghi","CULTURE_SRR21307379_Zaghi","CULTURE_SRR21307381_Zaghi"),])

barplot(Zaghi_iPSCder[c("AST","Endothelial","Excitotary_Neurons","Inhibitory_Neuron","Microglia","Oligodendrocytes","OPC" )],
        names.arg=c('Astro','Endo','Ex','Inh','Mcg','Olig','OPCs'),ylab = 'percent',ylim=c(0,0.8),main='Music,\n iPSC-derived neurons',cex.main=0.8)


################CIBERSORT###################
library(CIBERSORT)
types=read.table('/media/anna/Samsung_T5/SCZ_new/rnaseq_counts/input_cibersort/CA_for_decompose.txt')
types=types[,c("Astrocytes","Endothelia","Excitatory","Inhibitory","Microglia",       
               "Oligodendrocytes","OPCs" )]
dat=read.table('/media/anna/Samsung_T5/SCZ_new/rnaseq_counts/cultures_samples_salmon_abundance_ens_id.txt',sep='\t')

genes=intersect(rownames(types),rownames(dat))
dat=as.matrix(dat[genes,])
types_sel=as.matrix(types[genes,])
results <- cibersort(types_sel, dat)

write.table(results,'/media/anna/Samsung_T5/SCZ_new/decompose/cultures.txt',quote=F)

rr=read.table('/media/anna/Samsung_T5/SCZ_new/decompose/cultures.txt')

Li_iPSCder=colMeans(rr["CULTURE_SRR3339422_Li",])
Lu_iPSCder=colMeans(rr[c("CULTURE_SRR9676400_Lu","CULTURE_SRR9676401_Lu" ,"CULTURE_SRR9676402_Lu" ),])
PRJNA_iPSCder=colMeans(rr[c("CULTURE_SRR17648092_PRJNA798046", "CULTURE_SRR17648093_PRJNA798046", "CULTURE_SRR17648094_PRJNA798046"),])
Zaghi_iPSCder=colMeans(rr[c("CULTURE_SRR21307327_Zaghi","CULTURE_SRR21307379_Zaghi","CULTURE_SRR21307381_Zaghi"),])

barplot(Li_iPSCder[c("Astrocytes","Endothelia","Excitatory","Inhibitory","Microglia","Oligodendrocytes","OPCs" )],
        names.arg=c('Astro','Endo','Ex','Inh','Mcg','Olig','OPCs'),ylab = 'percent',ylim=c(0,0.8),main='CIBERSORT,\n iPSC-derived neurons',cex.main=0.8)









