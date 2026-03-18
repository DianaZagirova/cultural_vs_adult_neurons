library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library("DGEobj.utils")
library("edgeR")
library(ggiraph)
library(tximportData)
library(GenomicFeatures)
library(tximport)
library(DESeq2)

####################### RNAseq all groups #########################
files=list.files('/media/anna/Samsung_T5/culture_vs_postmortem/RNAseq_new/cell_types',full.names = T)
sample_names=list.files('/media/anna/Samsung_T5/culture_vs_postmortem/RNAseq_new/cell_types')
ss=sapply(1:length(sample_names),function(i){unlist(strsplit(sample_names[i],'.sf')[1])})
sample_names=ss
treatment=c(rep('culture',13),rep('iPSC',3),rep('NPC',9),rep('PM',16))

dir='/media/anna/Samsung_T5/culture_vs_postmortem/quants_expr/'
tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v46.csv"),sep='\t')
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
txi$abundance
pc=read.table('/media/anna/Samsung_T5/culture_vs_postmortem/protein_coding_for_expression_data.txt',sep='\t',header = T)

samples=as.data.frame(treatment)
colnames(samples)=c('condition')
rownames(samples)=sample_names
dds <-DESeqDataSetFromTximport(txi,colData = samples,design=~condition)
rld <- rlog(dds, blind=TRUE)
rr=assay(rld)
rr=as.data.frame(rr)
colnames(rr)=sample_names


acc_rous=rowMeans(rr[,c("PM_ACC_SRR21161849_Roussos","PM_ACC_SRR21161863_Roussos","PM_ACC_SRR21161867_Roussos",
                        "PM_ACC_SRR21161872_Roussos","PM_ACC_SRR21161896_Roussos")])
riz=rowMeans(rr[,c("PM_BA9_SRR5343878_Rizzardi","PM_BA9_SRR5343882_Rizzardi","PM_BA9_SRR5343886_Rizzardi",     
                   "PM_BA9_SRR5343889_Rizzardi","PM_BA9_SRR5343892_Rizzardi","PM_BA9_SRR5343896_Rizzardi")])
dlpfc_rous=rowMeans(rr[,c("PM_DLPFC_SRR21161802_Roussos","PM_DLPFC_SRR21161822_Roussos","PM_DLPFC_SRR21161913_Roussos",   
                          "PM_DLPFC_SRR21161944_Roussos","PM_DLPFC_SRR21161979_Roussos" )])
Li_iPSCder=rr[,"CULTURE_SRR3339422_Li"]
Lu_iPSCder=rowMeans(rr[,c("CULTURE_SRR9676400_Lu","CULTURE_SRR9676401_Lu" ,"CULTURE_SRR9676402_Lu" )])
PRJNA_iPSCder=rowMeans(rr[,c("CULTURE_SRR17648092_PRJNA798046", "CULTURE_SRR17648093_PRJNA798046", "CULTURE_SRR17648094_PRJNA798046")])
Zaghi_iPSCder=rowMeans(rr[,c("CULTURE_SRR21307327_Zaghi","CULTURE_SRR21307379_Zaghi","CULTURE_SRR21307381_Zaghi")])
Ciceri50_iPSCder=rowMeans(rr[,c("Ciceri_d50_rep1","Ciceri_d50_rep2","Ciceri_d50_rep3")])
Flamier_iPSC=rowMeans(rr[,c('iPSC_SRR27991498_Flamier','iPSC_SRR27991499_Flamier','iPSC_SRR27991500_Flamier',
                            'iPSC_SRR27991501_Flamier','iPSC_SRR27991502_Flamier','iPSC_SRR27991503_Flamier',
                           'iPSC_SRR27991504_Flamier')])
Lu_iPSC=rowMeans(rr[,c('iPSC_SRR9676397_Lu','iPSC_SRR9676398_Lu','iPSC_SRR9676399_Lu')])
Topol_NPC=rowMeans(rr[,c('NPC_SRR1686362_Topol','NPC_SRR1686363_Topol','NPC_SRR1686364_Topol',
                         'NPC_SRR1686365_Topol','NPC_SRR1686366_Topol','NPC_SRR1686367_Topol')])
Lu_NPC=rowMeans(rr[,c('NPC_SRR9676403_Lu','NPC_SRR9676404_Lu','NPC_SRR9676405_Lu')])

expr=cbind(acc_rous,dlpfc_rous,riz,Li_iPSCder,Lu_iPSCder,PRJNA_iPSCder,
           Zaghi_iPSCder,Ciceri50_iPSCder,Lu_iPSC,Topol_NPC,Lu_NPC,Flamier_iPSC)

cols=c('ACC_Roussos','DLPFC_Roussos','Rizzardi','Li_iPSCder','Lu_iPSCder','PRJNA_iPSCder','Zaghi_iPSCder',
       'Ciceri50_iPSCder','Lu_iPSC','Topol_NPC','Lu_NPC','Flamier_iPSC')
colnames(expr)=cols

boxplot(expr)
write.table(expr,'/media/anna/Samsung_T5/culture_vs_postmortem/rlog_all_cortex_all_blind_true.txt',sep='\t',quote = F)

idx=intersect(pc$id,rownames(expr))
expr=expr[idx,]
write.table(expr,'/media/anna/Samsung_T5/culture_vs_postmortem/rlog_protein_coding_cortex_all_blind_true.txt',sep='\t',quote = F)

########TPM################
files=list.files('/media/anna/Samsung_T5/culture_vs_postmortem/RNAseq_new/cell_types',full.names = T)
sample_names=list.files('/media/anna/Samsung_T5/culture_vs_postmortem/RNAseq_new/cell_types')
ss=sapply(1:length(sample_names),function(i){unlist(strsplit(sample_names[i],'.sf')[1])})
sample_names=ss
treatment=c(rep('culture',13),rep('iPSC',3),rep('NPC',9),rep('PM',16))

dir='/media/anna/Samsung_T5/culture_vs_postmortem/quants_expr/'
tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v46.csv"),sep='\t')

txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreAfterBar = TRUE)
rr=txi$abundance
colnames(rr)=sample_names

acc_rous=rowMeans(rr[,c("PM_ACC_SRR21161849_Roussos","PM_ACC_SRR21161863_Roussos","PM_ACC_SRR21161867_Roussos",
                        "PM_ACC_SRR21161872_Roussos","PM_ACC_SRR21161896_Roussos")])
riz=rowMeans(rr[,c("PM_BA9_SRR5343878_Rizzardi","PM_BA9_SRR5343882_Rizzardi","PM_BA9_SRR5343886_Rizzardi",     
                   "PM_BA9_SRR5343889_Rizzardi","PM_BA9_SRR5343892_Rizzardi","PM_BA9_SRR5343896_Rizzardi")])
dlpfc_rous=rowMeans(rr[,c("PM_DLPFC_SRR21161802_Roussos","PM_DLPFC_SRR21161822_Roussos","PM_DLPFC_SRR21161913_Roussos",   
                          "PM_DLPFC_SRR21161944_Roussos","PM_DLPFC_SRR21161979_Roussos" )])
Li_iPSCder=rr[,"CULTURE_SRR3339422_Li"]
Lu_iPSCder=rowMeans(rr[,c("CULTURE_SRR9676400_Lu","CULTURE_SRR9676401_Lu" ,"CULTURE_SRR9676402_Lu" )])
PRJNA_iPSCder=rowMeans(rr[,c("CULTURE_SRR17648092_PRJNA798046", "CULTURE_SRR17648093_PRJNA798046", "CULTURE_SRR17648094_PRJNA798046")])
Zaghi_iPSCder=rowMeans(rr[,c("CULTURE_SRR21307327_Zaghi","CULTURE_SRR21307379_Zaghi","CULTURE_SRR21307381_Zaghi")])
Ciceri50_iPSCder=rowMeans(rr[,c("Ciceri_d50_rep1","Ciceri_d50_rep2","Ciceri_d50_rep3")])

expr=cbind(acc_rous,dlpfc_rous,riz,Li_iPSCder,Lu_iPSCder,PRJNA_iPSCder,
           Zaghi_iPSCder,Ciceri50_iPSCder)

cols=c('acc_rous', 'dlpfc_rous', 'riz','Li_iPSCder','Lu_iPSCder','PRJNA_iPSCder','Zaghi_iPSCder',
       'Ciceri50_iPSCder')
colnames(expr)=cols
boxplot(log(expr+1),ylim=c(-0.01,8),cex=0.3)

write.table(expr,'/media/anna/Samsung_T5/culture_vs_postmortem/all_cortex_blind_true_TPM.txt',sep='\t',quote = F)

