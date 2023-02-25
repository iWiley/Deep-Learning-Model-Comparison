library(maftools)
tmp = read.table(file = "Data/RAW.XENA.TCGA.SNP.MuTect2.tsv", header=T, sep= "\t")
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
tmp = subset(tmp, substr(tmp$Tumor_Sample_Barcode, 1, 12) %in% OrginalData.TCGA$Name)
tmp$Name = substr(tmp$Tumor_Sample_Barcode, 1, 12)
groups = data.frame(
  Name = OrginalData.TCGA$Name,
  Group = ifelse(OrginalData.TCGA$prec.TIL < median(OrginalData.TCGA$prec.TIL), "Low TIL", "High TIL")
)
tmp = merge(tmp, groups, by = 'Name')
tmp$Name = NULL

tmp.h = subset(tmp, tmp$Group == "High TIL")
tmp.l = subset(tmp, tmp$Group != "High TIL")
clinc = data.frame(Tumor_Sample_Barcode = tmp$Tumor_Sample_Barcode, Group = tmp$Group)

maf.h <- read.maf(
  maf = tmp.h, 
  vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))
maf.l <- read.maf(
  maf = tmp.l, 
  vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))
maf = read.maf(
  maf = tmp, 
  vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))),
  clinicalData = clinc)

coOncoplot(
  m1=maf.h,
  m2=maf.l, 
  m1Name="High TIL",
  m2Name="Low TIL")

oncoplot(maf = maf.h)
oncoplot(maf = maf.l)

OncogenicPathways(maf = maf.h)
OncogenicPathways(maf = maf.l)

stage = clinicalEnrichment(maf=maf, clinicalFeature = 'Group')
pdf("cnv.pdf", width = 8, height = 6)
plotEnrichmentResults(stage)
dev.off()

#计算tmb值
tmb.h = tmb(maf = maf.h)
tmb.h$Group = "High TIL"
tmb.l = tmb(maf = maf.l)
tmb.l$Group = "Low TIL"

# 去除离群值
h.val = outlier_values <- boxplot.stats(tmb.h$total_perMB_log)$out 
tmb.h = subset(tmb.h, !(tmb.h$total_perMB_log %in% h.val))

l.val = outlier_values <- boxplot.stats(tmb.l$total_perMB_log)$out 
tmb.l = subset(tmb.l, !(tmb.l$total_perMB_log %in% l.val))

tmb = rbind(tmb.h, tmb.l)
wistu.website.r.tools::LoadFunc("Plot.FrameBox")
tmb$Value = tmb$total_perMB_log
Plot.FrameBox(tmb, method = "wilcox.test")
