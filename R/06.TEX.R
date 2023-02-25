source("00.Functions.R")
source("01.CalculateCutoff.R")

genes.TEX = "HLA-A
HLA-B
B2M
JAK1
STAT6
PDGFRB
IRF1
MET
RBM10
NFKBIE
PRKD2
EML4
LTB
LOX
WAS
PRF1
DDB2
IL7R
LATS2
CDH11
MYH11
CD79B
CBL
ZEB1
CCR7
P2RY8
IKZF3
BTK
MSI2
TCF4
ERG
PCDH17
SFMBT2
IRF4
SP140
ZNF521
AXIN2
TLL1
CR1
DCC"
genes.TEX = strsplit(genes.TEX, '\n')[[1]]

load("Data/TCGA.LIHC.Counts.RData")
load("Data/TCGA.LIHC.TPM.RData")
# TCGA.LIHC.Counts = TCGA.LIHC.TPM
OrginalData.TCGA <- read_excel("Data/OrginalData.TCGA.xlsx")
OrginalData.TCGA = subset(OrginalData.TCGA, OrginalData.TCGA$invalid == 0)
names = gsub("-", ".", OrginalData.TCGA$Name)
data.counts = TCGA.LIHC.Counts[,substr(colnames(TCGA.LIHC.Counts),14,14)  == "0"]
data.counts = data.counts[,substr(colnames(data.counts),1,12) %in% names]
data.counts = IDTran(data.counts)

# data.counts = limma::normalizeBetweenArrays(data.counts)

# save(data.counts, file = "data.counts")
# save(data.counts, file = "data.TPM")
# save(data.counts, file = "data.TPM.Nor")
# save(data.counts, file = "data.counts.Nor")

data.counts = data.frame(t(data.counts))
data.counts$Name = substr(rownames(data.counts), 1, 12)

data.Group = data.frame(
  Name = OrginalData.TCGA$Name, 
  TIL = OrginalData.TCGA$prec.TIL,
  TLS = OrginalData.TCGA$prec.TLS
)
data.Group$TIL = ifelse(data.Group$TIL < res.cut.TCGA$prec.TIL$estimate, "Low TIL", "High TIL")
data.Group$TLS = ifelse(data.Group$TLS < res.cut.TCGA$prec.TLS$estimate, "Low TLS", "High TLS")
data.Group$Name = gsub("-",".",data.Group$Name)
data.counts = merge(data.counts, data.Group, by = "Name")


wistu.website.r.tools::LoadFunc("ssGSEA")


data.counts = data.counts[!duplicated(data.counts$Name),]
rownames(data.counts) = data.counts$Name

gp = data.frame(
  row.names = rownames(data.counts),
  TLS = data.counts$TLS,
  TIL = data.counts$TIL
)
data.counts$Name = data.counts$TLS = data.counts$TIL = NULL
data.counts.ss = data.frame(t(data.counts))
result.ssGSEA = ssGSEA(list('TEX' = c(genes.TEX)), data.counts.ss)

result.ssGSEA = data.frame(t(result.ssGSEA))
result.ssGSEA = merge(result.ssGSEA, gp, by = "row.names")
rownames(result.ssGSEA) = result.ssGSEA$Row.names

# wistu.website.r.tools::LoadFunc("Plot.FrameBox")
# result.ssGSEA$Group = result.ssGSEA$TIL
# result.ssGSEA$Value = result.ssGSEA$TEX
# Plot.FrameBox(result.ssGSEA, method = "wilcox.test")
# result.ssGSEA$Group = result.ssGSEA$TLS
# Plot.FrameBox(result.ssGSEA, method = "wilcox.test")

# 
# data.counts = data.counts[,colnames(data.counts) %in% c(genes.TEX, "TIL", "TLS")]
# 
# data.counts.TIL = data.counts
# data.counts.TLS = data.counts
# 
# data.counts.TIL$Group = data.counts$TIL
# data.counts.TLS$Group = data.counts$TLS
# data.counts.TIL$Name = data.counts.TIL$TIL = data.counts.TIL$TLS = NULL
# data.counts.TLS$Name = data.counts.TLS$TIL = data.counts.TLS$TLS = NULL
# 
# rownames(data.counts.TIL) <- paste0("row_", seq(nrow(data.counts.TIL)))
# data.counts.TIL = data.counts.TIL[order(data.counts.TIL$Group),]
# 
# annotation_col = data.frame(
#   row.names = rownames(data.counts.TIL),
#   Group = data.counts.TIL$Group
# )
# colors = list(
#   Group = c("blue" , "yellow")
# )
# data.counts.TIL$Group = NULL
# data.counts.TIL = log2(data.counts.TIL + 1)
# library(pheatmap)
# library(ggplot2)
# 
# p <- pheatmap(
#   t(data.counts.TIL),
#   # scale = 'row',
#   trace = "none",
#   annotation_col = annotation_col,
#   # gaps_col = table(annotation_col)[1],
#   border = F,
#   cluster_row = T,
#   cluster_col = F,
#   show_colnames = F,
#   show_rownames = T,
#   silent = T
# )
# require(ggplotify)
# p = as.ggplot(p)
# p
# 
# 
# 
# 
# rownames(data.counts.TLS) <- paste0("row_", seq(nrow(data.counts.TLS)))
# data.counts.TLS = data.counts.TLS[order(data.counts.TLS$Group),]
# 
# annotation_col = data.frame(
#   row.names = rownames(data.counts.TLS),
#   Group = data.counts.TLS$Group
# )
# colors = list(
#   Group = c("blue" , "yellow")
# )
# 
# 
# data.counts.TLS$Group = NULL
# data.counts.TLS = log2(data.counts.TLS + 1)
# library(pheatmap)
# library(ggplot2)
# 
# p <- pheatmap(
#   t(data.counts.TLS),
#   # scale = 'row',
#   trace = "none",
#   annotation_col = annotation_col,
#   # gaps_col = table(annotation_col)[1],
#   border = F,
#   cluster_row = T,
#   cluster_col = F,
#   show_colnames = F,
#   show_rownames = T,
#   silent = T
# )
# require(ggplotify)
# p = as.ggplot(p)
# p


load("data.counts")
# load("data.TPM")
# load("data.TPM.Nor")
# load("data.counts.Nor")

data.counts = data.frame(t(data.counts))
dt = data.frame(
  row.names = rownames(data.counts),
  `PD-1` = data.counts$PDCD1,
  `PD-L1` = data.counts$CD274,
  CTLA4 = data.counts$CTLA4
)
dt$Name = substr(rownames(dt),1,12)
dt = dt[!duplicated(dt$Name),]
rownames(dt) = dt$Name
dt = merge(dt, gp,by= "row.names")


mergeData = function(data, mergeCol){
  data = data.frame(data)
  mergeCol = gsub("-", ".", mergeCol)
  rows = c()
  for (i in 1:length(rownames(data))) {
    row = data[i,]
    mcol = data[i, mergeCol]
    for (d in mergeCol) {
      row = row[,-which(colnames(row) == d)]
    }
    for (d in mergeCol) {
      rows = rbind(rows, c(unlist(row), d[[1]], data[i,d][[1]]))
    }
  }
  data = data.frame(rows)
  return(data)
}

data = mergeData(dt, mergeCol = c("TIL", "TLS"))
data = mergeData(data, mergeCol = c("PD-1", "PD-L1", "CTLA4"))

data$Row.names = NULL
data$V5 = NULL
colnames(data) = c("Group", "X", "Y")
data$Y = as.numeric(data$Y)
data$Y = log2(data$Y + 1)
data$X = gsub("\\.", "-", data$X)

data$X = as.factor(data$X)
data$Group = as.factor(data$Group)

data.TIL = subset(data, data$Group == "Low TIL" | data$Group == "High TIL")
data.TLS = subset(data, data$Group != "Low TIL" & data$Group != "High TIL")

ggplot(data = data.TIL, 
       aes(x = X, y = Y, fill = Group)) +
  stat_boxplot(width = .5)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  xlab(NULL)+
  ylab("Expression level")

ggplot(data = data.TLS, 
       aes(x = X, y = Y, fill = Group)) +
  stat_boxplot(width = .5)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  xlab(NULL)+
  ylab("Expression level")


dt$Group = dt$TIL
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')

dt = data.frame(
  row.names = rownames(data.counts),
  Value = data.counts$PDCD1
)
dt = merge(dt, gp,by= "row.names")
dt$Group = dt$TLS
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')

dt = data.frame(
  row.names = rownames(data.counts),
  Value = data.counts$CD274
)
dt = merge(dt, gp,by= "row.names")
dt$Group = dt$TIL
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')

dt = data.frame(
  row.names = rownames(data.counts),
  Value = data.counts$CD274
)
dt = merge(dt, gp,by= "row.names")
dt$Group = dt$TLS
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')


dt = data.frame(
  row.names = rownames(data.counts),
  Value = data.counts$CTLA4
)
dt = merge(dt, gp,by= "row.names")
dt$Group = dt$TIL
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')



dt = data.frame(
  row.names = rownames(data.counts),
  Value = data.counts$CTLA4
)
dt = merge(dt, gp,by= "row.names")
dt$Group = dt$TLS
dt$Value = log2(dt$Value + 1)
Plot.FrameBox(dt, method = 'wilcox.test')



