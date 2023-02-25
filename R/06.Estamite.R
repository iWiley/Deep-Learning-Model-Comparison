library(estimate)
source("00.Functions.R")
# We need to get the cutoff value
source("01.CalculateCutoff.R")

load("Data/TCGA.LIHC.TPM.RData")
library(readxl)
OrginalData.TCGA <- read_excel("Data/OrginalData.TCGA.xlsx")
OrginalData.TCGA = subset(OrginalData.TCGA, OrginalData.TCGA$invalid == 0)
names = gsub("-", ".", OrginalData.TCGA$Name)
TCGA.LIHC.Counts = TCGA.LIHC.TPM
data.counts = TCGA.LIHC.Counts[,substr(colnames(TCGA.LIHC.Counts),14,14)  == "0"]
data.counts = data.counts[,substr(colnames(data.counts),1,12) %in% names]
data.counts = log2(data.counts + 1)
wistu.website.r.tools::LoadFunc("IDTrans")
data.counts = IDTrans(data.counts)

write.table(data.counts, file = "data.counts.txt", sep = '\t', quote = F)
outputGCT("data.counts.txt", 'data.counts.gct')   
filterCommonGenes("data.counts.txt", 'data.counts.gct', id="GeneSymbol")

estimateScore("data.counts.gct", "results.data.counts.gct")

plotPurity(scores="results.data.counts.gct", samples="TCGA.CC.A5UC.01A.11R.A28V.07")

result.est = read.table('../results.data.counts.gct', skip = 2, header = T)
rownames(result.est) = result.est$NAME
result.est$NAME = result.est$Description = NULL
name = rownames(result.est)
result.est = sapply(result.est, as.numeric)
rownames(result.est) = name

write.csv(result.est, file = "result.est.csv")

result.est = data.frame(t(result.est))
result.est$Name = substr(rownames(result.est), 1, 12)

result.est = result.est[!duplicated(result.est$Name),]

til = data.frame(
    Name = gsub("-", ".", OrginalData.TCGA$Name),
    prec.TIL = OrginalData.TCGA$prec.TIL,
    prec.TLS = OrginalData.TCGA$prec.TLS
)

result.est = merge(result.est, til , by = 'Name')
result.est = result.est[!duplicated(result.est$Name),]
result.est$Group = ifelse(result.est$prec.TIL > median(result.est$prec.TIL), "High TIL", "Low TIL")

wistu.website.r.tools::LoadFunc("Plot.FrameBox")

result.est$Value = result.est$StromalScore
Plot.FrameBox(result.est, method = "wilcox.test")
# .022
result.est$Value = result.est$ImmuneScore
Plot.FrameBox(result.est, method = "wilcox.test")
# .025
result.est$Value = result.est$ESTIMATEScore
Plot.FrameBox(result.est, method = "wilcox.test")
# .025
result.est$Value = result.est$TumorPurity
Plot.FrameBox(result.est, method = "wilcox.test")