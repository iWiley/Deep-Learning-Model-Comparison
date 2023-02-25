warning("Please manually run the code line by line and save the image from the preview box.")
source("00.Functions.R")
# We need to get the cutoff value
source("01.CalculateCutoff.R")

OrginalData.TCGA <- read_excel("Data/OrginalData.TCGA.xlsx")
OrginalData.TCGA = subset(OrginalData.TCGA, OrginalData.TCGA$invalid == 0)

data.estimation = read.csv("Data/data.estimation.csv", header = T)
row.names(data.estimation) = data.estimation$X
data.estimation$X = NULL
data.estimation$Name = substr(rownames(data.estimation), 1 ,12)

data.estimation.IMM = data.frame(
  Name = gsub("-", ".", OrginalData.TCGA$Name), 
  TIL = OrginalData.TCGA$prec.TIL,
  TLS = OrginalData.TCGA$prec.TLS
)
data.estimation.IMM = merge(data.estimation.IMM, data.estimation, by = "Name")
data.estimation.IMM = data.estimation.IMM[!duplicated(rownames(data.estimation.IMM)),]

data.estimation.IMM$TIL = ifelse(data.estimation.IMM$TIL < res.cut.TCGA$prec.TIL$estimate, "Low TIL", "High TIL")
data.estimation.IMM$TLS = ifelse(data.estimation.IMM$TLS < res.cut.TCGA$prec.TLS$estimate, "Low TLS", "High TLS")

data.Enr.dif = c()
data.estimation.IMM.TILh = subset(data.estimation.IMM, data.estimation.IMM$TIL == "High TIL")
data.estimation.IMM.TILl = subset(data.estimation.IMM, data.estimation.IMM$TIL != "High TIL")
data.estimation.IMM.TLSh = subset(data.estimation.IMM, data.estimation.IMM$TLS == "High TLS")
data.estimation.IMM.TLSl = subset(data.estimation.IMM, data.estimation.IMM$TLS != "High TLS")
for (cn in 1:length(colnames(data.estimation.IMM))) {
  if (cn < 4) {
    next
  }
  ptil = wilcox.test(data.estimation.IMM.TILh[,cn], data.estimation.IMM.TILl[,cn])
  ptls = wilcox.test(data.estimation.IMM.TLSh[,cn], data.estimation.IMM.TLSl[,cn])
  data.Enr.dif = rbind(
    data.Enr.dif,
    c(
      row.names = colnames(data.estimation.IMM.TILh)[cn],
      TIL = ptil$p.value,
      TLS = ptls$p.value
    )
  )
}
data.Enr.dif = data.frame(data.Enr.dif)
data.Enr.dif = subset(data.Enr.dif, data.Enr.dif$TIL < .05)

data.estimation.IMM = data.estimation.IMM[!duplicated(data.estimation.IMM$Name),]
rownames(data.estimation.IMM) = data.estimation.IMM$Name
data.estimation.IMM = data.estimation.IMM[, colnames(data.estimation.IMM) %in% c(data.Enr.dif$row.names, "TLS", "TIL")]

data.estimation.IMM$T.cell.CD4...non.regulatory._QUANTISEQ = NULL

library(tidyr)
dt1 <- data.estimation.IMM %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample") %>% 
  pivot_longer(cols=4:26,
               names_to= "celltype",
               values_to = "Proportion") %>%   
  pivot_longer(cols=2:3,
               names_to= "Group")
dt1$Proportion = log2(dt1$Proportion+1)
dt1$Proportion = log2(dt1$Proportion+1)

dt1$value = factor(dt1$value, levels = c("High TIL", "Low TIL","High TLS", "Low TLS"))

pdf(file = "estimation.pdf", width = 10, height = 6)
ggplot(
    data = dt1,
    aes(x = celltype,y = Proportion, fill=value),
    xlab = "",
     palette = "jco",ylab = "Cell composition"
  ) +   
    # geom_violin(trim = T, scale = "width") +
    geom_boxplot(position = position_dodge(0.9), width = 0.4) +
    stat_compare_means(
      aes(group = value),
      label = "p.signif",
      method = "wilcox.test",
      hide.ns = T,
      size = 4.5 
    ) + 
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ))
dev.off()



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
result.est$TIL = ifelse(result.est$prec.TIL > median(result.est$prec.TIL), "High TIL", "Low TIL")
result.est$TLS = ifelse(result.est$prec.TLS > median(result.est$prec.TLS), "High TLS", "Low TLS")


Plot.FrameBox = function (data, title = "", method = "wilcox.test") 
{
    CheckPackage(c("ggplot2", "ggpubr", "ggthemes"))
    ggplot(data, aes(x = Group, y = Value, color = Group)) + 
        stat_boxplot(geom = "errorbar", width = 0.15) + geom_boxplot() + 
        geom_jitter(shape = 16, position = position_jitter(0.2)) + 
        labs(title = title, x = "", y = "") +
        stat_compare_means(method = method, lable = "p.format") + 
        ylim(-4,4)+
        scale_color_tableau()
}

# TIL
result.est$Group = result.est$TIL
result.est$Value = scale(result.est$StromalScore)
pdf(file = "est.TIL.StromalScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est, method = "wilcox.test")
dev.off()

result.est$Value = scale(result.est$ImmuneScore)
pdf(file = "est.TIL.ImmuneScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est, method = "wilcox.test")
dev.off()

result.est$Value = scale(result.est$ESTIMATEScore)
pdf(file = "est.TIL.ESTIMATEScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

result.est$Value = scale(result.est$TumorPurity)
pdf(file = "est.TIL.TumorPurity.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

# TLS
result.est$Group = result.est$TLS
result.est$Value = scale(result.est$StromalScore)
pdf(file = "est.TLS.StromalScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

result.est$Value = scale(result.est$ImmuneScore)
pdf(file = "est.TLS.ImmuneScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

result.est$Value = scale(result.est$ESTIMATEScore)
pdf(file = "est.TLS.ESTIMATEScore.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

result.est$Value = scale(result.est$TumorPurity)
pdf(file = "est.TLS.TumorPurity.pdf", width = 3, height = 6)
Plot.FrameBox(result.est)
dev.off()

# result.est$prec.TIL = result.est$prec.TLS = NULL
# result.est$Name = NULL

# result.est = pivot_longer(result.est, 
# cols=1:4,
# names_to= "Item",
# values_to = "Value")

# result.est = pivot_longer(result.est, 
# cols=1:2,
# names_to= "Group")


# result.est$Value = log2(result.est$Value + 1)
# pdf(file = "est.pdf", width = 6, height = 6)
# ggplot(
#     data = result.est,
#     aes(x = Item, y = Value, fill=value),
#     xlab = "",
#      palette = "jco",ylab = "Cell composition"
#   ) +   
#     # geom_violin(trim = T, scale = "width") +
#     geom_boxplot(position = position_dodge(0.9), width = 0.8) +
#     stat_compare_means(
#       aes(group = value),
#       label = "p.signif",
#       method = "wilcox.test",
#       hide.ns = T,
#       size = 4.5 
#     )
# dev.off()

