proxy = 'http://127.0.0.1:7890'
Sys.setenv("http_proxy" = proxy)
Sys.setenv("https_proxy" = proxy)
warning("Please manually run the code line by line and save the image from the preview box.")
source("00.Functions.R")
# We need to get the cutoff value
source("01.CalculateCutoff.R")

load("Data/TCGA.LIHC.Counts.RData")
OrginalData.TCGA <- read_excel("Data/OrginalData.TCGA.xlsx")
OrginalData.TCGA = subset(OrginalData.TCGA, OrginalData.TCGA$invalid == 0)
names = gsub("-", ".", OrginalData.TCGA$Name)
data.counts = TCGA.LIHC.Counts[,substr(colnames(TCGA.LIHC.Counts),14,14)  == "0"]
data.counts = data.counts[,substr(colnames(data.counts),1,12) %in% names]
data.counts = IDTran(data.counts)
data.counts = data.frame(t(data.counts))
data.counts$Name = substr(rownames(data.counts), 1, 12)
  
data.Group = data.frame(
  Name = OrginalData.TCGA$Name, 
  TIL = OrginalData.TCGA$prec.TIL
)
data.Group$TIL = ifelse(data.Group$TIL < res.cut.TCGA$prec.TIL$estimate, "Low TIL", "High TIL")
data.Group$Name = gsub("-",".",data.Group$Name)

data.counts = merge(data.counts, data.Group, by = "Name")
data.counts$Group = data.counts$TIL
data.counts$Name = data.counts$TIL = NULL

result.deg = edgeR(data.counts)
result.deg = subset(result.deg, result.deg$sig != "none")

# GO(rownames(result.deg), type = "BP")
# GO(rownames(result.deg), type = "CC")
# GO(rownames(result.deg), type = "MF")
# GO(rownames(result.deg), type = "DO")

wistu.website.r.tools::LoadFunc("Enrichment")
Enrichment = function(deg, genes = NULL, geneType = "auto", pCutoff = 0.05, organism = 'hsa', OrgDb="org.Hs.eg.db"){
  if (!is.null(genes)) {
    deg = subset(deg, rownames(deg) %in% genes)
  }
  logFC = data.frame(
    Genes = rownames(deg),
    logFC = deg$logFC
  )
  if (geneType == "auto") {
    if (substr(rownames(deg)[[1]], 1, 4) == "ENSG") {
      geneType = "ensg"
    } else {
      geneType = "symbol"
    }
  }
  geneType = ifelse(geneType == "ensg", "ENSEMBL", "SYMBOL")
  #开始ID转换，会有丢失
  
  gene = bitr(logFC$Genes, fromType=geneType, toType="ENTREZID",OrgDb=OrgDb) 
  
  colnames(gene)[which(colnames(gene) == geneType)] = "Genes"
  gene = merge(gene, logFC, by = "Genes")
  gene = gene[order(gene$logFC, decreasing = T),]
  logFC = gene$logFC
  names(logFC) = gene$ENTREZID
  
  re.kegg = gseKEGG(
    logFC,
    keyType  = 'ncbi-geneid',
    organism = organism,
    pvalueCutoff = pCutoff
  )
  kegg = nrow(re.kegg@result)
  if(kegg == 0) {
    cat("给定的基因列表在KEGG数据库中未富集到任何通路。\n")
  } else {
    cat(paste("给定的基因列表在KEGG数据库中富集到", kegg ,"条通路。\n"))
  }
  
  if (organism == 'mmu') {
    organism = 'mouse'
  }
  
  re.reactome = 
    gsePathway(
      logFC,
      organism = ifelse(organism == "hsa", "human", 
                        ifelse(organism == "mmu" | organism == "mm", "mouse",
                               stop("输入物种错误，要求为hsa或mmu\n"))),
      pvalueCutoff = pCutoff)
  
  reactome = nrow(re.reactome@result)
  if (reactome == 0) {
    cat("给定的基因列表在Reactome数据库中未富集到任何通路。\n")
  } else {
    cat(paste("给定的基因列表在Reactome数据库中富集到", reactome ,"条通路。\n"))
  }
  
  if (organism == "hsa") {
    re.gobp = gseGO(logFC,
                    pvalueCutoff = pCutoff,
                    OrgDb = OrgDb,
                    ont = "BP",
                    keyType  = 'ENTREZID')
    re.gocc = gseGO(logFC,
                    pvalueCutoff = pCutoff,
                    OrgDb = OrgDb,
                    ont = "CC",
                    keyType  = 'ENTREZID')
    re.gomf = gseGO(logFC,
                    pvalueCutoff = pCutoff,
                    OrgDb = OrgDb,
                    ont = "MF",
                    keyType  = 'ENTREZID')
    
    cat(paste0("给定的基因列表在GO数据库中的结果：BP.",
               nrow(re.gobp@result),", CC.",
               nrow(re.gocc@result),", MF.",
               nrow(re.gomf@result),"。\n"))
    
    library(msigdf)
    library(dplyr)
    
    c2<- msigdf.human %>%
      filter(category_code == "c2") %>% 
      select(geneset, symbol) %>% 
      as.data.frame
    
    c2.g = bitr(c2$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    colnames(c2.g)[1] = "symbol"
    c2 = merge(c2, c2.g, by = "symbol")
    c2 = c2[,-1]
    re.GSEA = GSEA(logFC,
                   TERM2GENE= c2,
                   minGSSize= 10,
                   maxGSSize= 500,
                   pvalueCutoff= pCutoff,
                   verbose= FALSE,
                   eps= 0)
    cat(paste("给定的基因列表在MsigDb数据库中富集到", nrow(re.GSEA@result) ,"条通路。\n"))
    
    cat("使用gseaplot2()来可视化结果。\n")
    return(list(
      KEGG = re.kegg, 
      Reactome = re.reactome,
      GOBP = re.gobp,
      GOCC = re.gocc,
      GOMF = re.gomf,
      MSigDb = re.GSEA))
  }
  
  cat("使用gseaplot2()来可视化结果。\n")
  return(list(
    KEGG = re.kegg, 
    Reactome = re.reactome,
    GOBP = re.gobp,
    GOCC = re.gocc,
    GOMF = re.gomf))
}

re = Enrichment(result.deg)


library(ggridges)
library(ggplot2)
library(enrichplot)

pdf(file = "P.KEGG.pdf", width = 10, height = 2.5)
ridgeplot(re$KEGG,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()

pdf(file = "P.Reactome.pdf", width = 10, height = 2.5)
ridgeplot(re$Reactome,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()

pdf(file = "P.GOBP.pdf", width = 10, height = 2.5)
ridgeplot(re$GOBP,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()

pdf(file = "P.GOCC.pdf", width = 10, height = 2.5)
ridgeplot(re$GOCC,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()

pdf(file = "P.GOMF.pdf", width = 10, height = 2.5)
ridgeplot(re$GOMF,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()

pdf(file = "P.MSigDb.pdf", width = 10, height = 2.5)
ridgeplot(re$MSigDb,
          showCategory= 5,
          label_format = 70,
          fill= "p.adjust",
          decreasing= T)
dev.off()
