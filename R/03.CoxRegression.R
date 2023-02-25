warning("Please manually run the code line by line and save the image from the preview box.")
source("00.Functions.R")
CheckPackage(c("readxl"))

# First we need to get the cutoff value
source("01.CalculateCutoff.R")
data.cox = OrginalData.TCGA
data.cox$Name = NULL
data.cox$Group = NULL
data.cox$NewTumorEvent = data.cox$NewTumorTime = NULL
data.cox$TIL = ifelse(data.cox$prec.TIL > res.cut.TCGA$prec.TIL$estimate, "High TIL", "Low TIL")
data.cox$TLS = ifelse(data.cox$prec.TLS > res.cut.TCGA$prec.TLS$estimate, "High TLS", "Low TLS")
data.cox$prec.TLS = data.cox$prec.TIL = NULL
data.cox$Age = ifelse(data.cox$Age < 65, "<65", ">=65")
data.cox$M = ifelse(data.cox$M == "M0" , "M0", "M1 / Mx")
data.cox$N = ifelse(data.cox$N == "N0" , "N0", "N1 / Nx")
data.cox$T = ifelse(data.cox$T == "T1" | data.cox$T == "T2" , "T1 - T2", "T3 - T4")
data.cox$Grade = ifelse(data.cox$Grade == 'G1' | data.cox$Grade == 'G2', "G1-G2", "G3-G4")
data.cox$Stage = NULL
data.cox$invalid = NULL
result.cox = Cox(data.cox, type = "mix", plot = T)

data.cox.XJH = OrginalData.XJH
data.cox.XJH$Name = NULL
data.cox.XJH$Group = NULL
data.cox.XJH$TIL = ifelse(data.cox.XJH$prec.TIL > res.cut.XJH$prec.TIL$estimate, "High TIL", "Low TIL")
data.cox.XJH$TLS = ifelse(data.cox.XJH$prec.TLS > res.cut.XJH$prec.TLS$estimate, "High TLS", "Low TLS")
data.cox.XJH$prec.TLS = data.cox.XJH$prec.TIL = NULL
data.cox.XJH$Age = ifelse(data.cox.XJH$Age < 65, "<65", ">=65")
data.cox.XJH$M = ifelse(data.cox.XJH$M == "M0" , "M0", "M1 / Mx")
data.cox.XJH$N = ifelse(data.cox.XJH$N == "N0" , "N0", "N1 / Nx")
data.cox.XJH$T = ifelse(data.cox.XJH$T == "T1" | data.cox.XJH$T == "T2" , "T1 - T2", "T3 - T4")
data.cox.XJH$Grade = ifelse(data.cox.XJH$Grade == 1 | data.cox.XJH$Grade == 2, "G1-G2", "G3-G4")
result.cox.XJH = Cox(data.cox.XJH, type = "mix", plot = T)
