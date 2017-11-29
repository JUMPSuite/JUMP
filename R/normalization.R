library(limma)
library(MASS)
library(FNN)

##################
## Data loading ##
##################
# Input file is a "raw_..._scan.txt" file generated from itraq.pl main script
## Extract a suffix from the input file name
suffix = tail(unlist(strsplit(inputFile, "/")), 1);
suffix = gsub("raw_", "", suffix)
suffix = gsub("_psm_nonzero.txt", "", suffix)

df = read.table(inputFile, header = TRUE, skip = 1, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "")
names(df) = names(df)[-1]
df = df[, -dim(df)[2]]
colnames(df)[which(colnames(df) == "Outfile")] = "PSM"

####################################
## Normalization at the PSM level ##
####################################
## Unique PSM-level data (= template)
newColNames = NULL
psmData = unique(df[, c(which(colnames(df) == "PSM"),
                        which(colnames(df) == "calcMH"),
                        which(colnames(df) == "RT"),
                        which(colnames(df) == "K_y1"),
                        which(colnames(df) == "R_y1"),
                        which(colnames(df) == "charge"),
                        which(colnames(df) == "endingAA"),
                        which(colnames(df) == "nLabels"),
                         grep("sig", colnames(df)))])
for (i in 1:dim(psmData)[2]) {
  if (length(grep('sig', colnames(psmData)[i])) > 0) {
    s = unlist(strsplit(colnames(psmData)[i], "[.]"))
    newColNames[i] = paste(s[1], " ", "(", tail(s, n = 1), ")", sep = "")
  } else if (length(grep('calcMH', colnames(psmData)[i])) > 0) {
    psmData[, i] = psmData[, i] / as.numeric(psmData$charge)
    newColNames[i] = "mz"
  } else {
    newColNames[i] = colnames(psmData)[i]
  }
}
colnames(psmData) = newColNames

## Normalization
colSig = c(9:dim(psmData)[2])
if (doNormalization == 1) {
  ## Pre-filtering using the defined 'noiseLevel' and 'SNratio'  
  rowSel = NULL
  for (i in colSig) {
    rowSel[[i - 8]] = psmData[, i] > noiseLevel * SNratio
  }
  ## Transform data,
  ## 1. For each PSM, divide reporter intensities by the PSM-mean
  ## 2. Take log2
  psmMean = rowMeans(psmData[, colSig])
  psmData[, colSig] = log2(psmData[, colSig] / psmMean)
  lowest = vector()
  highest = vector()
  for (i in 9:dim(psmData)[2]) {
    sel = sort(psmData[rowSel[[i - 8]], i])
    truncInd = round(length(sel) * pctRemoval / (100 * 2))
    lowest[i - 8] = sel[truncInd]
    highest[i - 8] = sel[length(sel) - truncInd - 1]
  }
  lowest = as.numeric(lowest)
  highest = as.numeric(highest)
  rowSel = which(apply(t(t(psmData[, colSig]) >= lowest & t(psmData[, colSig]) <= highest), 1, sum)
                 == length(colSig))
  if (methodNormalization == 1) { ## trimmed-mean
    meanIntensity = colMeans(psmData[rowSel, colSig])
    normFactor = meanIntensity - mean(meanIntensity)
    psmData[, colSig] = psmData[, colSig] - rep(normFactor, each = nrow(psmData))
    
  } else if (methodNormalization == 2) {  ## trimmed-median
    medianIntensity = apply(psmData[rowSel, colSig], 2, median)
    normFactor = medianIntensity - mean(medianIntensity)
    psmData[, colSig] = psmData[, colSig] - rep(normFactor, each = nrow(psmData))
  }
  psmData[, colSig] = 2^(psmData[, colSig])
  psmData[, colSig] = psmData[, colSig] * rep(psmMean)
}

canRemoveInterference = 0
if (interference_removal == 1) {
  ######################################################
  ## Interference removal (controlled by a parameter) ##
  ######################################################
  
  ## Build linear models for peptides ending with K
  infoK = NULL
  for (charge in 2:3) {
    for (nLabel in 2:3) {
      rowSel = as.numeric(which((psmData$endingAA == "K") &
                                  (psmData$K_y1 > 1000) & 
                                  (psmData$charge == charge) & 
                                  (psmData$nLabels == nLabel) &
                                  is.finite(apply(psmData[, colSig], 1, sum))))
      if (length(rowSel) < 10) {
        next
      }
      x = as.numeric(psmData$K_y1[rowSel])
      y = as.numeric(apply(psmData[, colSig], 1, sum)[rowSel])
      fit = rlm(x, y, maxit = 100)
      tmpInfoK = cbind(rowSel, psmData$RT[rowSel], psmData$mz[rowSel], 
                       matrix(rep(c(charge, nLabel, fit$coefficients), each = length(rowSel)), nrow = length(rowSel)))
      infoK = rbind(infoK, tmpInfoK)
    }
  }
  colnames(infoK) = NULL
  
  ## Build linear models for peptides ending with R
  infoR = NULL
  for (charge in 2:3) {
    for (nLabel in 1:2) {
      rowSel = as.numeric(which((psmData$endingAA == "R") & 
                                  (psmData$R_y1 > 1000) & 
                                  (psmData$charge == charge) & 
                                  (psmData$nLabels == nLabel) &
                                  is.finite(apply(psmData[, colSig], 1, sum))))
      if (length(rowSel) < 10) {
        next
      }
      x = as.numeric(psmData$R_y1[rowSel])
      y = as.numeric(apply(psmData[, colSig], 1, sum)[rowSel])
      fit = rlm(x, y, maxit = 100)
      tmpInfoR = cbind(rowSel, psmData$RT[rowSel], psmData$mz[rowSel], 
                       matrix(rep(c(charge, nLabel, fit$coefficients), each = length(rowSel)), nrow = length(rowSel)))
      infoR = rbind(infoR, tmpInfoR)
    }
  }
  colnames(infoR) = NULL
  
  ## Interference correction
  if (!is.null(infoK) & !is.null(infoR)) {
    for (i in 1:dim(psmData)[1]) {
      corrected = 0
      noiseSaturated = 0
      ## Note that PSMs in infoK and infoR will not be corrected since they do not have 'noisy' peaks in the scans
      if (sum(infoK[, 1] == i) > 0 | sum(infoR[, 1] == i) > 0) {
        next
      }
      ## PSMs ending other than K or R will be ignored
      if (psmData$endingAA[i] != "K" & psmData$endingAA[i] != "R") {
        next
      }
      ## Scans with too many charges or TMT-labeling sites will not be corrected
      if (psmData$charge[i] > 3) {
        next
      }
      if (psmData$endingAA[i] == "K" & psmData$nLabels[i] > 3) {
        next
      }
      if (psmData$endingAA[i] == "R" & psmData$nLabels[i] > 2) {
        next
      }
      ## Scans with zero intensities of both K_y1 and R_y1 will be ignored
      if (psmData$K_y[i] == 0 & psmData$R_y1[i] == 0) {
        next
      }
      ## Correct the interference in each scan
      if (psmData$endingAA[i] == "K") {
        ## Obtain "signal slope"
        signalSlope = infoK[infoK[, 4] == psmData$charge[i] & infoK[, 5] == psmData$nLabels[i], ncol(infoK)]
        if (length(signalSlope) == 0) {
          next
        } else {
          signalSlope = as.numeric(signalSlope[1])
        }
        ## Estimate "noise slope" using KNN (from FNN package)
        kk = min(100, round(dim(infoR)[1] / 10))
        rowkNN = as.numeric(get.knnx(infoR[, c(2, 3)], matrix(c(psmData$RT[i], psmData$mz[i]), nrow = 1), k = kk)$nn.index)
        noiseSlope = mean(infoR[rowkNN, dim(infoR)[2]])
        ## Estimate signal, noise and a correction factor
        estSignal = signalSlope * psmData$K_y1[i]
        estNoise = noiseSlope * psmData$R_y1[i]
        corrFactor = sum(psmData[i, colSig]) / (estSignal + estNoise)
        finalNoise = estNoise * corrFactor / length(colSig) * 2
      } else if (psmData$endingAA[i] == "R") {
        ## Obtain "signal slope"
        signalSlope = infoR[infoR[, 4] == psmData$charge[i] & infoR[, 5] == psmData$nLabels[i], ncol(infoR)]
        if (length(signalSlope) == 0) {
          next
        } else {
          signalSlope = as.numeric(signalSlope[1])
        }
        ## Estimate "noise slope" using KNN (from FNN package)
        kk = min(100, round(dim(infoK)[1] / 10))
        rowkNN = as.numeric(get.knnx(infoK[, c(2, 3)], matrix(c(psmData$RT[i], psmData$mz[i]), nrow = 1), k = kk)$nn.index)
        noiseSlope = mean(infoK[rowkNN, dim(infoK)[2]])
        ## Estimate signal, noise and a correction factor
        estSignal = signalSlope * psmData$R_y1[i]
        estNoise = noiseSlope * psmData$K_y1[i]
        corrFactor = sum(psmData[i, colSig]) / (estSignal + estNoise)
        finalNoise = estNoise * corrFactor / length(colSig) * 2
      }
      ## Manipulation of noise
      if (finalNoise >= 0.5 * min(psmData[i, colSig])) {
        finalNoise = 0.5 * min(psmData[i, colSig])
      }
      ## Correction of reporter intensities by subtracting the "finalNoise"
      psmData[i, colSig] = psmData[i, colSig] - finalNoise
      canRemoveInterference = 1
    }
  } else {
    canRemoveInterference = 0
  }
}
psmData[, colSig] = log2(psmData[, colSig]) ## Log2-transform of normalized/unnormalized intensity

## Correlation matrix of psmData (log2 intensity)
rho = cor(psmData[, colSig])
rho = as.numeric(rho)

## Write the normalized data to a file
psmDataFile = paste("norm", suffix, "psm.txt", sep = "_")
psmDataFile = paste(saveDir, "/", psmDataFile, sep = "")
write.table(psmData[, c(1, colSig)], file = psmDataFile, sep = ";", row.names = FALSE, quote = FALSE)