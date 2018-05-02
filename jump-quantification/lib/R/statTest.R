library(limma)

##########################
## Define a subroutine  ##
##########################
limmaTest <- function(inputData, outputLevel, fdrMethod, saveDir, suffix, comparison, comparisonName, columnNames, selColNames, design, contMatrix) {
  fit = lmFit(inputData[, which(colnames(inputData) %in% selColNames)], design)
  fit = contrasts.fit(fit, contMatrix)
  fit = eBayes(fit)
  result = topTable(fit, genelist = inputData[[outputLevel]], n = nrow(inputData), adjust = fdrMethod, sort = "none")
  ## Change column names of the result table
  if (nGroups == 2) {
    colnames(result)[which(names(result) == "logFC")] = paste("Log2Fold(", columnNames, ")", sep = "")
  } else if (nGroups > 2) {
    ind = grep("group", colnames(result))
    for (k in 1:length(ind)) {
      colnames(result)[ind[k]] = paste("Log2Fold(", columnNames[k], ")", sep = "")
    }
  }
  result$B = NULL
  result$AveExpr = NULL
  result$t = NULL
  result$F = NULL
  colnames(result)[which(names(result) == "P.Value")] = paste("p-value", "_", comparisonName[[1]], "_", comparison[[1]], sep = "")
  colnames(result)[which(names(result) == "adj.P.Value")] = paste("FDR", "_", comparisonName[[1]], "_", comparison[[1]], sep = "")
  colnames(result)[which(names(result) == "adj.P.Val")] = paste("FDR", "_", comparisonName[[1]], "_", comparison[[1]], sep = "")
  colnames(result)[1] = outputLevel
  if (outputLevel == "PSM") {
    suffix = paste(suffix, "psm.txt", sep = "_")
  } else if (outputLevel == "Peptide") {
    suffix = paste(suffix, "pep.txt", sep = "_")
  } else if (outputLevel == "Protein") {
    suffix = paste(suffix, "prot.txt", sep = "_")
  } else if (outputLevel == "Site") {
    suffix = paste(suffix, "site.txt", sep = "_")
  }
  resultFile = paste(comparisonName, suffix, sep = "_")
  resultFile = paste(saveDir, "/", resultFile, sep = "")
  write.table(result, file = resultFile, sep = "\t", row.names = FALSE, quote = FALSE)
}

## From isobar package
fitCauchy <- function(x) {
  cauchy.fit <- function(theta,x){
    -sum(dcauchy(x,location=theta[1],scale=theta[2],log=TRUE),na.rm=T)
  }
  good <- !is.na(x) & !is.nan(x)
  x = x[good]
  x = x[x > quantile(x, 0.25) & x < quantile(x, 0.75)]
  theta.start <- c(median(x),IQR(x)/2)
  res <- nlminb(theta.start,cauchy.fit, x = x,lower=c(-10,1e-20),upper=c(10,10))
}

cauchyTest <- function(inputData, outputLevel, fdrMethod, saveDir, suffix, comparisonName, columnNames, reporters) {
  ## Assumption: there are only two groups, i.e. two reporters
  log2FC = inputData[[reporters[1]]] - inputData[[reporters[2]]]
  fit = fitCauchy(log2FC)
  pval = sapply(log2FC, function(r) {
    if (is.null(fit) || is.na(log2FC))
      return (NA)
    pcauchy(r, location = fit$par[1], scale = fit$par[2], lower.tail = r < fit$par[1])
  })
  pval = 2 * pval
  fdr = p.adjust(pval, method = fdrMethod)
  
  result = data.frame(cbind(inputData[[outputLevel]], log2FC, pval, fdr))
  colnames(result)[1] = outputLevel
  colnames(result)[2] = paste("Log2Fold(", columnNames, ")", sep = "")
  colnames(result)[3] = paste("p-value", "_", comparisonName[[1]], "_", comparison[[1]], sep = "");
  colnames(result)[4] = paste("FDR", "_", comparisonName[[1]], "_", comparison[[1]], sep = "")
  if (outputLevel == "PSM") {
    suffix = paste(suffix, "psm.txt", sep = "_")
  } else if (outputLevel == "Peptide") {
    suffix = paste(suffix, "pep.txt", sep = "_")
  } else if (outputLevel == "Protein") {
    suffix = paste(suffix, "prot.txt", sep = "_")
  } else if (outputLevel == "Site") {
    suffix = paste(suffix, "site.txt", sep = "_")
  }
  resultFile = paste(comparisonName, suffix, sep = "_")
  resultFile = paste(saveDir, "/", resultFile, sep = "")
  write.table(result, file = resultFile, sep = "\t", row.names = FALSE, quote = FALSE)
}

##################
## Data loading ##
##################

## input arguments
## 1. saveDir: indicates the path of 'normalized intensity' files
## 2. suffix: for the names of files to be read and written
## 3. fdrMethod: defines the multiple testing correction method
## 4. comparisons: defines the comparisons and reporters used
## 5. comparisonNames: defined the names of comparison studies

psmFile = paste(saveDir, "/", "norm_", suffix, "_psm.txt", sep = "")
pepFile = paste(saveDir, "/", "norm_", suffix, "_pep.txt", sep = "")
protFile = paste(saveDir, "/", "norm_", suffix, "_prot.txt", sep = "")
siteFile = paste(saveDir, "/", "norm_", suffix, "_site.txt", sep = "")
psmData = read.table(psmFile, header = TRUE, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "")
pepData = read.table(pepFile, header = TRUE, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "")
if (file.exists(protFile)) {
  protData = read.table(protFile, header = TRUE, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "") 
}
if (file.exists(siteFile)) {
  siteData = read.table(siteFile, header = TRUE, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "")
}

reporters = NULL
samples = NULL
for (i in 1:dim(psmData)[2] - 1) {
  s = unlist(strsplit(colnames(psmData)[i + 1], "[.]"))
  reporters[i] = s[1]
  samples[i] = tail(s, n = 1)
}
colnames(psmData)[-1] = samples
colnames(pepData)[-1] = samples
if (exists("protData")) {
  colnames(protData)[-1] = samples  
}
if (exists("siteData")) {
  colnames(siteData)[-1] = samples  
}

######################
## Comparison study ##
######################

for (i in 1:length(comparisons)) {
  ## Retrieve comparison group information
  comparison = unlist(strsplit(gsub("\\s", "", comparisons[[i]]), ":"))
  nGroups = length(comparison)
  groups = list()
  nSamples = 0
  compSamples = NULL
  for (g in 1:nGroups) {
    groups[[g]] = unlist(strsplit(comparison[g], ","))
    nSamples = nSamples + length(groups[[g]])
    compSamples = c(compSamples, groups[[g]])
  }
  
  ## Generate a design matrix (which contains the information of comparison)
  subColNames = colnames(psmData)[which(colnames(psmData) %in% compSamples)]
  design = matrix(0, nrow = nSamples, ncol = nGroups)
  for (g in 1:nGroups) {
    design[which(subColNames %in% groups[[g]]), g] = 1
  }
  colnames(design) = paste("group", seq(1, nGroups), sep = "")
  
  ## Generate a contrast matrix and new column names for the LIMMA result table
  contVec = NULL
  newColNames = NULL
  combMatrix = combn(seq(1, nGroups), 2)
  for (j in 1:ncol(combMatrix)) {
    contVec = c(contVec, paste(paste("group", combMatrix[1, j], sep = ""), paste("group", combMatrix[2, j], sep = ""), sep = "-"))
    newColNames = c(newColNames, paste(comparison[combMatrix[1, j]], "/", comparison[combMatrix[2, j]], sep = ""))
  }
  contMatrix = makeContrasts(contrasts = contVec, levels = design)
  
  if (nGroups == 2 & max(colSums(design)) == 1) {
    ## Cauchy test
    cauchyTest(psmData, "PSM", fdrMethod, saveDir, suffix, comparisonNames[i], newColNames, compSamples)
    cauchyTest(pepData, "Peptide", fdrMethod, saveDir, suffix, comparisonNames[i], newColNames, compSamples)
    if (exists("protData")) {
      cauchyTest(protData, "Protein", fdrMethod, saveDir, suffix, comparisonNames[i], newColNames, compSamples)  
    }
    if (exists("siteData")) {
      cauchyTest(siteData, "Site", fdrMethod, saveDir, suffix, comparisonNames[i], newColNames, compSamples)  
    }
  } else if (nGroups > 2 && max(colSums(design)) == 1) {
    stop("For the comparison of multiple groups, replicates are required")
  } else {
    ## LIMMA running
    ## Statistical testing is performed to the "compSamples"
    limmaTest(psmData, "PSM", fdrMethod, saveDir, suffix, comparisons[i], comparisonNames[i], newColNames, compSamples, design, contMatrix)
    limmaTest(pepData, "Peptide", fdrMethod, saveDir, suffix, comparisons[i], comparisonNames[i], newColNames, compSamples, design, contMatrix)
    if (exists("protData")) {
      limmaTest(protData, "Protein", fdrMethod, saveDir, suffix, comparisons[i], comparisonNames[i], newColNames, compSamples, design, contMatrix)  
    }
    if (exists("siteData")) {
      limmaTest(siteData, "Site", fdrMethod, saveDir, suffix, comparisons[i], comparisonNames[i], newColNames, compSamples, design, contMatrix)  
    }
    
  }
}
