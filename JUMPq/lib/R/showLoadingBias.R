## Show the loading biases over reporter ions
## based on trimmed-mean of log2(intensity)

## Input arguments
## 1. inputFile: the path and name of "raw..._scan.txt" file
## 2. noiseLevel: currently, it is hard-coded as 1000
## 2. SNratio: signal-to-noise ratio for the loading-bias correction (specified in the jump -q parameter file)
## 3. pctRemoval: percentage of highest and lowest intensities to be trimmed for the correction

df = read.table(inputFile, header = TRUE, skip = 1, sep = ";", row.names = NULL, stringsAsFactors = FALSE, comment.char = "")
names(df) = names(df)[-1]
df = df[, -dim(df)[2]]

## Unique PSM-level data (= template)
reporters = NULL
sampleLabels = NULL
psmData = unique(df[, c(which(colnames(df) == "Outfile"), grep("^sig", colnames(df)))])
for (i in 1:dim(psmData)[2] - 1) {
  s = strsplit(colnames(psmData)[i + 1], "[.]")
  reporters[i] = s[[1]][1]
  sampleLabels[i] = tail(s[[1]], n = 1)
}
colnames(psmData)[-1] = reporters

## Loading-biases before normalization will be shown using trimmed-mean of reporter intensities
## Pre-filtering based on the defined 'noiseLevel' and 'SNratio'
selInd = matrix(FALSE, nrow = dim(psmData)[1], ncol = dim(psmData)[2] - 1)
for (i in 1:(dim(psmData)[2] - 1)) {
  selInd[, i] = psmData[, i + 1] > noiseLevel * SNratio
}

## Transform data,
## 1. For each PSM, divide reporter intensities by the PSM-mean
## 2. Take log2
psmMean = rowMeans(psmData[, -1])
psmData[, -1] = log2(psmData[, -1] / psmMean)

## Filtering the most variable intensity values from each reporter
lowest = vector()
highest = vector()
for (i in 1:(dim(psmData)[2] - 1)) {
  sel = sort(psmData[selInd[, i], i + 1])
  truncInd = round(length(sel) * pctRemoval / (100 * 2))
  if (truncInd == 0) {
    lowest[i] = NA
    highest[i] = NA
    next
  }
  lowest[i] = sel[truncInd]
  highest[i] = sel[length(sel) - truncInd - 1]
}
lowest = as.numeric(lowest)
highest = as.numeric(highest)
selInd = which(apply(t(t(psmData[, -1]) >= lowest & t(psmData[, -1]) <= highest), 1, sum) == dim(psmData[, -1])[2])

## Calculate the loading-bias
nRows = length(selInd)
meanIntensity = as.numeric(2 ^ colMeans(psmData[selInd, -1]) * 100)
sdVal = as.numeric(apply(psmData[selInd, -1], 2, sd))
## To express SD as percentage in raw-scale, take the mean deviation of 2^sd and 2^(-sd)
## For example, sdVal = 0.2 for a specific reporter which represent +/-0.2 in log2-scale
## 2^0.2 = 1.1487, 2^(-0.2) = 0.87
## positive variation = 2^0.2 - 1 = 0.1487 (assuming the center is ideally 1 (=2^0))
## negative variation = 1 - 2^(-0.2) = 0.13
## mean variation (in raw-scale) = (0.1487 + 0.13) / 2 = 0.1394
## Let the standard deviation of the reporter be 13.94% (= 0.1394 * 100)
sdIntensity = as.numeric(((2^sdVal - 1) + (1 - 2^(-sdVal))) / 2 * 100)
semIntensity = as.numeric(sdIntensity / sqrt(nRows))
if (length(selInd) > 0) {
  canShowLoadingBias = 1
} else {
  canShowLoadingBias = 0
}
