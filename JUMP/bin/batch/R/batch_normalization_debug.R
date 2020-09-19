# inputs:
method='InternalStandard'
raw_exp='raw_exp.txt'
batch_vector='batch_vector.txt'
output='raw_exp.txt'
internal_standard='internal_standard_vector.txt' # required if (method eq 'InternalStandard')
#----------------------------------------------------------------------------
library(limma)
#source('commandLineArgs.R')

# normalize Batche effects By Internal Standard
normalizeBatchesByInternalStandard=function(tb,bch,std) {
	norm=tb
	goldStd=as.character(std[std[,3]=='standard',2])
	for (i in 1:nrow(std)) {
		if (std[i,3] == 'non_standard') {
			bchStandard=as.character(std[i,2])
			f=tb[,goldStd]-tb[,bchStandard]
			crt=as.character(std[i,1])
			#smp=grep(crt,bch)# bug here
			smp=which(bch==crt)# corrected on 6/19/19
			norm[,smp]=tb[,smp]+f
		}
	}
	return(norm)
}

#-----------------------------------------------------------------------------
# To parse R parameters
# Credit to Nurcan Tuncbag from MIT
#args=(commandArgs(TRUE))
for (e in commandArgs(T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  var = ta[[1]][1]
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp,',')[[1]]
    }
    assign(var,temp)
    cat("assigned ",var," the value of |",temp,"|\n")
  } else {
    var_fields = strsplit(var,'-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var,TRUE)
    cat("assigned ",var," the value of TRUE\n")
  }
}

#-----------------------------------------------------------------------------

tb=read.table(raw_exp,head=T,sep="\t",row.names=1)
smallValue=1
tb=log10(tb+smallValue)

batch_vector=read.table(batch_vector,head=F)
bch=as.vector(batch_vector[,1])

if (method == 'InternalStandard') {
	standard_vector=read.table(internal_standard,head=F)
	norm=normalizeBatchesByInternalStandard(tb,bch,standard_vector)
} else{ if (method == 'LinearModel') {
	norm=removeBatchEffect(tb,bch)
}}

write.table(10^norm,file=output,sep="\t",quote=F)
