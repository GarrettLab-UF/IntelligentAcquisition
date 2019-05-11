##remove all R objects from memory
rm(list = ls())

##Simply highlight all and click run
##Jeremy Koelmel - 04/22/2019

########################Other Parameters, can manually change########################
#which inclusion list format are you using?
waters=TRUE #are you using waters (may not be applicable to all instruments)
thermoQE=FALSE #thermo? (formatted for Q-Exactive, some instruments require different formats)
thermoLumos=FALSE #thermo? (formatted for Lumos Fusion, some instruments require different formats)
#the percentile of peak areas will be calculated for samples, and used to determine is the samples > average(blanks)+c*stdev(blanks)
PERCENTILE<-0.75
#CAMERA: do you do adduct searching? This is an extra step and is not necessary for inclusion list generation
CameraAdductSearch<-FALSE
# True or false for adduct searching, does this adduct appear in your chromatography? (eg. metDNA library)
HCO2<-TRUE
CH3COO<-TRUE
NH4<-TRUE
#Exact Mass Search? To reduce the number of included ion only to those in your MS/MS library
ExactMass_LibrarySearch<-TRUE
#How narrow a window do you want for inclusion (minutes)
RTWindow<-0.2
#How to split you data files in the case of multiple inclusion lists? (ONLY CHOOSE ONE)
#Have the first list be the top most abundant ions, second list second top most, etc.
SplitByAbundance<-FALSE
#Have each list evenly distributed throughout the retention time
#eg. RT 1,2,3,4,5,6 if split into two lists would be: 1,3,5 and 2,4,6
SplitByRT<-TRUE

#Parameters that can either be changed manually or using the pop-boxes after selecting all and running
PPM<-15
SNTHR <-10
MINFRAC <-.1
PeakMin<-5
PeakMax<-15
#spell this exactly as such "positive" or "negative"
POLARITY<-"positive"
blankFilterConstant<-5
maxInclusionListSize <-100
Threshold <- 2000

###########################################Inclusion Step##################################################
##Step 1: XCMS Peak picking
##Step 2: Filter by blanks and intensity threshold
##Step 3: Filter by exact mass matches to metabolite library / lipid library (if desired)
##Step 4: Generates inclusion lists

if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets")}
if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk")}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
require(gWidgets)
require(gWidgetstcltk)
options(guiToolkit="tcltk") 
Install = ginput(message="first time running? \nShould we attempt to install all packages?\n (software will check if packages exist before installing) \ninput: y or n", title="Install Packages?",icon="question")

if (Install == "y") {
  install.packages("installr")
  library(installr)
  updateR()
  ## Install packages (more than are loaded...)
  #install.packages("colorspace")
  #install.packages("ggplot2")
  source("http://bioconductor.org/biocLite.R")
  if(!requireNamespace("xcms")){biocLite("xcms")}
  #if(!requireNamespace("pandaR")){biocLite("pandaR", version = "3.8")}
  #if(!requireNamespace("faahKO")){biocLite("faahKO")}
  #if(!requireNamespace("MSnbase")){biocLite("MSnbase")}
  #if(!requireNamespace("pander")){biocLite("pander")}
  if(!requireNamespace("CAMERA")){biocLite("CAMERA")}
  #if(!requireNamespace("limma")){biocLite("limma")}
}


## load packages
library(xcms)
# library(faahKO)
# library(RColorBrewer)
# library(pander)
# library(magrittr)
library(MSnbase)
library(CAMERA)
library(dplyr)

#inputs: Import directory with mzML files
d.in <- choose.dir(caption="Directory with mzML files (do not select a file)\nshould contain atleast 1 file ending in neg.mzML or pos.mzML\nas well as 1-3 files'blank' somewhere in the name \n and ending in _neg or _pos.mzML")
if (CameraAdductSearch==TRUE){
  AdductFile<-choose.files(caption="import .csv file with adducts to search against \nto be used in CAMERA annotation",multi=FALSE)
}
if (ExactMass_LibrarySearch==TRUE) {
  DirExactMassLibrary<-choose.files(caption="import .csv file with metabolites/lipids to search against \nto be used in exact mass matching \nand to reduce the size of the inclusion list",multi=FALSE)
}
ExportName<-ginput(message="What should the export file be called? \nDo not include extension (.csv), periods, or special characters", title="default",icon="question")
setwd(d.in)
POLARITY = ginput(message="What is the polarity? \ninput example: positive OR negative, spell exactly, case sensitive", title="polarity (positve OR negative)",icon="question")
default = ginput(message="use default parameters? \ninput: y or n", title="default",icon="question")


######################interface for parameters#######################
if (default=="n"){
  PPM = ginput(message="What is mass accuracy in ppm for peak picking? \ninput: numeric, example 15", title="mass accuracy (ppm)",icon="question")
  PeakMin = ginput(message="What is your minimum peak width (seconds)? \ninput: numeric, example 2", title="min Peak Width (seconds)",icon="question")
  PeakMax = ginput(message="What is your maximum peak width (seconds)? \ninput: numeric, example 30", title="max Peak Width (seconds)",icon="question")
  SNTHR = ginput(message="What do you want for the signal to noise thresold? \nDefault 10, input: numeric", title="signal to noise threshold",icon="question")
  MINFRAC = ginput(message="What fraction of samples must have peaks? \ninput example: 0.5", title="signal to noise threshold",icon="question")
  blankFilterConstant = ginput(message="What is the number (c) to multiply to the blank standard deviation \nfor which average(sample)-(average(blank)+c*stdev(blank)) must be > 0? \ndefault = 10", title="Blank Subtraction",icon="question")
  Threshold = ginput(message="What is the maximum extracted chromatogram intensity to be consider a peak \nall peaks below this intensity will be removed", title="Threshold Peak Filter",icon="question")
  maxInclusionListSize = ginput(message="What is the maximum number of ions contained on each inclusion list?", title="Exclusion List Max Size",icon="question")
  PPM = as.numeric(PPM)
  PeakMin = as.numeric(PeakMin)
  PeakMax = as.numeric(PeakMax)
  SNTHR = as.numeric(SNTHR)
  MINFRAC = as.numeric(MINFRAC)
  blankFilterConstant = as.numeric(blankFilterConstant)
  maxInclusionListSize = as.numeric(maxInclusionListSize)
  }



# ------------------xcms OLD------------------
# 1.  peak detection
files <- dir(d.in, pattern = '(?i)mzml', recursive = TRUE, full.names = TRUE)####(?i)????????????????????????;"recursive"????????????????????????????????????
numOfSamples <- length(list.files(path=d.in, pattern="(NEG)+|(POS)+", ignore.case=FALSE))
numOfBlanks <- length(list.files(path=d.in, pattern="+blank", ignore.case=FALSE))
if((numOfSamples-numOfBlanks)<1){
  tkmessageBox(title = "An error has occured!",
               message = "Please have atleast one sample (should have 'POS' or 'NEG' in the name\n and not have 'blank' in the name", icon = "error", type = "ok")
}
#Error testing
#error message
if ((numOfSamples<2)||(length(files)<2)){
  tkmessageBox(title = "An error has occured!",
               message = "Please have 2 or more .mzML files in the input folder(s)\nMake sure files have 'pos' or 'neg' in the name\nMake sure blanks have 'blank' in the name", icon = "error", type = "ok")
}

xset <- xcmsSet(files, method = 'centWave', ppm = PPM, snthr = SNTHR, peakwidth = c(PeakMin,PeakMax))

# 2.   grouping 1
xset <- group(xset, minfrac = MINFRAC)

# RT correction obiwarp
xset2 <- retcor(xset, method = 'obiwarp') ## error, change to below--error too, delete QC10022
xset2 <- retcor(xset, method = 'peakgroups', plottype = 'deviation') ##  NO plottype Error!

xset2<-xset


# 3.  grouping 2
xset2 <- group(xset2, bw = 10, mzwid = 0.015, minfrac = MINFRAC)

# 4.   filling gaps
#xset3 <- fillPeaks(xset2)
# groupmat <- xcms::groups(xset3)


# 6.   --------------- CAMERA annotation -----------------

xa <- xsAnnotate(xset2, polarity= POLARITY, nSlaves = 1)
xa <- groupFWHM(xa)
xa <- findIsotopes(xa)
if (CameraAdductSearch==TRUE) {
  rules.camera <- read.csv(AdductFile)
  xa <- findAdducts(xa, rules = rules.camera, polarity = polarity)
  xa <- findAdducts(xa, rules = NULL, polarity = polarity)
}

peaklist.anno <- cbind('name' = groupnames(xset2), #change into xset2 since not using filling gaps
                        getPeaklist(xa, intval = "into")) # into is peakarea without baseline correct, intb with baseline correct(result NA because of the incorrect baseline removement), maxo is intensity
colnames(peaklist.anno)[c(2, 5)] <- c('mzmed', 'rtmed')
#empty matrix to fill with values from peak list
PeakList<-matrix(0,nrow(peaklist.anno),ncol(peaklist.anno)+4)
PeakList[1:nrow(peaklist.anno),1:ncol(peaklist.anno)]<-as.matrix(peaklist.anno)
#select column names with blank
Blanks<-select(peaklist.anno,contains("blank"))

#select samples
Samples<-select(peaklist.anno,matches("(POS)|(NEG)"),-contains("blank"))

#create an average and stdev of blanks and filter by blanks
if (ncol(Blanks)!=0) {
  if (ncol(Blanks)>2) {
  for (i in 1:nrow(Blanks)) { 
    Blanks[i,is.na(Blanks[i,])]<-0
    PeakList[i,ncol(PeakList)-3]<-sum(Blanks[i,])/length(Blanks[i,]) #average blank
    PeakList[i,ncol(PeakList)-2]<-sd(Blanks[i,])/length(Blanks[i,]) #standard deviation blank
    Samples[i,is.na(Samples[i,])]<-0
    PeakList[i,ncol(PeakList)-1]<-(sum(Samples[i,])/length(Samples[i,]))
    SamplePercentile<-as.numeric(quantile(Samples[i,], PERCENTILE))
    #percentile(samples)-(average(blanks)+c*standardDeviation(blanks))
    PeakList[i,ncol(PeakList)]<-SamplePercentile-(as.numeric(PeakList[i,(ncol(PeakList)-3)])+blankFilterConstant*as.numeric(PeakList[i,(ncol(PeakList)-2)]))
  }
  } else {
    #in the case there is only one blank
    for (i in 1:nrow(Blanks)) {
      Blanks[i,is.na(Blanks[i,])]<-0
      #average blanks
      PeakList[i,ncol(PeakList)-3]<-sum(Blanks[i,])/length(Blanks[i,])
      Samples[i,is.na(Samples[i,])]<-0
      #average samples
      PeakList[i,ncol(PeakList)-1]<-(sum(Samples[i,])/length(Samples[i,])) 
      #sample percentile
      SamplePercentile<-as.numeric(quantile(Samples[i,], PERCENTILE))
      # in the case of less than 3 samples, instead of standard deviation use the average and look for a FC difference between samples and blanks
      PeakList[i,ncol(PeakList)]<-SamplePercentile-(as.numeric(ncol(PeakList)-3)+blankFilterConstant*as.numeric(ncol(PeakList)-3))
    }
  }
} else {
  for (i in 1:nrow(Samples)) {
    Samples[i,is.na(Samples[i,])]<-0
    PeakList[i,ncol(PeakList)-1]<-(sum(Samples[i,])/length(Samples[i,])) 
  }
}

peaklist.anno<-cbind(peaklist.anno,PeakList[,ncol(PeakList)-1])
names(peaklist.anno)[ncol(peaklist.anno)]<-"Avg_Samples"

###########################subset data by threshold & Blanks#####################
if (ncol(Blanks)!=0) {
  #Maybe issue, changed PeakList[,20] to PeakList[,ncol(PeakList)]
  peaklist.anno.filtered<-peaklist.anno[as.logical(as.numeric(as.numeric(PeakList[,ncol(PeakList)-1])>as.numeric(Threshold))*as.numeric(as.numeric(PeakList[,ncol(PeakList)])>0)),]
} else {
  peaklist.anno.filtered<-peaklist.anno[as.logical(as.numeric(as.numeric(PeakList[,ncol(PeakList)-1])>as.numeric(Threshold))),]
}

#############subset by exact mass hits##############################
if (ExactMass_LibrarySearch==TRUE) {
  ExactMassLibrary<-read.csv(DirExactMassLibrary)
  
  # Reduce the library to certain adducts and a given polarity
  if (POLARITY=="negative") {
    ExactMassLibrary<-ExactMassLibrary[ExactMassLibrary[,4]=="negative",]
    if (HCO2==FALSE) {
      ExactMassLibrary<-ExactMassLibrary[ExactMassLibrary[,2]!="[M+HCO2]-",]    
    }
    if (CH3COO==FALSE) {
      ExactMassLibrary<-ExactMassLibrary[ExactMassLibrary[,2]!="[M+CH3COO]-",]    
    }
  }
  if (POLARITY=="positive") {
    ExactMassLibrary<-ExactMassLibrary[ExactMassLibrary[,4]=="positive",]
    if (NH4==FALSE) {
      ExactMassLibrary<-ExactMassLibrary[ExactMassLibrary[,2]!="[M+NH4]+",]    
    }
  }
  
  ExactMassLibrary <- as.matrix(ExactMassLibrary)
  #Exact mass matching for Precursor ExactMassLibraryrary to feature table m/z's
  #Need to create a constant to multiply to  mass to get a search tolerance in Da's per mass
  PPM_CONST <- (10^6 + PPM) / 10^6
  
  NumExactMassLibraryMZ <- as.numeric(ExactMassLibrary[,3])
  peaklist.anno.filtered<-cbind(peaklist.anno.filtered,0)
  for (i in 1:nrow(peaklist.anno.filtered)) {
    NumData <- as.numeric(as.character(peaklist.anno.filtered[i,2]))
    DaTolerance<-NumData*PPM_CONST - NumData
    TempID <- ExactMassLibrary[(NumData-DaTolerance < NumExactMassLibraryMZ) & (NumExactMassLibraryMZ < NumData+DaTolerance), 5]
    TempID<-as.character(paste(TempID,collapse=" & "))
    peaklist.anno.filtered[i,ncol(peaklist.anno.filtered)]<-TempID
  }
  names(peaklist.anno.filtered)[ncol(peaklist.anno.filtered)]<-"Metabolite_ID"
  peaklist.anno.filtered.anno<-peaklist.anno.filtered[peaklist.anno.filtered[,ncol(peaklist.anno.filtered)]!="",]
  #sort the filtered and annotated table by most abundant ions or RT
  if (SplitByAbundance==TRUE) {
    peaklist.anno.filtered.anno$Avg_Samples<-as.numeric(as.character(peaklist.anno.filtered.anno$Avg_Samples))
    peaklist.anno.filtered.anno<-arrange(peaklist.anno.filtered.anno, desc(Avg_Samples))
  } else {
    peaklist.anno.filtered.anno$rtmed<-as.numeric(as.character(peaklist.anno.filtered.anno$rtmed))
    peaklist.anno.filtered.anno<-arrange(peaklist.anno.filtered.anno, rtmed)
  }
}


#sort the filtered table by most abundant ions or RT
if (SplitByAbundance==TRUE) {
  peaklist.anno.filtered$Avg_Samples<-as.numeric(as.character(peaklist.anno.filtered$Avg_Samples))
  peaklist.anno.filtered<-arrange(peaklist.anno.filtered, desc(Avg_Samples))
} else {
  peaklist.anno.filtered$rtmed<-as.numeric(as.character(peaklist.anno.filtered$rtmed))
  peaklist.anno.filtered<-arrange(peaklist.anno.filtered, rtmed)
}

#determine inclusion list format based on instrument type 
if (waters==TRUE) {
  #Inclusion List Format
  peaklist.anno.filtered.incl<-matrix(-1,nrow(peaklist.anno.filtered),8)
  peaklist.anno.filtered.incl[,1]<-peaklist.anno.filtered[,2]
  peaklist.anno.filtered.incl[,5]<-peaklist.anno.filtered[,5]
  if (ExactMass_LibrarySearch==TRUE) {  
    peaklist.anno.filtered.anno.incl<-matrix(-1,nrow(peaklist.anno.filtered.anno),8)
    peaklist.anno.filtered.anno.incl[,1]<-peaklist.anno.filtered.anno[,2]
    peaklist.anno.filtered.anno.incl[,5]<-peaklist.anno.filtered.anno[,5]
  }
}
if (thermoQE==TRUE) {
  #header for for thermo formatted exclusion list (Q-Exactive)
  peaklist.anno.filtered.incl = matrix("",nrow(peaklist.anno.filtered),9)
  colnames(peaklist.anno.filtered.incl) = c("Mass [m/z]", "Formula [M]", "Formula type", "Species", "CS [z]", "Polarity", "Start [min]", "End [min]", "Comment")
  #Fill in m/z and RT values for first exclusion list
  peaklist.anno.filtered.incl[,1] = peaklist.anno.filtered[,2] #m/z
  peaklist.anno.filtered.incl[,7] = (as.numeric(peaklist.anno.filtered[,5])/60)-(RTWindow/2) #rt min (minutes)
  peaklist.anno.filtered.incl[,8] = (as.numeric(peaklist.anno.filtered[,5])/60)+(RTWindow/2) #rt max (minutes)
  if (ExactMass_LibrarySearch==TRUE) {  
    #header for for thermo formatted exclusion list (Q-Exactive)
    peaklist.anno.filtered.anno.incl = matrix("",nrow(peaklist.anno.filtered.anno),9)
    colnames(peaklist.anno.filtered.anno.incl) = c("Mass [m/z]", "Formula [M]", "Formula type", "Species", "CS [z]", "Polarity", "Start [min]", "End [min]", "Comment")
    #Fill in m/z and RT values for first exclusion list
    peaklist.anno.filtered.anno.incl[,1] = peaklist.anno.filtered.anno[,2] #m/z
    peaklist.anno.filtered.anno.incl[,7] = (as.numeric(peaklist.anno.filtered.anno[,5])/60)-(RTWindow/2) #rt min (minutes)
    peaklist.anno.filtered.anno.incl[,8] = (as.numeric(peaklist.anno.filtered.anno[,5])/60)+(RTWindow/2) #rt max (minutes)
  }
}
#***BETA, NOT SURE OF EXACT FORMAT (used exclusion format)
if (thermoLumos==TRUE) {
  #Inclusion List Format
  peaklist.anno.filtered.incl<-matrix("",nrow(peaklist.anno.filtered),4) #create empty matrix for lumos (n rows by 3 columns)
  colnames(peaklist.anno.filtered.incl) = c("m/z","Name","t start (min)","stop (min)")
  peaklist.anno.filtered.incl[,1]<-peaklist.anno.filtered[,2] #add m/z values
  #add start and end RT (min)
  peaklist.anno.filtered.incl[,3]<-(as.numeric(peaklist.anno.filtered[,5])/60)-(RTWindow/2) #rt min (minutes)
  peaklist.anno.filtered.incl[,4]<-(as.numeric(peaklist.anno.filtered[,5])/60)+(RTWindow/2) #rt max (minutes)
  if (ExactMass_LibrarySearch==TRUE) {  
    #Inclusion List Format
    peaklist.anno.filtered.anno.incl<-matrix("",nrow(peaklist.anno.filtered.anno),4) #create empty matrix for lumos (n rows by 3 columns)
    colnames(peaklist.anno.filtered.anno.incl) = c("m/z","Name","t start (min)","stop (min)")
    peaklist.anno.filtered.anno.incl[,1]<-peaklist.anno.filtered.anno[,2] #add m/z values
    #add start and end RT (min)
    peaklist.anno.filtered.anno.incl[,3]<-(as.numeric(peaklist.anno.filtered.anno[,5])/60)-(RTWindow/2) #rt min (minutes)
    peaklist.anno.filtered.anno.incl[,4]<-(as.numeric(peaklist.anno.filtered.anno[,5])/60)+(RTWindow/2) #rt max (minutes)
  }
}

#the number of lists which the inclusion list will be split into
nFilteredLists.filtered<-ceiling(nrow(peaklist.anno.filtered.incl)/maxInclusionListSize)
StepSize.filtered<-nFilteredLists.filtered
#How large will each individual inclusion list be (may remove some start and end ions)
SplitSize.filtered<-floor(nrow(peaklist.anno.filtered.incl)/StepSize.filtered)
dir.create("InclusionLists_Filtered")

if (ExactMass_LibrarySearch==TRUE) {
  nFilteredLists.filtered.anno<-ceiling(nrow(peaklist.anno.filtered.anno.incl)/maxInclusionListSize)
  StepSize.filtered.anno<-nFilteredLists.filtered.anno
  SplitSize.filtered.anno<-floor(nrow(peaklist.anno.filtered.anno.incl)/StepSize.filtered.anno)
  dir.create("InclusionLists_Annotated")
}

#Export the inclusion lists after filtering (sorted by most abundant ions)
if (SplitByAbundance==TRUE) {
  start<-1
  for (i in 1:nFilteredLists.filtered) {
    if (i == nFilteredLists.filtered) {
      tempPeakList<-peaklist.anno.filtered.incl[start:nrow(peaklist.anno.filtered.incl),]
      if (waters==TRUE){
        write.table(tempPeakList, paste("InclusionLists_Filtered/",ExportName,"_Filtered_Incl_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
      } else {
        write.table(tempPeakList, paste("InclusionLists_Filtered/",ExportName,"_Filtered_Incl_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      }  
    } else {
      tempPeakList<-peaklist.anno.filtered.incl[start:(i*SplitSize.filtered),]
      start<-i*SplitSize.filtered+1
      if (waters==TRUE){
        write.table(tempPeakList, paste("InclusionLists_Filtered/",ExportName,"_Annotated_Filtered_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
      } else {
        write.table(tempPeakList, paste("InclusionLists_Filtered/",ExportName,"_Annotated_Filtered_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      }  
    }
  }
  #export the inclusion list after exact mass searching
  if (ExactMass_LibrarySearch==TRUE) {
    start<-1
    for (i in 1:nFilteredLists.filtered.anno) {
      if (i == nFilteredLists.filtered.anno) {
        tempPeakList<-peaklist.anno.filtered.anno.incl[start:nrow(peaklist.anno.filtered.anno.incl),]
        if (waters==TRUE){
          write.table(tempPeakList, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
        } else {
          write.table(tempPeakList, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
        }
      } else {
        tempPeakList<-peaklist.anno.filtered.anno.incl[start:(i*maxInclusionListSize),]
        start<-i*maxInclusionListSize+1
        if (waters==TRUE){
          write.table(tempPeakList, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
        } else {
          write.table(tempPeakList, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
        }  
      }
    }
  }
}

#Export the inclusion lists after filtering (split so each inclusion list covers the same retention time ranges)
#This will optimally reduce the density of ions on the list which have the same retention time
if (SplitByRT==TRUE) {
  for (i in 1:nFilteredLists.filtered){
    temporaryInc<-peaklist.anno.filtered.incl[1:SplitSize.filtered,]
    a<-i
    for (x in 1:SplitSize.filtered) {
      temporaryInc[x,]<-peaklist.anno.filtered.incl[a,]
      a<-a+StepSize.filtered
    }
    if (waters==TRUE){
      write.table(temporaryInc, paste("InclusionLists_Filtered/",ExportName,"_Filtered_Incl_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
    } else {
      write.table(temporaryInc, paste("InclusionLists_Filtered/",ExportName,"_Filtered_Incl_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    }  
  }
  if (ExactMass_LibrarySearch==TRUE) {
    for (i in 1:nFilteredLists.filtered.anno){
      temporaryInc<-peaklist.anno.filtered.anno.incl[1:SplitSize.filtered.anno,]
      a<-i
      for (x in 1:SplitSize.filtered.anno) {
        temporaryInc[x,]<-peaklist.anno.filtered.anno.incl[a,]
        a<-a+StepSize.filtered.anno
      }
      if (waters==TRUE){
        write.table(temporaryInc, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".txt",sep=""), sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
      } else {
        write.table(temporaryInc, paste("InclusionLists_Annotated/",ExportName,"_Annotated_Incl_",i,".csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      }  
    }
  }
}

write.csv(peaklist.anno.filtered, paste(ExportName,"_filtered.csv",sep=""))

if (ExactMass_LibrarySearch==TRUE) {
write.csv(peaklist.anno.filtered.anno, paste(ExportName,"_annotated.csv",sep=""))
}
##Jeremy Koelmel - 04/10/2019
##Step 1: XCMS Peak picking
##Step 2: Camera Annotation
##Step 3: Filter by blanks and threshold (change to maximum or 3rd quartile)
##Step 4: Filter by exact mass matches to metabolite library
##Step 5: Generates inclusion lists
