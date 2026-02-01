# Script to compare ground truth selection tables against detector output selection tables to
# compute various detector performance evaluation metrics
# Note: Assumes ground truth selection tables & prediction selection tables each correspond to a single audio file and
# that the selection table file names contain the corresponding audio file name

library(stringr)
library(warbleR)
library(tuneR)
library(dplyr)
library(ggplot2)
library(pracma)

# Settings -------------------
# Path to ground truth selection tables
anotDir = 'C:/Users/rec297/Documents/GitHub/IISER_YC_AdvancedBioacoustics/Modules/DetClass_ML/Data/Testing/GroundTruthAnnotations/Marine'
# name of column containing ground truth labels
labCol = 'Species' 
# Path to detector output selection tables
detDir = 'C:/Users/rec297/Documents/GitHub/IISER_YC_AdvancedBioacoustics/Modules/DetClass_ML/Data/CustomModel/Predictions'
# name of column containing detection labels
detCol = "Species Code" 
# Path to audio files
audioDir = 'C:/Users/rec297/Documents/GitHub/IISER_YC_AdvancedBioacoustics/Modules/DetClass_ML/Data/Testing/AudioFiles/Marine'
# File extension of audio files
fileExt = '.wav'
# Search nested directories for audio files?
recur = TRUE
# Path to save performance metrics & plots
saveDir = 'C:/Users/rec297/Documents/GitHub/IISER_YC_AdvancedBioacoustics/Modules/DetClass_ML/Data/CustomModel'
# provide class mapping if labels in annotations and predictions don't match; list(c('annotation label', 'prediction label'))
classMap = list(c('Killer Whale','KillerWhale_1000'))

# duration of detector bins (s); if -1, will estimate detector bin size from detector output
binSize = -1
# minimum temp1oral overlap (% of annotation duration) to count a detection window as containing a manually annotated call
minO = 0.5
# either specify channel(s) of interest, or leave as an empty list to analyze all available channels
chans = 1 

# Calculations ------------------------- 

annotation_files <- dir(anotDir,pattern='.txt')
detection_files <- list.files(detDir, pattern = ".txt")
audio_files <- dir(audioDir,pattern=fileExt,recursive=recur)

# Get durations of files in list - ONLY WORKS FOR WAVE/MP3 FILES
durs = numeric()
for (i in 1:numel(audio_files)){
  durs = rbind(durs, (readWave(paste(audioDir,'/',audio_files[i],sep=""),header=TRUE)$samples/readWave(paste(audioDir,'/',audio_files[i],sep=""),header=TRUE)$sample.rate))}

allTimeBins = data.frame(File=character(),Time=numeric())
allALabels = data.frame(File=character(),Channel=numeric(),Time=numeric(),Label=character())
allDLabels = data.frame(File=character(),Channel=numeric(),Time=numeric(),Label=character(),Score=numeric())

for (i in 1:length(audio_files)) { # step through each audio file and evaluate TP/FP/TN/FN
  
  # Find corresponding ground truth annotations and detector output selection table
  matchingDets = str_which(detection_files,str_remove(audio_files[i],fileExt))
  matchingAnnots = str_which(annotation_files,str_remove(audio_files[i],fileExt))
  
  if (!isempty(matchingAnnots)){
    ann_file <- paste(anotDir, annotation_files[matchingAnnots], sep = "/")
    annotations <-   read.table(ann_file, header = TRUE, sep = "\t",check.names = FALSE,quote="",fill=TRUE)
    
    if (dim(annotations)[1]>0){
      # If two views are present in the selection tables, remove annotation type 'waveform'
      if ("View" %in% colnames(annotations) & "Spectrogram 1" %in% annotations$View & "Waveform 1" %in% annotations$View){
        annotations <- subset(annotations, annotations$View == "Spectrogram 1")
      }
      if (!("Delta Time (s)" %in% colnames(annotations))){
        annotations$"Delta Time (s)" = annotations$"End Time (s)" - annotations$"Begin Time (s)"
      }
      # Get rid of rows missing an annotation label in the indicated column
      missAn = which(annotations[,labCol]=='')
      if (!isempty(missAn)){
        annotations = annotations[-missAn,]}
      
      # Sort into ascending start time order, reset row numbers
      annotations = annotations[order(annotations$"Begin Time (s)"),]
      rownames(annotations) = 1:nrow(annotations)
    }
  } else {
    annotations = data.frame()
  }
  
  if (!isempty(matchingDets)){
    # Load detections
    detec_file <-  paste(detDir,detection_files[matchingDets],sep = "/")
    detections <- read.table(detec_file, header = TRUE, sep = "\t",check.names = FALSE,quote="",fill=TRUE)
    
    if (dim(detections)[1]>0){
      # If two views are present in the selection tables, remove annotation type 'waveform'
      if (("View" %in% colnames(detections) & "Spectrogram 1" %in% detections$View & "Waveform 1" %in% detections$View)){
        detections <- subset(detections, detections$View == "Spectrogram 1")
      } # add Delta Time column if not already present
      if (!("Delta Time (s)" %in% colnames(detections))){
        detections$"Delta Time (s)" = detections$"End Time (s)" - detections$"Begin Time (s)"
      }
      if (!("Score" %in% colnames(detections))){
        detections$Score = detections$Confidence
      }
      noCallInd = str_which(detections[[detCol]],'nocall')
      if (!isempty(noCallInd)){
        detections = detections[-noCallInd,]}
      
      if (nrow(detections)>0){
        # Sort into ascending start time order, reset row numbers
        detections = detections[order(detections$"Begin Time (s)"),]
        rownames(detections) = 1:nrow(detections)
      } else {
        detections = data.frame()
      }
    }
  } else {
    detections = data.frame()
  }
  
  
  # Establish time bins consistent with how the detector saw the data 
  if (binSize>0){
    timeBins = seq(0,durs[i],by=binSize)
  }else{
    timeBins = seq(0,durs[i],by=round(min(detections$"End Time (s)"-detections$"Begin Time (s)"),digits=2))}
  timeBins = data.frame(File=rep(audio_files[i],length(timeBins)),Time=timeBins)
  allTimeBins = rbind(allTimeBins,timeBins)
  
  # Trim any annotations that exceed the last full time bin seen by the detector
  # (there may be a tiny bit of data at the end of the file which is not seen by the detector because it's not long enough to be a full time bin)
  tooLong = which(annotations$"End Time (s)">timeBins$Time[nrow(timeBins)])
  annotations$"End Time (s)"[tooLong] = timeBins$Time[nrow(timeBins)]
  
  # Determine full set of channels any annotations or detections exist in
  if (isempty(chans)){
    allChans = sort(unique(c(annotations$Channel,detections$Channel)))
  } else { allChans = chans}
  
  for (j in 1:length(allChans)){
    
    # Find annotations and detections in this channel
    anInd = which(annotations$Channel==allChans[j])
    detInd = which(detections$Channel==allChans[j])
    
    # for each time bin, note existing annotation and/or detection labels
    for (k in 1:(nrow(timeBins)-1)){
      
      # find any annotations which sufficiently overlap with this bin
      whichAnInds = which(annotations$"Begin Time (s)"[anInd]<timeBins$Time[k+1] & annotations$"End Time (s)"[anInd]>timeBins$Time[k]) 
      Alabels = character()
      if (length(whichAnInds)>0){
        for (l in 1:length(whichAnInds)){
          # for each annotation, determine if overlap is sufficient
          overlap = (min(timeBins$Time[k+1],annotations$"End Time (s)"[anInd[whichAnInds[l]]]) - max(annotations$"Begin Time (s)"[anInd[whichAnInds[l]]],timeBins$Time[k]))/annotations$"Delta Time (s)"[anInd[whichAnInds[l]]]
          if (overlap >=minO){
            Alabels = c(Alabels,annotations[anInd[whichAnInds[l]],labCol])
          }
        }
        temp1 = data.frame(File=rep(audio_files[i],length(Alabels)),Channel=rep(j,length(Alabels)),Time=rep(timeBins$Time[k],length(Alabels)),Label=Alabels)
        allALabels = rbind(allALabels,temp1)
        rm(temp1)
      }
      
      # find any detections which sufficiently overlap with this bin
      whichDetInds = which(detections$"Begin Time (s)"[detInd]<timeBins$Time[k+1] & detections$"End Time (s)"[detInd]>timeBins$Time[k]) 
      Dlabels = character()
      Dscores = numeric()
      if (length(whichDetInds)>0){
        Dlabels = c(Dlabels,detections[detInd[whichDetInds],detCol])
        Dscores = c(Dscores,detections$Score[detInd[whichDetInds]])
        
        temp2 = data.frame(File=rep(audio_files[i],length(Dlabels)),Channel=rep(j,length(Dlabels)),Time=rep(timeBins$Time[k],length(Dlabels)),Label=Dlabels,Scores=Dscores)
        allDLabels = rbind(allDLabels,temp2)
        rm(temp2)
      }
      
    }
  }
  
}

# Tally TP/FP/TN/FN for each class and compute performance metrics
thresh = c(0.1,0.15,0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)

# If necessary, update labels according to class mapping
if (length(classMap)>0){
  for (i in 1:length(classMap)){
    ind = str_which(allALabels$Label,classMap[[i]][1])
    allALabels$Label[ind] = classMap[[i]][2]
  }}

allLabels = sort(unlist(unique(c(allALabels$Label,allDLabels$Label))))

for (i in 1:length(allLabels)){
  
  metMat = matrix(nrow=length(thresh),ncol=11)
  colnames(metMat) = c("Thresh",'nCalls','nTP','nFP','nTN','nFN','A','P','R','F1','FPR')
  metMat[,2] = length(which(allALabels$Label==allLabels[i]))
  
  for (j in 1:length(thresh)){
    
    # indices with this detection label and a scores exceeding the confidence threshold
    goodDInds = which(allDLabels$Score>=thresh[j] & allDLabels$Label==allLabels[i]) 
    # indices with this annotation label
    goodAInds = which(allALabels$Label==allLabels[i]) 
    # detections with this label and a confidence score exceeding the threshold that match ground truth
    TP = nrow(inner_join(allDLabels[goodDInds,],allALabels,by=join_by(File,Channel,Time,Label))) 
    # ground truth annotations that do NOT match a detection with this label and confidence score exceeding threshold
    FN = nrow(setdiff(allALabels[goodAInds,],allDLabels[goodDInds,1:4])) 
    # detections with this label and a confidence score exceeding the threshold that do NOT match ground truth
    FP = nrow(setdiff(allDLabels[goodDInds,1:4],allALabels)) 
    # bins that do not contain an annotation of this label and also do not have a detection of this label exceeding the confidence score threshold
    TN = nrow(setdiff(setdiff(allTimeBins,allALabels[goodAInds,c(1,3)]),allDLabels[goodDInds,c(1,3)]))
    
    metMat[j,3] = TP
    metMat[j,4] = FP
    metMat[j,5] = TN
    metMat[j,6] = FN
    metMat[j,7] = round((TP+TN)/(TP+TN+FP+FN),3) #Accuracy
    metMat[j,8] = round(TP/(TP+FP),3) # Precision
    metMat[j,9] = round(TP/(TP+FN),3) # Recall
    metMat[j,10] = round((2*TP)/((2*TP)+FP+FN),3) # F1 Score
    metMat[j,11] = round(FP/(FP+TN),3) # FPR
  }
  
  metMat = as.data.frame(metMat)
  metMat$Thresh = thresh
  write.table(metMat,paste(saveDir,'/',str_remove(allLabels[i],' '),'_PerformanceMetrics.txt',sep=""),sep="\t",row.names=FALSE)
  cat(paste('Label: ',allLabels[i],
            '\nAccuracy: ',as.character(min(metMat$A,na.rm=TRUE)*100),'-',as.character(max(metMat$A,na.rm=TRUE)*100),
            '\nPrecision: ',as.character(min(metMat$P,na.rm=TRUE)*100),'-',as.character(max(metMat$P,na.rm=TRUE)*100),
            '\nRecall: ',as.character(min(metMat$R,na.rm=TRUE)*100),'-',as.character(max(metMat$R,na.rm=TRUE)*100),
            '\nF1: ',as.character(min(metMat$F1,na.rm=TRUE)*100),'-',as.character(max(metMat$F1,na.rm=TRUE)*100),
            '\nFPR: ',as.character(min(metMat$FPR,na.rm=TRUE)*100),'-',as.character(max(metMat$FPR,na.rm=TRUE)*100),'\n',sep=""))
  
  ## Plot performance curves for each species
  # PR curve vs confidence score
  ggplot(metMat,aes(label=Thresh))+
    geom_point(aes(x=R,y=P))+
    geom_path(aes(x=R,y=P))+
    geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
         x='Recall',
         y='Precision')
  ggsave(filename=paste(saveDir,'/',str_remove(allLabels[i],' '),'_PR_conf.png',sep=""))
  
  # ROC curve vs conf
  ggplot(metMat,aes(label=Thresh))+
    geom_point(aes(x=FPR,y=R))+
    geom_path(aes(x=FPR,y=R))+
    geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste('ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
         x='FPR',
         y='Recall')
  ggsave(filename=paste(saveDir,'/',str_remove(allLabels[i],' '),'_ROC.png',sep=""))
  
  
  # Plot precision vs thresh
  ggplot(metMat)+
    geom_point(aes(x=Thresh,y=P))+
    geom_path(aes(x=Thresh,y=P))+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
         x='Threshold',
         y='Precision')
  ggsave(filename=paste(saveDir,'/',str_remove(allLabels[i],' '),'_PvThresh.png',sep=""))
  
  # Plot recall vs thresh
  ggplot(metMat)+
    geom_point(aes(x=Thresh,y=R))+
    geom_path(aes(x=Thresh,y=R))+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
         x='Threshold',
         y='Recall')
  ggsave(filename=paste(saveDir,'/',str_remove(allLabels[i],' '),'_RvThresh.png',sep=""))
}


# Compute model performance across all classes
metMat = matrix(nrow=length(thresh),ncol=11)
colnames(metMat) = c("Thresh",'nCalls','nTP','nFP','nTN','nFN','A','P','R','F1','FPR')
metMat[,2] = length(allALabels$Label)

for (j in 1:length(thresh)){
  
  # indices in allDLabels where detection scores exceeded confidence threshold
  goodInds = which(allDLabels$Score>=thresh[j]) 
  # detections with this label and a confidence score exceeding the threshold that match ground truth
  TP = nrow(inner_join(allDLabels[goodInds,],allALabels,by=join_by(File,Channel,Time,Label))) 
  # ground truth annotations that do NOT match a detection with this label and confidence score exceeding threshold
  FN = nrow(setdiff(allALabels,allDLabels[goodInds,1:4])) 
  # detections with this label and a confidence score exceeding the threshold that do NOT match ground truth
  FP = nrow(setdiff(allDLabels[goodInds,1:4],allALabels)) 
  # bins that do not contain an annotation of this label and also do not have a detection of this label exceeding the confidence score threshold
  TN = nrow(setdiff(setdiff(allTimeBins,allALabels[,c(1,3)]),allDLabels[goodInds,c(1,3)]))
  
  metMat[j,3] = TP
  metMat[j,4] = FP
  metMat[j,5] = TN
  metMat[j,6] = FN
  metMat[j,7] = round((TP+TN)/(TP+TN+FP+FN),3) #Accuracy
  metMat[j,8] = round(TP/(TP+FP),3) # Precision
  metMat[j,9] = round(TP/(TP+FN),3) # Recall
  metMat[j,10] = round((2*TP)/((2*TP)+FP+FN),3) # F1 Score
  metMat[j,11] = round(FP/(FP+TN),3) # FPR
}

metMat = as.data.frame(metMat)
metMat$Thresh = thresh
write.table(metMat,paste(saveDir,'/Overall_PerformanceMetrics.txt',sep=""),sep="\t",row.names=FALSE)
cat(paste('Label: ',allLabels[i],
          '\nAccuracy: ',as.character(min(metMat$A,na.rm=TRUE)*100),'-',as.character(max(metMat$A,na.rm=TRUE)*100),
          '\nPrecision: ',as.character(min(metMat$P,na.rm=TRUE)*100),'-',as.character(max(metMat$P,na.rm=TRUE)*100),
          '\nRecall: ',as.character(min(metMat$R,na.rm=TRUE)*100),'-',as.character(max(metMat$R,na.rm=TRUE)*100),
          '\nF1: ',as.character(min(metMat$F1,na.rm=TRUE)*100),'-',as.character(max(metMat$F1,na.rm=TRUE)*100),
          '\nFPR: ',as.character(min(metMat$FPR,na.rm=TRUE)*100),'-',as.character(max(metMat$FPR,na.rm=TRUE)*100),'\n',sep=""))

## Plot performance curves
# PR curve vs confidence score
ggplot(metMat,aes(label=Thresh))+
  geom_point(aes(x=R,y=P))+
  geom_path(aes(x=R,y=P))+
  geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
       x='Recall',
       y='Precision')
ggsave(filename=paste(saveDir,'/Overall_PR_conf.png',sep=""))

# ROC curve vs conf
ggplot(metMat,aes(label=Thresh))+
  geom_point(aes(x=FPR,y=R))+
  geom_path(aes(x=FPR,y=R))+
  geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title=paste('ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
       x='FPR',
       y='Recall')
ggsave(filename=paste(saveDir,'/Overall_ROC.png',sep=""))


# Plot precision vs thresh
ggplot(metMat)+
  geom_point(aes(x=Thresh,y=P))+
  geom_path(aes(x=Thresh,y=P))+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
       x='Threshold',
       y='Precision')
ggsave(filename=paste(saveDir,'/Overall_PvThresh.png',sep=""))

# Plot recall vs thresh
ggplot(metMat)+
  geom_point(aes(x=Thresh,y=R))+
  geom_path(aes(x=Thresh,y=R))+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
       x='Threshold',
       y='Recall')
ggsave(filename=paste(saveDir,'/Overall_RvThresh.png',sep=""))
