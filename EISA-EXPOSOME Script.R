#' @title EISA-EXPOSOME Script
#' 

#-----library------#
library(tidyverse)
library(patchwork)
library(xcms)
library(doParallel)

# Rawdata and database file path
rawDatafile <- "G:/test/sample.mzXML"
dbPath <- "G:/test/std_200.xlsx"


# Read rawdata
rawData <- read_files(rawDatafile)

# Targeted extraction of features
featuretable <- targetExtractpeaks(rawData = rawData,
                                   dbData = dbPath,
                                   mz_tol = 0.01,
                                   rt_tol = 1000,   # If rt is not in the database, set rt_tol to greater than or equal to 1000.
                                   rt_window = 60,  # Use to calculate zigzag and peakcor
                                   minpoints = 3,
                                   m = 50,          # Set according to the scan rate
                                   zigzagthre = 0.2,
                                   intthre = 1000,
                                   cores = 4)

# Matching MS/MS,calculating peak cor and marking peak type
Output_table <- matchMS2score(rawData = rawData,
                              target_featuretable = featuretable,
                              dbData = dbPath,
                              mz_tol = 0.01,   # extract eics
                              rt_window = 60,  # use to calculate zigzag and peakcor
                              ms2_tol = 0.02,
                              cores = 4)

# Extracting features using xcms or using other peak extraction software
featuretable_xcms <- pickpeak.XCMS(rawDatafile,ppm=10,mzdiff = 0.01,snthresh = 3,noise = 100,peakwidth = c(1,30))

# Sorting of matching features
filter_result <- selectfeature(Output_table = Output_table,
                               featuretable_xcms,
                               SSM_weight = 0.3,
                               MFR_weight = 0.7,
                               zigzagthre_type = 0.2,
                               candidate_num = 2,
                               mztol_xcms = 0.01,
                               rttol_xcms = 30)

setwd("G:/test")#set workpath to store pdf

ploteics(rawData = rawData,
         Output_table = filter_result,
         dbFile = dbPath,
         rt_tol = 30,
         mz_tol = 0.01)


##-----------------------------function-----------------------------------------------##

read_files <- function(filesPath) {
  options (warn = -1)
  if(grepl("mzML|mzXML|mzData", filesPath)){
    rawData <- MSnbase::readMSData(files = filesPath ,
                                   msLevel. = 1,
                                   mode = "onDisk")
  }else{
   stop("Please check if the format of the raw file is mzML|mzXML|mzData.\n") 
  }
  message("Reading EISA MS RawData...")
  mzs <- xcms::mz(rawData)
  rts <- xcms::rtime(rawData)
  intensities <- xcms::intensity(rawData)
  rawDatas <- list(rawData,mzs,rts,intensities)
  message("DONE.")
  return(rawDatas)
}

targetExtractpeaks <- function(rawData,
                               dbData,
                               mz_tol= 0.01,
                               rt_tol= 1000,   #if database do not have rt,please set rt_tol larger than 1000
                               rt_window = 60,  ##use to calculate zigzag and peakcor
                               intthre= 1000,
                               zigzagthre= 0.2,
                               start_elution = 90,
                               end_elution = 900,
                               m = 20,
                               minpoints = 3,
                               cores = 4) {
  if(str_detect(dbData,".xlsx")){
    dbData <- openxlsx::read.xlsx(dbData) 
  }else if(str_detect(dbData,".csv")){
    dbData <- read.csv(dbData, header = TRUE) 
  }
  # check database file
  must_colnames <- c("PrecursorMZ","ProductMZ","Intensity","NAME","ID")
  if(!all(must_colnames %in% colnames(dbData))){
    stop("please check your database file colnames!!\n")
  }
  id <- na.omit(unique(dbData$ID)) 
  # Calculate the number of cores
  # detectCores() 
  cl <- cores
  # Initiate cluster
  registerDoParallel(cl)
  message("Detecting peaks...")
  result <- foreach::foreach(i = id,
                             .combine = rbind,
                             .packages = "tidyverse",
                             .export = c("peakDectAndZigzagCal","removeshapepeak","zigzag_index",
                                         "targetExtractEIC","get_eic","peakDetection"))%dopar%{
                                           lib_data <- dbData %>% filter(ID == i)
                                           mzmed <- unique(lib_data$PrecursorMZ)
                                           # Targeted extraction of the precursor ion EIC
                                           if(rt_tol >= 1000){
                                             mz_range <- matrix(c(mzmed - mz_tol, mzmed + mz_tol), ncol = 2) 
                                             precursorEIC <- get_eic(object = rawData,
                                                                     mz = mz_range,
                                                                     rt = NULL,
                                                                     aggregationFun = "max")
                                             colnames(precursorEIC) <- c("rt", sprintf("%.4f", mzmed))
                                           }else{
                                             rtmed <- unique(lib_data$RT)  #time is min
                                             precursorEIC <- targetExtractEIC(rawData,
                                                                              mzmed,
                                                                              rtmed,
                                                                              rt_tol = rt_tol,
                                                                              mz_tol = mz_tol)
                                           }
                                           # Remove the time period during which the chromatographic gradient begins to rise and begin to fall  
                                           # for example remove rt < 90s and rt > 900s
                                           precursorEIC <- precursorEIC %>% filter(between(rt,start_elution,end_elution))
                                           # Start detecting peaks and filtering peaks
                                           peaks <- peakDectAndZigzagCal(precursorEIC,
                                                                         zigzagthre = zigzagthre,
                                                                         m = m,
                                                                         rt_window = rt_window,
                                                                         intthre = intthre,
                                                                         minpoints = minpoints)
                                           # output feature informations
                                           if(nrow(peaks)>0){
                                             info <- lib_data[1,] %>% select(-c("ProductMZ","Intensity"))
                                             results <- cbind(info[rep(1,nrow(peaks)),],peaks)
                                           }else{
                                             results <- NULL
                                           }
                                         }
  stopImplicitCluster()
  message("DONE.")
  print(head(result))
  return(result)
}


targetExtractEIC <- function(rawData,
                             mzs,
                             rt,
                             mz_tol=0.01,
                             rt_tol=30){
  mz_range <- matrix(c(mzs - mz_tol, mzs + mz_tol), ncol = 2)
  rt_range <- rt*60 + c(-rt_tol, rt_tol)
  #extraction of the precursor ion EIC
  eics <- get_eic(object = rawData,
                  mz = mz_range,
                  rt = rt_range,
                  aggregationFun = "max")
  colnames(eics) <- c("rt", sprintf("%.4f", mzs))
  return(eics)
}

#' @title get_eic
#' 
#' @description 
#' Similar to the chromatogram function in XCMS for targeted extraction of chromatographic peaks
#' 
#' @param aggregationFun 
#' 
#' @param mz
#' 
#' @param rt
#' 
#' @return A data.frame contains the retention time and intensity
#'  
get_eic <- function(object, mz,rt=NULL,
                    aggregationFun = "max"){
  mzs <- object[[2]]
  rts <- object[[3]]
  inten <- object[[4]]
  if(is.null(rt)==TRUE){
    rt_index <- 1:length(rts)
  }else{
    rt_index <- which(rts>=rt[1]&rts<=rt[2])
  }
  new_mzs <- mzs[rt_index]
  new_inten <- inten[rt_index]
  # combine mz and intensity
  mzandinten <- lapply(1:length(rt_index),function(x){
    data <- data.frame(mz=new_mzs[x],intensity=new_inten[x])
  })
  # Extract the intensity corresponding to the matched mz points and construct the eic
  EICs <- lapply(1:nrow(mz), function(x){
    data <- mz[x,]
    eic <- lapply(1:length(mzandinten), function(i){
      mz_inten <- mzandinten[i][[1]]
      mz_index <- which(mz_inten[,1]>=data[1]&mz_inten[,1]<=data[2])
      if(length(mz_index)==0){
        result <- 0
      }
      if(length(mz_index)==1){
        result <- mz_inten[mz_index,2]
      }
      if(length(mz_index)>1){
        result <- mz_inten[mz_index,]
        if(aggregationFun == "max"){result <- max(result[,2])}
        if(aggregationFun == "mean"){result <- mean(result[,2])}
        if(aggregationFun == "min"){result <- min(result[,2])}
        if(aggregationFun == "sum"){result <- sum(result[,2])}
      }
      return(result)
    }) %>% do.call(rbind,.)
  }) %>% do.call(cbind,.)
  EICs <- cbind.data.frame(rts[rt_index],EICs)
  return(EICs)
}

#' find_peaks function come from github: "https://github.com/stas-g/findPeaks"

peakDetection <- function(eicData,
                          m = 20,
                          intthre=1000){
  find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
  }
  extractedPeaks <- data.frame(trOfPeak= numeric(), peakHeight= numeric(), stringsAsFactors=FALSE)
  peak_detection_result <- find_peaks(eicData[,2], m)
  extractedPeaks <- as.data.frame(cbind(as.matrix(eicData[,1][peak_detection_result]), 
                                        as.matrix(eicData[,2][peak_detection_result])))
  colnames(extractedPeaks) <- c("trOfPeak","peakHeight")
  extractedPeaks <- extractedPeaks[which(extractedPeaks$peakHeight>=intthre),]
  return(extractedPeaks)
}

# zigzag formula
# "https://link.springer.com/article/10.1007/s11306-020-01738-3/tables/2"


zigzag_index <- function(eicdata){
  I <- eicdata[,2]
  N <- nrow(eicdata)
  epi <- max(eicdata)-mean(I)
  FZ <- lapply(2:(N-1), function(x){
    s <- (2*I[x]-I[x-1]-I[x+1])^2
    s
  }) %>% unlist()
  a <- sum(FZ)/(N*epi^2)
  return(a)
}


peakDectAndZigzagCal <- function(eicData,
                                 rt_window = 60,
                                 intthre = 1000,
                                 zigzagthre = 0.6,
                                 m = 20,
                                 minpoints = 3){
  peaks <- peakDetection(eicData,m=m)
  if(nrow(peaks)==0){
    if(max(eicData[,2])>=intthre){
      peakHeight <- max(eicData[,2])
      rt <- eicData[which(eicData[,2]==peakHeight),1]
      peak_eic <- eicData[which(eicData[,1] > (rt-rt_window/2) & eicData[,1] < (rt+rt_window/2)),]
      # The minimum number of points is determined first, and then the zigzag is calculated, 
      # if it is less than the point threshold, there is no need to calculate
      filter_peak <- removeshapepeak(peak_eic,minpoints = minpoints)
      if(max(filter_peak[,2] == 0)){
        results <- data.frame(trOfPeak=0,peakHeight=0,zigzag=100)
      }else{
        #snr <- SNR_cal(filter_peak)
        zigzag <- round(zigzag_index(filter_peak),digits = 4)
        results <- data.frame(trOfPeak=rt,peakHeight=peakHeight,zigzag=zigzag)
      }
    }else{
      results <- data.frame(trOfPeak=0,peakHeight=0,zigzag=100)
    }
  }else{
    if(nrow(peaks)>0){
      results <- lapply(1:nrow(peaks), function(x){
        rt <- peaks[x,1]
        peak_eic <- eicData[which(eicData[,1] > (rt-rt_window/2) & eicData[,1] < (rt+rt_window/2)),]
        # The minimum number of points is determined first, and then the zigzag is calculated, 
        # if it is less than the point threshold, there is no need to calculate zigzag
        filter_peak <- removeshapepeak(peak_eic,minpoints = minpoints)
        if(max(filter_peak[,2]) == 0){
          results <- data.frame(trOfPeak=0,peakHeight=0,zigzag=100)
        }else{
          #snr <- SNR_cal(filter_peak)
          zigzag <- round(zigzag_index(filter_peak),digits = 4)
          results <- data.frame(trOfPeak=rt,peakHeight=peaks[x,2],zigzag=zigzag)
        }
        results
      }) %>% do.call(rbind,.)
    }
  }
  #Filter results where the Zigzag index is greater than the threshold
  results <- results[which(results$zigzag<=zigzagthre),] 
  return(results)
}


find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#' @description
#' Detecting whether the data points that make up 
#' a chromatographic peak satisfy a minimum threshold value

removeshapepeak <- function(eic,minpoints = 3){
  int <- eic[,2]
  # data points with an intensity greater than 0
  point_index <- which(int > 0)
  if(length(point_index) < minpoints){
    eic[,2] <- rep(0,length(int))
  }
  if(length(point_index) == minpoints){
    d <- diff(point_index)
    if(length(which(d ==1)) != (minpoints-1)){
      eic[,2] <- rep(0,length(int)) 
    }
  }
  if(length(point_index) > minpoints){
    n <- length(point_index)-(minpoints-1)
    for (i in 1:n) {
      points <- point_index[i:(i+(minpoints-1))]
      d <- diff(points)
      if(length(which(d ==1))==(minpoints-1)){
        if_ok <- 1
        break
      }else{
        if_ok <- 0
      }
    } 
    if(if_ok == 0){
      eic[,2] <- rep(0,length(int))
    } 
  }
  return(eic)
}

########################################################Annotation###################################################
#' @title matchMS2score  
#' 
#' @description Extract the mass spectra at the chromatographic peak, 
#'              match them to MS/MS in the database, calculate the MFR, SSM, Cor scores 
#'              and output the chromatographic spectra of all ions and the mass spectra of MS/MS
#'  
#' @return   output matched table and graph of the chromatograms and mass spectra
matchMS2score <- function(rawData,
                          target_featuretable,
                          dbData,
                          mz_tol = mz_tol,
                          rt_window = rt_window,
                          ms2_tol = ms2match_tol,
                          cores = 1){
  
  if(str_detect(dbData,".xlsx")){
    dbData <- openxlsx::read.xlsx(dbData) 
  }else if(str_detect(dbData,".csv")){
    dbData <- read.csv(dbData, header = TRUE) 
  }
  # For extraction mass spectrometry
  rts <- rawData[[3]]
  mzs <- rawData[[2]]
  intensities <- rawData[[4]]
  # Remove features with a height of 0
  if(length(which(target_featuretable$peakHeight==0))>0){
    target_featuretable <- target_featuretable[-which(target_featuretable$peakHeight==0),]
  }
  # Associated precursor and fragment ions
  id <- unique(target_featuretable$ID) %>% na.omit()
  # Calculate the number of cores
  cl <- cores
  # Initiate cluster
  registerDoParallel(cl)
  message("Starting annotation...")
  result <- foreach::foreach(i=id,
                             .combine = rbind,
                             .packages = "tidyverse",
                             .export = c("MS2match","getDP","CalpeakCor","targetExtractEIC","get_eic"))%dopar%{
                               feature <- target_featuretable[which(target_featuretable$ID==i),]
                               ms2_lib <- data.frame(dbData[which(dbData$ID==i),]) %>% dplyr::select("ProductMZ","Intensity")
                               names(ms2_lib) <- c("mz","intensity")
                               mzmed <- unique(feature$PrecursorMZ)
                               # match fragments
                               match_results <- lapply(1:nrow(feature), function(x){
                                 data <- feature[x,]
                                 # First,match rt
                                 rt_index <- which(rts==data$trOfPeak)
                                 # Pseudo-secondary spectra corresponding to retention times
                                 ms2_eisa <- data.frame(mz=mzs[rt_index],
                                                        intensity=intensities[rt_index])
                                 names(ms2_eisa) <- c("mz","intensity")
                                 # remove fragments mz larger than precursor mz
                                 ms2_eisa <- ms2_eisa %>% filter(ms2_eisa$mz<=mzmed)
                                 # match
                                 match_result <- MS2match(ms2_eisa,ms2_lib,mz_diff= ms2_tol)
                                 apex_intensity <- data.frame(mz_exp=mzmed,intensity_exp=data$peakHeight,mz.diff=0)
                                 abs_eisa <- rbind(apex_intensity,match_result[,c(1,2,5)])
                                 # Judge if there is an ion match.
                                 if(is.null(match_result)==FALSE){
                                   # Eisa MS/MS data normalisation
                                   match_result[,2] <- match_result[,2]/max(match_result[,2])*100
                                   if(nrow(ms2_lib)>2){
                                     SSM <- getDP(match_result$intensity_exp,
                                                  as.numeric(match_result$intensity_lib)) %>% round(.,2)
                                     # Fragment match score
                                     MFR <- paste0(nrow(match_result),"/",nrow(ms2_lib))
                                   }else{
                                     SSM <- 0
                                     MFR <- paste0(nrow(match_result),"/",nrow(ms2_lib))
                                   }
                                   result <- cbind(name=unique(data$NAME),abs_eisa,SSM=SSM,
                                                   MFR=MFR,apex_rt=round(data$trOfPeak,2),
                                                   zigzag=data$zigzag,DB.ID=data$ID)
                                   # Calculation of peak correlation between precursor and fragment ions
                                   peak_cor <- CalpeakCor(rawData = rawData,
                                                          result = result,
                                                          rt_window = rt_window,
                                                          mz_tol = mz_tol)
                                   output <- cbind(result,peak_cor=peak_cor)
                                 }else{
                                   #There are no matching fragment ions
                                   #SSM =0
                                   # MFR =paste0("0","/",nrow(ms2_lib))
                                   #peak_cor <- 0
                                   #output <- cbind(name=unique(data$NAME),abs_eisa,SSM=SSM,MFR=MFR,apex_rt=round(data$trOfPeak,2),zigzag=data$zigzag,ID=data$ID,peak_cor=peak_cor)
                                   output <- NULL
                                 }
                                
                               }) %>% do.call(rbind,.)
                               return(match_results)
                             }
  stopImplicitCluster()
  message("DNOE.")
  return(result)
}

#calculate peak cor and plot eics
#require("patchwork")
CalpeakCor <- function(rawData,
                       result,
                       rt_window,
                       mz_tol){
  allmz <- result$mz_exp
  apex_rt <- unique(result$apex_rt)/60
  eics <- targetExtractEIC(rawData = rawData,
                           mzs = allmz,
                           rt = apex_rt,
                           mz_tol = mz_tol,
                           rt_tol = rt_window/2) %>% as.data.frame()
  cor_score <- cor(eics[,-1])[,1]%>% round(.,3)
  return(cor_score)
}

#' Title
#' MS2match
#' exp_spec and lib_spec must contain m/z and intensity
#' 
MS2match <- function(exp_spec,lib_spec,mz_diff){
  ms2_exp <- exp_spec
  ms2_lib <- lib_spec
  match_result <- lapply(1:nrow(ms2_lib),function(x){
    data <- ms2_lib[x,]
    fragMZ1 <- data$mz
    fragMZ2 <- ms2_exp$mz
    mz.error <- which(abs(fragMZ1-fragMZ2)<=mz_diff)
    if(length(mz.error)>0){
      match_ms2 <- cbind.data.frame(ms2_exp[mz.error,],data,
                                    round((data$mz-ms2_exp[mz.error,][,1]),4))
      names(match_ms2) <- c("mz_exp","intensity_exp","mz_lib","intensity_lib","mz.diff")
      match_ms2 <- match_ms2[which(abs(match_ms2$mz.diff)==min(abs(match_ms2$mz.diff))),]
      match_ms2 <- match_ms2[which(match_ms2$intensity_exp==max(match_ms2$intensity_exp)),]
    }else{
      match_ms2 <- NULL
    }
  }) %>% do.call(rbind,.)
  return(match_result)
}

#' Title
#' DP
#' 
getDP <- function (exp.int, lib.int) 
{
  exp.weight <- lapply(exp.int, function(x) {
    1/(1 + x/(sum(exp.int) - 0.5))
  }) %>% unlist()
  lib.weight <- lapply(lib.int, function(x) {
    1/(1 + x/(sum(lib.int) - 0.5))
  }) %>% unlist()
  x <- exp.weight * exp.int
  y <- lib.weight * lib.int
  return(sum(x * y)^2/(sum(x^2) * sum(y^2)))
}


selectfeature <- function(Output_table,
                          featuretable_xcms,
                          SSM_weight = 0.3,
                          MFR_weight = 0.7,
                          zigzagthre_type = 0.2,
                          candidate_num = 2,
                          mztol_xcms = 0.01,
                          rttol_xcms =30){
  
  if(mode(featuretable_xcms)=="character"){
    featuretable_xcms <-  openxlsx::read.xlsx(featuretable_xcms)
  }
  xcms_mz <- featuretable_xcms$mz
  xcms_rt <- featuretable_xcms$rt
  id <- unique(Output_table$DB.ID)
  # start filter
  results <- lapply(id, function(x){
    Output_data <- Output_table[which(Output_table$DB.ID==x),]
    # different peaks have different rt
    apex__rt <- unique(Output_data$apex_rt)
    mark_result <- data.frame()
    for (i in apex__rt) {
      Output_data.1 <- Output_data[which(Output_data$apex_rt==i),]
      #calculate weigthed score
      MFR = as.numeric(strsplit(Output_data.1$MFR[1],"/") [[1]][1])/as.numeric(strsplit(Output_data.1$MFR[1],"/") [[1]][2])
      SSM = Output_data.1$SSM[1]
      if(is.null(SSM) | is.na(SSM)){
        SSM = c(0)
      }
      weigthed_score <- round(SSM_weight*SSM + MFR_weight*MFR,3) 
      #if pick by xmcs and mark level
      mz_exp <- Output_data.1$mz_exp[1]
      rt_exp <- Output_data.1$apex_rt[1]
      fit_index <- which(abs(mz_exp-xcms_mz)<=mztol_xcms & abs(abs(rt_exp-xcms_rt)<=rttol_xcms))
      if(length(fit_index)>0){
        type <- c(1)
      }else{
        ifelse(Output_data.1$zigzag[1]<=zigzagthre_type,type <- c(2),type <- c(3))
      }
      #cbind
      a <- cbind(Output_data.1,weigthed_score,type)
      mark_result <- rbind(mark_result,a)
    }
    # sort level then w_score
    sort_w_score <- mark_result[order(mark_result$weigthed_score,decreasing = TRUE),]
    sort_result <- rbind(sort_w_score[which(sort_w_score$type==1),],
                         sort_w_score[which(sort_w_score$type==2),],
                         sort_w_score[which(sort_w_score$type==3),])
    # Output reference features based on candidate_num parameters
    apex__rt <- unique(sort_result$apex_rt)
    if(length(apex__rt)>candidate_num){
      # Compare the weight scores of the last and +1 bits of the candidate feature number
      w_s.1 <- sort_result[which(sort_result$apex_rt == apex__rt[candidate_num]),]
      w_s.2 <- sort_result[which(sort_result$apex_rt == apex__rt[candidate_num+1]),]
      if(unique(w_s.1$weigthed_score) == unique(w_s.2$weigthed_score)){
        feature_num <- candidate_num+1 
      }else{
        feature_num <- candidate_num
      }
      sort_result <- lapply(1:feature_num, function(x){
        w_s <- sort_result[which(sort_result$apex_rt==apex__rt[x]),]
      }) %>% do.call(rbind,.)
    }
    sort_result
  }) %>% do.call(rbind,.)
  return(results)
}

#

pickpeak.XCMS <- function(rawData,ppm,mzdiff,snthresh,noise,peakwidth){
  rawData <- MSnbase::readMSData(rawData,msLevel. = 1,mode = "onDisk")
  cwp <- CentWaveParam(snthresh = snthresh, ppm = ppm,noise = noise,mzdiff = mzdiff,
                       peakwidth = peakwidth)
  peaks <- findChromPeaks(rawData, param = cwp)
  featuretable <- chromPeaks(peaks) %>% as.data.frame()
  return(featuretable)
}


#' Title
#' ploteics
#'

ploteics <- function(rawData,Output_table,dbFile,rt_tol,mz_tol){
  # create workdir to store pictures
  if(!dir.exists("./EISA-EXPOSOME-plotEICs")){
    dir.create("./EISA-EXPOSOME-plotEICs")
  }
  if(str_detect(dbFile,".xlsx")){
    database_data <- openxlsx::read.xlsx(dbFile) 
  }else if(str_detect(dbFile,".csv")){
    database_data <- read.csv(dbFile, header = TRUE) 
  }
  # get annotation table id
  id <- unique(Output_table$DB.ID)
  # Start the cycle
  for (x in id) {
    print(paste("plot ID is ",x))
    data <- Output_table[which(Output_table$DB.ID==x),]
    database_data.1 <- database_data[which(database_data$ID==x),]
    # lib ms/ms
    mz_lib <- database_data.1$ProductMZ
    int_lib <- as.numeric(database_data.1$Intensity)*-1
    lib_msms <- data.frame(mz = mz_lib,int = int_lib,msms_type = "lib")
    # There are multiple features within a rt window
    f_rt <- unique(data$apex_rt)
    for (i in f_rt) {
      data2 <- data[which(data$apex_rt==i),]
      int_exp <- data2$intensity_exp[-1]/max(data2$intensity_exp[-1])*100
      exp_msms <- data.frame(mz = data2$mz_exp[-1],int = int_exp,msms_type = "exp")
      msms <- rbind(exp_msms,lib_msms)
      eics <- targetExtractEIC(rawData = rawData,mzs = data2$mz_exp,rt=unique(data2$apex_rt)/60,
                               mz_tol = mz_tol,rt_tol = rt_tol)%>% as.data.frame()
      #precursor ion eic
      p_eic <- eics[,c(1,2)]
      p_eic <- reshape2::melt(p_eic,id="rt")  #data.frame
      colnames(p_eic) <- c("rt","mz","intens")
      p0 <- p_eic %>% 
        ggplot(aes(x = rt,
                   y = intens)) +
        geom_line(aes(color = mz,
                      group = mz)) +
        geom_point()+
        guides(color = "none") +
        labs(x = "Retention time [s]",
             y = "Intensity",
             title = paste("Precursor mz: ", data2$mz_exp[1],"  zigzag: ",unique(data2$zigzag))) +
        #ggforce::facet_zoom(y = intens < median(result$intensity_exp),zoom.size = 1) +
        theme_minimal()
      
      #p+f
      eics <- reshape2::melt(eics,id="rt")  #data.frame
      colnames(eics) <- c("rt","mz","intens")
      p1 <- eics %>% 
        ggplot(aes(x = rt,
                   y = intens)) +
        geom_line(aes(color = mz,
                      group = mz)) +
        
        guides(color = "none") +
        labs(x = "Retention time [s]",
             y = "Intensity",
             title = paste0("precursor and fragments eics")) +
        #ggforce::facet_zoom(y = intens < median(result$intensity_exp),zoom.size = 1) +
        theme_minimal()
      #ms/ms
     p2 <- ggplot(data=msms,mapping=aes(x=mz,y=int))+
        xlim(10,max(msms$mz))+
        ylim(-106,106)+
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + # add y=0
        geom_point(aes(color = msms$msms_type), size = 1) +
        geom_segment(aes(xend = mz, yend = 0,colour = msms_type), linewidth = 1.2, lineend = "butt")+
        geom_text_repel(
          aes(x = mz, y = int, label = round(mz,4)),
          direction = "y",  # 主要沿y轴方向调整
          nudge_y = ifelse(msms$int > 0, 0.5, -0.5),  # 根据y值正负调整初始偏移
          segment.size = 0.2,
          # segment.color = NA,
          size = 3,
          box.padding = 0.5  # 标签周围的填充
        )+ 
        labs(x = "mz",y = "Relative intensity",title = "Mass spectra of matched fragment ions")+
        theme_bw(base_size = 20,)+
        theme(legend.position = "none")+
        scale_color_manual(values = c("#E41A1C","#377EB8"))+
        annotate("text",x= 25,y= -105,label="Reference",colour="#377EB8",size = 6)+
        annotate("text",x= 25,y= 105,label="Experimental",colour="#E41A1C",size = 6)+
        annotate("text",x= 25,y= 93,label=paste0("MFR:",unique(data2$MFR)),colour="black",size = 5)+
        annotate("text",x= 25,y= -93,label=paste0("SSM:",unique(data2$SSM)),colour="black",size = 5)
      (p0+p1)/p2
      ggsave(paste0("./EISA-EXPOSOME-plotEICs/",unique(data$name),"-",data2$mz_exp[1],"-",i,".pdf"),width = 10, height = 6, dpi = 300)
    } 
  }
  
}
