# library
packages <- c('shiny','shinydashboard','xcms','doParallel','tidyverse','reshape2',
               'patchwork','shinyjs','shinycssloaders','RaMS')
lapply(packages, require, character.only = TRUE) 
#--------------------------------Function---------------------------------------#

read_files <- function(filesPath) {
  
  if(grepl("mzML|mzXML|mzData", filesPath)){
    rawData <- MSnbase::readMSData(files = filesPath ,
                                   msLevel. = 1,
                                   mode = "onDisk")
  }else{
    stop("Please check if the format of the raw file is mzML|mzXML|mzData.\n") 
  }
  rawData <- MSnbase::readMSData(files = filesPath ,
                                 msLevel. = 1,
                                 mode = "onDisk")
  mzs <- xcms::mz(rawData)
  rts <- xcms::rtime(rawData)
  intensities <- xcms::intensity(rawData)
  rawDatas <- list(rawData,mzs,rts,intensities)
  return(rawDatas)
}

targetExtractpeaks <- function(rawData,
                               dbData,
                               mz_tol= 0.01,
                               rt_tol= 1000,    # If database do not have rt,please set rt_tol larger than 1000
                               rt_window = 60,  # Use to calculate zigzag and peakcor
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
    stop("please check your database colnames!!")
  }
  id <- na.omit(unique(dbData$ID)) 
  # Calculate the number of cores
  # detectCores() 
  cl <- cores
  # Initiate cluster
  registerDoParallel(cl)
  result <- foreach::foreach(i=id,
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
  dbData <- openxlsx::read.xlsx(dbData)
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
                                 # Read if there is an ion match.
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
                                 #return(output)  
                               }) %>% do.call(rbind,.)
                               return(match_results)
                             }
  stopImplicitCluster() 
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
  #start filter
  results <- lapply(id, function(x){
    Output_data <- Output_table[which(Output_table$DB.ID==x),]
    #different peaks have different rt
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
    #sort level then w_score
    sort_w_score <- mark_result[order(mark_result$weigthed_score,decreasing = TRUE),]
    sort_result <- rbind(sort_w_score[which(sort_w_score$type==1),],
                         sort_w_score[which(sort_w_score$type==2),],
                         sort_w_score[which(sort_w_score$type==3),])
    #Output reference features based on candidate_num parameters
    apex__rt <- unique(sort_result$apex_rt)
    if(length(apex__rt)>candidate_num){
      w_s.1 <- sort_result[which(sort_result$apex_rt==apex__rt[candidate_num]),]
      w_s.2 <- sort_result[which(sort_result$apex_rt==apex__rt[candidate_num+1]),]
      if(unique(w_s.1$weigthed_score)==unique(w_s.2$weigthed_score)){
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
  #featuretable <- featuretable[,-c(2,3,5,6)]
  return(featuretable)
}
########################################shiny#######################################################
#ui
ui <- dashboardPage(
  skin = "green",
  title = "EISA-EXPOSOME",
  dashboardHeader(title = "EISA-EXPOSOME"),
      #Sidebar#######################################################################
  dashboardSidebar(
    sidebarMenu(
      menuItem("Related literature", tabName = "literature", icon = icon("home")),
      #menuItem("Introduction", tabName = "introduction", icon = icon("th")),
      #menuItem("Database", tabName = "Database", icon = icon("cog")),
      #menuItem("Targeted Peak", tabName = "TargetedPeak", icon = icon("chart-line")),
      #menuItem("Annotation", tabName = "Annotation", icon = icon("credit-card")),
      menuItem("EISA-EXPOSOME", tabName = "EISA-EXPOSOME", icon = icon("list-alt")),
      #menuItem("Construct database", tabName = "ConstructDB", icon = icon("credit-card")),
      menuItem("Bioinformatics tools", tabName = "tools", icon = icon("bar-chart"))
    )),
  dashboardBody(
    tabItems(
      #First 
      #tabItem(tabName = "introduction"),
      #Second #######################################################################
      tabItem(tabName = "Database",
              fluidRow(
                  box( titlePanel("Database for EISA-EXPOSOME"),
                       div(style="width:fit-content;width:-webkit-fit-content;width:-moz-fit-content;font-size:150%;margin-top:20px",
                           HTML("<b>EISA-EXPOSOME</b> is a computational tool for exposome targered extraction and annotation based on databases.
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The database contains the following parts:
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>PrecursorMZ</b></a>: After 
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>ProductMZ</b></a>: After 
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Intensity</b></a>: After 
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>RT</b></a>: After 
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>NAME</b></a>: After
                   <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>ID</b> </a>: After")),
                 mainPanel(img(src="database_example.png", height = 330, width = 800)),
                 width = 6,
                 height = 700),
                 box(
                   titlePanel("Database MS/MS (public)"),
                   checkboxGroupInput("db","You can find candidates for experimental and/or predicted spectra.",c("T3DB","NIST14","MONA","HMDB","Norman")),
                   actionButton("view.DB", h3("View"), width = "100%"),
                   numericInput(inputId = "dbmz", label = "compound mz", value = 100,width = 1000),
                   numericInput(inputId = "mz.error", label = "mz window", value = 0.01,width = 1000),
                   actionButton("search.0", h3("Search"), width = "100%"),
                   height = 700)
                
              ),
              fluidRow(
                box(titlePanel("Database Lists (public)"),
                    dataTableOutput("datalist"),
                    width = 12) 
              )),
      #second -end
      
      #Third pick peak################################################################
      tabItem(tabName = "TargetedPeak",
              fluidRow(
               box(
                  id = "Parameters.of.MS1.3",
                  h3("Parameters of EISA-EXPOSOME (MS1)"),
                  numericInput(inputId = "MS1deltaMZ.3", label = "Delta m/z of MS1", value = 0.01),
                  numericInput(inputId = "MS1deltaTR.3", label = "Delta tR of MS1", value = 30),
                  numericInput(inputId = "zizagThre.3", label = "zigzag threshold", value = 0.4),
                  numericInput(inputId = "intThre.3", label = "Intensity threshold", value = 1000),
                  width = 6,
                  height = 400
                ),
                box( id = "Database.input.3",
                     #background = "aqua",
                     h3("Database input"),
                     fileInput(inputId = "db.path.3", label = "Database", multiple = F,accept = ".xlsx"),
                     h3("EISA MS1 data input"),
                     fileInput(inputId = "msRawData.3", label = "RawData", multiple = F,accept = ".mzXML"),
                     ###Start actionButton
                     #h3("If all parameters were set, click the button.", align = "center"),
                     actionButton("start.3", h3("Start"), width = "100%"),
                     width = 6,
                     height = 400
                )
              ),
              fluidRow(
                box(title = "Featuretable",
                    dataTableOutput("dynamic"),
                    downloadButton("down","Download feature table")),
                box(title = "Extract ion chromatograms(EIC)",
                    numericInput(inputId = "ID", label = "Enter the ID of the feature to output its chromatogram", value = 2,width = "100%"),
                    selectInput("datapoints", label = "show data points",choices = list("Y" = "o", "N" = "l"),selected = "Y"),
                    selectInput("col", label = "Specify the colour of the drawing",choices = list("red" = "red", "blue" = "blue","green"="green"),selected = "red"),
                    actionButton("plot.3", h3("Plot"), width = "100%"),  #开始绘制
                    shinycssloaders::withSpinner(plotOutput("plot", width = "1000px",height = "600px"))))),
      #Third pick peak - end
      
      # #Forth Annotation################################################################
      # #annotation  ui
      # tabItem(tabName = "Annotation",
      #         fluidRow(
      #           tabBox(
      #             ##设置EISA MS1 数据参数
      #             tabPanel("Parameters of EISA-EXPOSOME (MS1)",
      #                      #"First tab content",
      #                      h3("Parameters of EISA-EXPOSOME (MS1)"),
      #                      numericInput(inputId = "MS1deltaMZ", label = "Delta m/z of MS1", value = 0.01),
      #                      numericInput(inputId = "MS1deltaTR", label = "Delta tR of MS1", value = 30),
      #                      numericInput(inputId = "rt_window.2", label = "Delta rt window for calculating zigzag and peakcor and ploting EIC", value = 60),
      #                      numericInput(inputId = "zizagThre", label = "zigzag threshold", value = 0.4),
      #                      numericInput(inputId = "intThre", label = "Intensity threshold", value = 1000),
      #                      numericInput(inputId = "MS2deltaMZ", label = "Delta m/z of MS/MS", value = 0.025),
      #                      selectInput(inputId = "cores", label = "The number of cores used for parallel computing", choices = c(1:8),selected = 2)),
      #             ##设置xcms MS1 数据参数
      #             tabPanel("Parameters of XCMS", 
      #                      h3("Parameters of XCMS"),
      #                      numericInput(inputId = "ppm.1", label = "ppm", value = 10),
      #                      numericInput(inputId = "noise.1", label = "noise", value = 100),
      #                      numericInput(inputId = "snthresh.1", label = "snthresh", value = 3),
      #                      numericInput(inputId = "mzdiff.1", label = "mzdiff", value = -0.001),
      #                      numericInput(inputId = "peakwidth.11", label = "peakwidth(min)", value = 1),
      #                      numericInput(inputId = "peakwidth.22", label = "peakwidth(max)", value = 30)),
      #             id = "Parameters.of.MS1",
      #             width = 6,
      #             height = 630
      #           ),
      #           box( id = "Database.input",
      #                h3("Database input"),
      #                fileInput(inputId = "db.path", label = "Database", multiple = F,accept = ".xlsx"),
      #                h3("EISA MS1 data input"),
      #                fileInput(inputId = "msRawData", label = "RawData", multiple = F,accept = ".mzXML"),
      #                h3("Feature input (Processed by xcms, mzmine etc)"),
      #                fileInput(inputId = "feature.xcms", label = "Featuretable(1. mz; 2. rt)", multiple = F),
      #                ###start 
      #                #h3("If all parameters were set, click the button.", align = "center"),
      #                actionButton("start", h3("Start"), width = "100%"),
      #                #textOutput("text.station"),
      #                width = 6,
      #                height = 630
      #           )
      #         ),
      #         fluidRow(
      #           box(title = "AnnotationTable",
      #               dataTableOutput("AnnotationTable"),
      #               downloadButton("down.annotation","Download annotation table")),
      #           box(title = "Extract ion chromatograms(EIC)",
      #               numericInput(inputId = "ID.1", label = "Enter the ID of the feature to output its chromatogram", value = 2,width = "100%"),
      #               numericInput(inputId = "ID.3", label = "Enter the serial number to plot the different chromatographic peaks of the above IDs", value = 1,width = "100%"),
      #               actionButton("plot.start", h3("Plot"), width = "100%"),  #start plot
      #               shinycssloaders::withSpinner(plotOutput("ploteics", width = "1400px",height = "900px"))))
      # ),
      
      #Forth Annotation - end###########################################################
      #tabItem
      ##EISA-EXPOSOME################################################################
      tabItem(tabName = "EISA-EXPOSOME",
              fluidRow(
                tabBox(
                  ##Set EISA MS1 data extraction parameters
                  tabPanel("Parameters of EISA-EXPOSOME (MS1)",
                           #"First tab content",
                           h3("Parameters of EISA-EXPOSOME (MS1)"),
                           column(5,
                                  numericInput(inputId = "MS1deltaMZ.5", label = "Delta m/z of MS1 (Da)", value = 0.01)),
                           column(5,
                                  numericInput(inputId = "zizagThre.5", label = "zigzag threshold", value = 0.2)),
                           column(5,
                                  numericInput(inputId = "intThre.5", label = "Intensity threshold", value = 1000)),
                           column(5,
                                  numericInput(inputId = "MS2deltaMZ.5", label = "Delta m/z of MS/MS (Da)", value = 0.025)),
                           column(5,
                                  numericInput(inputId = "Lower_RT_Limit", label = "Time to start elution (sec)", value = 90)),
                           column(5,
                                  numericInput(inputId = "Upper_RT_Limit", label = "Time to end elution (sec)", value = 900)),
                           column(5,
                                  selectInput(inputId = "cores.5", label = "The number of cores used for parallel computing", choices = c(1:8),selected = 4)),
                           column(5,
                                  numericInput(inputId = "rt_window.3", label = "Delta rt window for calculating zigzag and peakcor and ploting EIC (sec)", value = 60)),
                           column(5,
                                  numericInput(inputId = "pickpeak_window", label = "RT Window for picking peaks (sec)", value = 1000)),
                           column(5,
                                  numericInput(inputId = "steps", label = "Steps for picking peaks", value = 50)),
                           column(5,
                                  numericInput(inputId = "minpoints", label = "Minimum number of data points to form a peak", value = 3))),
                           #column(7,h3("Important notes : "))),
                  #Set filter parameters
                  tabPanel("Parameters of EISA-EXPOSOME to filter result",
                           #"First tab content",
                           h3("Parameters of EISA-EXPOSOME to filter result"),
                           numericInput(inputId = "SSM_weigth", label = "SSM_weigth", value = 0.3),
                           numericInput(inputId = "MFR_weight", label = "MFR_weight", value = 0.7),
                           numericInput(inputId = "Delta_mz_xcms", label = "Delta m/z of MS1 to fit xcms (Da)", value = 0.01),
                           numericInput(inputId = "Delta_rt_xcms", label = "Delta rt of MS1 to fit xcms (sec)", value = 30),
                           numericInput(inputId = "Zigzagthre", label = "zigzagthre to mark type2 and type3", value = 0.2),
                           numericInput(inputId = "Candidate_num", label = "The number of candidate features to output", value = 2)),
                  ##Set xcms MS1 data extraction parameters
                  tabPanel("Parameters of XCMS", 
                           h3(" Parameters of XCMS"),
                           numericInput(inputId = "ppm.5", label = "ppm", value = 10),
                           numericInput(inputId = "noise.5", label = "noise", value = 100),
                           numericInput(inputId = "snthresh.5", label = "snthresh", value = 3),
                           numericInput(inputId = "mzdiff.5", label = "mzdiff", value = -0.001),
                           numericInput(inputId = "peakwidth.5", label = "peakwidth(min)", value = 1),
                           numericInput(inputId = "peakwidth.6", label = "peakwidth(max)", value = 30)),
                  id = "Parameters.of.MS1.5",
                  width = 6,
                  height = 580
                ),
                box( id = "files input",
                     h3("Database input"),
                     status = "danger",
                     fileInput(inputId = "db.path.5", label = "Database", multiple = F,accept = ".xlsx"),
                     h3("EISA MS1 data input"),
                     fileInput(inputId = "msRawData.5", label = "RawData", multiple = F,accept = ".mzXML"),
                     h3("Feature input (Processed by xcms, mzmine etc)"),
                     fileInput(inputId = "feature.xcms.5", label = "Featuretable(1. mz; 2. rt)", multiple = F),
                     #start
                     actionButton("zhu", h3("Start"), width = "100%"),
                     
                     width = 6,
                     height = 550
                )
              ),
              fluidRow(
                tabBox(
                  tabPanel(title = "FeatureTable",
                           dataTableOutput("FeatureTable.5"),
                           downloadButton("down.FeatureTable.5","Download feature table")),
                  tabPanel(title = "AnnotationTable",
                           dataTableOutput("AnnotationTable.5"),
                           downloadButton("down.AnnotationTable.5","Download annotation table")),
                  tabPanel(title = "Filter-AnnotationTable",
                           dataTableOutput("FilterAnnotationTable.5"),
                           downloadButton("down.FilterAnnotationTable.5","Download Filter-AnnotationTable"))
                       ),
                tabBox(
                  tabPanel(title = "Chromatographic peaks and Mass spectrometry",
                           column(5,
                                  textInput(inputId = "ID.6", label = "Enter the ID of the feature to output its chromatogram", value = 2,width = "100%")),
                           column(7,
                                  sliderInput(inputId = "ID.5",label = "Enter the serial number to plot the different chromatographic peaks and Mass spectra of above ID",min = 1,max = 10,value = 1)),
                        
                           actionButton("plot.start.5", h3("Plot"), width = "100%"),  #start plot
                           shinycssloaders::withSpinner(plotOutput("ploteics.5", width = "1400px",height = "900px"))),
                  tabPanel(title = "Extract ion chromatograms(EIC)",
                           textInput(inputId = "ID.7", label = "Enter the ID of the feature to output its chromatogram", value = 2,width = "100%"),
                           actionButton("plot.eic.start", h3("Plot"), width = "100%"),
                           shinycssloaders::withSpinner(plotOutput("ploteics.8", width = "1400px",height = "900px")))
                       ),
               
                ) 
              ),
      ##EISA-EXPOSOME  - end ################################################################
      #tabItem
      tabItem(tabName = "literature", titlePanel("Literature related to EISA technology"),
              div(style="width:fit-content;width:-webkit-fit-content;width:-moz-fit-content;font-size:200%;margin-top:20px",
                  HTML("<b>EISA</b> is a computational tool for exposome targered extraction and annotation based on databases."),
                  HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The literature contains the following parts:"),
                  HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>(1) Quantitative multiple fragment monitoring with enhanced in-source fragmentation/annotation mass spectrometry</b></a>: <a href='https://www.nature.com/articles/s41596-023-00803-0' target='_blank'>View</a>"),
                  HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>(2) Single Quadrupole Multiple Fragment Ion Monitoring Quantitative Mass Spectrometry</b></a>: <a href='https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c01246' target='_blank'>View</a>"),
                  HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>(3) Enhanced in-Source Fragmentation Annotation Enables Novel Data Independent Acquisition and Autonomous METLIN Molecular Identification</b></a>: <a href='https://pubs.acs.org/doi/10.1021/acs.analchem.0c00409' target='_blank'>View</a>"),
                  HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>(4) Proteomics with Enhanced In-Source Fragmentation/Annotation: Applying XCMS-EISA Informatics and Q-MRM High-Sensitivity Quantification</b></a>: <a href='https://pubs.acs.org/doi/10.1021/jasms.1c00188' target='_blank'>View</a> "),
                  #img(src="xcmseisa.jpg", height = 330, width = 800) 
                  )),
      ##########################Bioinformatics tools###############################
      tabItem(tabName ="tools",
              fluidRow(
                box(
                  h3("MS1 data input and pickpeak by XCMS(3.16.1)"),
                  fileInput(inputId = "msRawData.4", label = "RawData", multiple = F,accept = ".mzXML"),
                  numericInput(inputId = "ppm", label = "ppm", value = 10),
                  numericInput(inputId = "noise", label = "noise", value = 100),
                  numericInput(inputId = "snthresh", label = "snthresh", value = 3),
                  numericInput(inputId = "mzdiff", label = "mzdiff", value = -0.001),
                  numericInput(inputId = "peakwidth.1", label = "peakwidth(min)", value = 1),
                  numericInput(inputId = "peakwidth.2", label = "peakwidth(max)", value = 30),
                  actionButton("xcms.start", h3("Start"), width = "100%"),
                  width = 4,
                  height = 710
                ),
                box( title = "FeatureTable",
                     dataTableOutput("FeatureTable.xcms"),
                     downloadButton("down.FeatureTable.xcms","Download FeatureTable.xcms"),
                     width = 8
                ))),
      ##Bioinformatics tools -end
      #construct database
      tabItem(tabName ="ConstructDB",
              fluidRow(
                box(
                  h3("Files input and set parameters"),
                  fileInput(inputId = "information", label = "std-information", multiple = F,accept = ".xlsx"),
                  fileInput(inputId = "ms2RawData", label = "ms2RawData", multiple = F,accept = ".mzXML"),
                  selectInput(inputId = "remove_ions", label = "Whether to perform spectral cleaning", choices = c("Yes","No"),selected = "Yes"),
                  selectInput(inputId = "type", label = "type of database", choices = c("EISA-EXPOSOME","Compound Discovery (undetermined)"),selected = "EISA-EXPOSOME"),
                  numericInput(inputId = "num_ions", label = "Number of retained fragment ions", value = 5),
                  actionButton("construct.start", h3("Start"), width = "100%"),
                  width = 4,
                  height = 600
                ),
                box( title = "Std-database (EISA-EXPOSOME)",
                     dataTableOutput("Std-database"),
                     downloadButton("Std-database.2","Download Std-database"),
                     width = 8
                )))
      ))#end  dashboardBody
)

#server
server <- function(input, output) {
  set.seed(122)
  options(shiny.maxRequestSize=1000*1024^2) #   file size < =1g
  ###################################database###########################################################
  # observeEvent(input$view.DB, {
  #   
  #   database_lists <- input$db  # T3DB NIST MONA GNPS...
  #   
  #   #db <- paste0("/root/srv/shiny-server/EISA-EXPOSOME/www/", x ,".xlsx")
  #   dblist <- openxlsx::read.xlsx("/srv/shiny-server/EISA-EXPOSOME/www/NIST14.xlsx")
  #   output$datalist <- renderDataTable(dblist, options = list(pageLength = 40))
  # })
  
  ###################################Third pick peak####################################################
  observeEvent(input$start.3, {
    withCallingHandlers({
      shinyjs::html("text.3", "")
      #if database input
      if (length(input$db.path.3$datapath)==0){
        dbFile <- "/srv/shiny-server/EISA-EXPOSOME/www/std_200.xlsx"      }
      else {
        dbFile <- input$db.path.3$datapath
      }
      #if rawdata input 
      if (length(input$msRawData.3$datapath)==0){
        msRawData <- "/srv/shiny-server/EISA-EXPOSOME/www/test.mzXML"
        mzXML.files.3 <- read_files(msRawData)}
      else {
        msRawData <- input$msRawData.3 
        ext <- tools::file_ext(msRawData$datapath)
        req(file)
        validate(need(ext == "mzXML", "Please upload a mzXML file"))
        #mzXML.files.3 <- paste0(msRawData.MetEx, "/", grep('.mzXML', dir(msRawData), value = TRUE))
        mzXML.files.3 <- read_files(msRawData$datapath)
        }
      #packageStartupMessage(mzXML.files.3[mzXML.files.3.i])
      featuretable <- targetExtractpeaks(rawData = mzXML.files.3,
                                         dbData = dbFile,
                                         mz_tol = input$MS1deltaMZ.3,
                                         rt_tol = input$MS1deltaTR.3,
                                         rt_window = 60, # use to calculate zigzag and plot eic
                                         m = 20,
                                         zigzagthre = input$zizagThre.3,
                                         intthre = input$intThre.3)
      
      #download  file
      output$down <- downloadHandler(
        filename = function() {
          paste("featureTable-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          openxlsx::write.xlsx(featuretable, file)
        }
      )
      #end-download  file
      shinyjs::show(id = "parameter.hide.button.3")
    },
    #output featuretable
    #Sys.sleep(20),
    output$dynamic <- renderDataTable(featuretable, options = list(pageLength = 10)),
    #plot EIC 
    observeEvent(input$plot.3, {
      #if input id is in ID
      ID <- unique(featuretable$ID)
      if(input$ID %in% ID){
        data <- featuretable[which(featuretable$ID==input$ID),]
        eic <- targetExtractEIC(rawData = mzXML.files.3,mzs = unique(data$PrecursorMZ),rt=unique(data$RT),
                                mz_tol = input$MS1deltaMZ.3,rt_tol = input$MS1deltaTR.3)
        output$plot <- renderPlot(plot(x=eic[,1],y=eic[,2],type=input$datapoints,lwd = 3,col= input$col,
                                       xlab ="rt",ylab="intensity",main= paste0("mz:",unique(data$PrecursorMZ))))   
      }else{
        output$plot <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type=input$datapoints,lwd = 3,col= input$col,
                                       xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
      }
    }),
    #end -ploteic
    message = function(m) {
      shinyjs::html(id = "text.3", html = m$message, add = F)
    })
  })
  ##Third pick peak-end
  
  # ###################################Forth annotation#########################################################
  # observeEvent(input$start, {
  #   withCallingHandlers({
  #     shinyjs::html("text", "")
  #     
  #     if (length(input$db.path$datapath)==0){
  #       dbFile <- "/srv/shiny-server/EISA-EXPOSOME/www/std_200.xlsx"      }
  #     else {
  #       dbFile <- input$db.path$datapath
  #     }
  #    
  #     if (length(input$feature.xcms$datapath)==0){
  #       #
  #       xcms.File <- NULL#
  #     }
  #     else {
  #       xcms.File <- input$feature.xcms$datapath
  #     }
  #     #packageStartupMessage(mzXML.files.3[mzXML.files.3.i])
  #     withProgress(message = 'Annotating !! ',
  #                  detail = 'This may take a while...', value = 0, 
  #                  expr = {
  #                    for (i in 1:4) {
  #                      incProgress(1/4)
  #                      if(i==1){
  #                        #if rawdata input 
  #                        if (length(input$msRawData$datapath)==0){
  #                          msRawData <- "/srv/shiny-server/EISA-EXPOSOME/www/test.mzXML"
  #                          mzXML.files <- read_files(msRawData)
  #                        }
  #                        else {
  #                          msRawData <- input$msRawData 
  #                          ext <- tools::file_ext(msRawData$datapath)
  #                          req(file)
  #                          validate(need(ext == "mzXML", "Please upload a mzXML file"))
  #                          #mzXML.files <- paste0(msRawData.MetEx, "/", grep('.mzXML', dir(msRawData), value = TRUE))
  #                          mzXML.files <- read_files(msRawData$datapath)
  #                        }
  #                      }
  #                      if(i==2){
  #                        featuretable <- targetExtractpeaks(rawData = mzXML.files,
  #                                                           dbData = dbFile,
  #                                                           mz_tol = input$MS1deltaMZ,
  #                                                           rt_tol = input$MS1deltaTR,
  #                                                           rt_window = input$rt_window.2,
  #                                                           m = 20,
  #                                                           zigzagthre = input$zizagThre,
  #                                                           intthre = input$intThre)
  #                      }
  #                      if(i==3){
  #                        Output_table <- matchMS2score(rawData = mzXML.files,
  #                                                      target_featuretable = featuretable,
  #                                                      dbData = dbFile,
  #                                                      mz_tol = input$MS1deltaMZ, # extract eics
  #                                                      rt_window = input$rt_window.2, 
  #                                                      ms2_tol = input$MS2deltaMZ,
  #                                                      cores = 1)
  #                      }
  #                    }
  #                  })
  #     #MARK LEVEL
  #     if(is.null(xcms.File)==FALSE){
  #       Output_table <- markAnnotationLevel(Output_table = Output_table,
  #                                           xcms.featuretable = xcms.File)
  #     }
  #     ##MRAK LEVEL END
  #     #download file
  #     output$down.annotation <- downloadHandler(
  #       filename = function() {
  #         paste("annotationTable-", Sys.Date(), ".xlsx", sep="")
  #       },
  #       content = function(file) {
  #         openxlsx::write.xlsx(Output_table, file)
  #       }
  #     )
  #     #end-download file
  #     shinyjs::show(id = "parameter.hide.button")
  #   },
  #   #Output featuretable
  #   #Sys.sleep(20),
  #   output$AnnotationTable <- renderDataTable(Output_table, options = list(pageLength = 20)),
  #   #plotEIC 
  #   observeEvent(input$plot.start, {
  #     
  #     ID <- unique(Output_table$DB.ID)
  #     if(input$ID.1 %in% ID){
  #       data <- Output_table[which(Output_table$DB.ID==input$ID.1),]
  #       database_data <- openxlsx::read.xlsx(dbFile)
  #       database_data.1 <- database_data[which(database_data$ID==input$ID.1),]
  #       mz_lib <- database_data.1$ProductMZ
  #       int_lib <- as.numeric(database_data.1$Intensity)*-1
  #       lib_msms <- data.frame(mz = mz_lib,int = int_lib,msms_type = "lib")
  #       f_rt <- unique(data$apex_rt)
  #       #choose figure
  #       if(input$ID.3 %in% 1:length(f_rt)){
  #         #
  #         data2 <- data[which(data$apex_rt==f_rt[input$ID.3]),]
  #         int_exp <- data2$intensity_exp[-1]/max(data2$intensity_exp[-1])*100
  #         exp_msms <- data.frame(mz = data2$mz_exp[-1],int = int_exp,msms_type = "exp")
  #         msms <- rbind(exp_msms,lib_msms)
  #         eics <- targetExtractEIC(rawData = mzXML.files,mzs = data2$mz_exp,rt=unique(data2$apex_rt)/60,
  #                                  mz_tol = input$MS1deltaMZ,rt_tol = input$rt_window.2/2)%>% as.data.frame()
  #         #precursor ion eic
  #         p_eic <- eics[,c(1,2)]
  #         p_eic <- reshape2::melt(p_eic,id="rt")  #data.frame
  #         colnames(p_eic) <- c("rt","mz","intens")
  #         p0 <- p_eic %>% 
  #           ggplot(aes(x = rt,
  #                      y = intens)) +
  #           geom_line(aes(color = mz,
  #                         group = mz)) +
  #           geom_point()+
  #           guides(color = "none") +
  #           labs(x = "Retention time [s]",
  #                y = "Intensity",
  #                title = paste0("Precursor mz: ", data2$mz_exp[1],"  zigzag(P): ",unique(data2$zigzag))) +
  #           #ggforce::facet_zoom(y = intens < median(result$intensity_exp),zoom.size = 1) +
  #           theme_minimal()
  #         
  #         #p+f
  #         eics <- reshape2::melt(eics,id="rt")  #data.frame
  #         colnames(eics) <- c("rt","mz","intens")
  #         p1 <- eics %>% 
  #           ggplot(aes(x = rt,
  #                      y = intens)) +
  #           geom_line(aes(color = mz,
  #                         group = mz)) +
  #           
  #           guides(color = "none") +
  #           labs(x = "Retention time [s]",
  #                y = "Intensity",
  #                title = paste0("precursor and fragments eics")) +
  #           #ggforce::facet_zoom(y = intens < median(result$intensity_exp),zoom.size = 1) +
  #           theme_minimal()
  #         
  #         p2 <- ggplot(data=msms,mapping=aes(x=mz,y=int))+
  #           ylim(-100,100)+
  #           geom_hline(yintercept = 0, color = "black", linewidth = 1) + # add y=0
  #           geom_point(aes(color = msms_type), size = 1) +         
  #           geom_bar(aes(fill = msms_type), stat="identity", width=0.2)+
  #           labs(x = "mz",y = "Intensity",title = "Mass spectra of matched fragment ions")+
  #           theme_bw(base_family = "Times") +
  #           geom_text(aes(label = round(mz,4), vjust = -0.5, hjust = 0.5), show.legend = TRUE)+
  #           annotate("text", label = paste0("  MFR : ",unique(data2$MFR),"  SSM : ",unique(data2$SSM)),x=min(mz_lib)+20,y = 90, size = 8, colour = "red")
  #         output$ploteics <- renderPlot((p0+p1)/p2)          
  #         #Sys.sleep(5)
  #       }else{
  #         output$ploteics <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type="l",lwd = 3,col= "red",
  #                                            xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
  #       }
  #     }else{
  #       output$ploteics <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type="l",lwd = 3,col= "red",
  #                                          xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
  #     }        
  #   }),
  #   #end -ploteic
  #   message = function(m) {
  #     shinyjs::html(id = "text.3", html = m$message, add = F)
  #   })
  # })
  # ##Forth Annotation-end
  
  ###################################Filter-annotation#########################################################
  observeEvent(input$zhu, {
      # Determine if there is a database input
      if (length(input$db.path.5$datapath)==0){
        dbFile <- "/srv/shiny-server/EISA-EXPOSOME/www/std_200.xlsx"      }
      else {
        dbFile <- input$db.path.5$datapath
      }
      # Determine if there is a featuretable input
      if (length(input$feature.xcms.5$datapath)==0){
        #
        xcms.File <- pickpeak.XCMS(rawData = input$msRawData.5$datapath,
                                   ppm = input$ppm.5,
                                   mzdiff = input$mzdiff.5 ,
                                   snthresh = input$snthresh.5,
                                   noise = input$noise.5,
                                   peakwidth = c(input$peakwidth.5,input$peakwidth.6))
      }
      else {
        xcms.File <- input$feature.xcms.5$datapath
      }
      
      withProgress(message = 'Annotating !! ',
                   detail = 'This may take a while...', value = 0, 
                   expr = {
                     for (i in 1:5) {
                       incProgress(1/5)
                       if(i==1){
                         #if rawdata input 
                         if (length(input$msRawData.5$datapath)==0){
                           msRawData <- "/srv/shiny-server/EISA-EXPOSOME/www/test.mzXML"
                           mzXML.files <- read_files(msRawData)
                         }
                         else {
                           msRawData <- input$msRawData.5 
                           ext <- tools::file_ext(msRawData$datapath)
                           req(file)
                           validate(need(ext == "mzXML", "Please upload a mzXML file"))
                           #mzXML.files <- paste0(msRawData.MetEx, "/", grep('.mzXML', dir(msRawData), value = TRUE))
                           mzXML.files <- read_files(msRawData$datapath)
                         }
                       }
                       if(i==2){
                         featuretable <- targetExtractpeaks(rawData = mzXML.files,
                                                            dbData = dbFile,
                                                            mz_tol = input$MS1deltaMZ.5,
                                                            rt_tol = input$pickpeak_window, #no rt,set larger than 1000
                                                            rt_window = input$rt_window.3,
                                                            start_elution = input$Lower_RT_Limit,
                                                            end_elution = input$Upper_RT_Limit,
                                                            minpoints = input$minpoints,
                                                            m = input$steps,
                                                            zigzagthre = input$zizagThre.5,
                                                            intthre = input$intThre.5,
                                                            cores = 4)
                         
                         print("finish extraing features ")  
                       }
                       if(i==3){
                         Output_table <- matchMS2score(rawData = mzXML.files,
                                                       target_featuretable = featuretable,
                                                       dbData = dbFile,
                                                       mz_tol = input$MS1deltaMZ.5, # extract eics
                                                       rt_window = input$rt_window.3,##use to calculate zigzag and peakcor
                                                                         
                                                       ms2_tol = input$MS2deltaMZ.5,
                                                       cores = 4)
                         print("finish annotating") 
                       }
                       if(i==4){
                         print("start filtering")
                         filter_result <- selectfeature(Output_table = Output_table,
                                                        featuretable_xcms = xcms.File,
                                                        SSM_weight = input$SSM_weigth,
                                                        MFR_weight = input$MFR_weight,
                                                        zigzagthre_type = input$Zigzagthre,
                                                        candidate_num = input$Candidate_num,
                                                        mztol_xcms = input$Delta_mz_xcms,
                                                        rttol_xcms = input$Delta_rt_xcms)
                         print("finish filtering")
                       }
                       
                     }
                   })
      
      ##
      # output featuretable
      #Sys.sleep(20),
      output$FeatureTable.5 <- renderDataTable(featuretable, options = list(pageLength = 20))
      
      output$AnnotationTable.5 <- renderDataTable(Output_table, options = list(pageLength = 20))
      
      output$FilterAnnotationTable.5 <- renderDataTable(filter_result, options = list(pageLength = 20))
      # download files
      # feature
      output$down.FeatureTable.5 <- downloadHandler(
        filename = function() {
          paste("FeatureTable-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          openxlsx::write.xlsx(featuretable, file)
        }
      )
      #annotation
      output$down.AnnotationTable.5 <- downloadHandler(
        filename = function() {
          paste("annotationTable-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          openxlsx::write.xlsx(Output_table, file)
        }
      )
      #filter annotation
      output$down.FilterAnnotationTable.5 <- downloadHandler(
        filename = function() {
          paste("filterannotationTable-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          openxlsx::write.xlsx(filter_result, file)
        }
      )
      #end-downloadfile
    # plot EIC peak
      Output_table <- filter_result 
      #Output_table <- Output_table[-which(Output_table==" "),]
      #Output_table$DB.ID <- as.numeric(Output_table$DB.ID)
      ID <- unique(Output_table$DB.ID)
    #
    observeEvent(input$plot.start.5, {
      # Read if the ID is in the table before allowing it to run.
      if(input$ID.6 %in% ID){
        data <- Output_table[which(Output_table$DB.ID==input$ID.6),]
        database_data <- openxlsx::read.xlsx(dbFile)
        database_data.1 <- database_data[which(database_data$ID==input$ID.6),]
        # MS/MS in library
        mz_lib <- database_data.1$ProductMZ
        int_lib <- as.numeric(database_data.1$Intensity)*-1
        lib_msms <- data.frame(mz = mz_lib,int = int_lib,msms_type = "lib")
        # Use retention time to distinguish between different peaks of the same m/z
        f_rt <- unique(data$apex_rt)
        # Choose which chart to output
        if(input$ID.5 %in% 1:length(f_rt)){
          #
          data2 <- data[which(data$apex_rt==f_rt[input$ID.5]),]
          int_exp <- data2$intensity_exp[-1]/max(data2$intensity_exp[-1])*100
          exp_msms <- data.frame(mz = data2$mz_exp[-1],int = int_exp,msms_type = "exp")
          msms <- rbind(exp_msms,lib_msms)
          eics <- targetExtractEIC(rawData = mzXML.files,
                                   mzs = data2$mz_exp,
                                   rt=unique(data2$apex_rt)/60,
                                   mz_tol = input$MS1deltaMZ.5,
                                   rt_tol = input$rt_window.3/2)%>% as.data.frame()
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
                 title = paste0("Precursor mz: ", data2$mz_exp[1],"  zigzag(P): ",unique(data2$zigzag))) +
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
          
          p2 <- ggplot(data=msms,mapping=aes(x=mz,y=int))+
            ylim(-100,100)+
            geom_hline(yintercept = 0, color = "black", linewidth = 1) + # add y=0
            geom_point(aes(color = msms_type), size = 1) +        
            geom_bar(aes(fill = msms_type), stat="identity", width=0.2)+
            labs(x = "mz",y = "Intensity",title = "Mass spectra of matched fragment ions")+
            theme_bw(base_family = "Times") +
            geom_text(aes(label = round(mz,4), vjust = -0.5, hjust = 0.5), show.legend = TRUE)+
            annotate("text", label = paste0("  MFR : ",unique(data2$MFR),"  SSM : ",unique(data2$SSM),
                                            "  W_score : ",unique(data2$weigthed_score), "  type : ",unique(data2$type)),
                     x=min(mz_lib)+20,y = 90, size = 4, colour = "red")
          output$ploteics.5 <- renderPlot((p0+p1)/p2)          
          #Sys.sleep(5)
        }else{
          output$ploteics.5 <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type="l",lwd = 3,col= "red",
                                             xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
        }
      }else{
        output$ploteics.5 <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type="l",lwd = 3,col= "red",
                                           xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
      }        
    })
    #plot eic
    observeEvent(input$plot.eic.start, {
      if(input$ID.7 %in% ID){
        data <- Output_table[which(Output_table$DB.ID==input$ID.7),]
        mzmed <- unique(data$mz_exp[1])
        mz_range <- matrix(c(mzmed - input$MS1deltaMZ.5, mzmed + input$MS1deltaMZ.5), ncol = 2) 
        precursorEIC <- get_eic(object = mzXML.files,
                                mz = mz_range,
                                rt = NULL,
                                aggregationFun = "max") %>% as.data.frame()
        colnames(precursorEIC) <- c("rt", sprintf("%.4f", mzmed))
        eics <- precursorEIC
        p_eic <- eics[,c(1,2)]
        p_eic <- reshape2::melt(p_eic,id="rt")  #data.frame
        colnames(p_eic) <- c("rt","mz","intens")
        p <- p_eic %>% 
                 ggplot(aes(x = rt,
                            y = intens)) +
                 geom_line(aes(color = mz,
                               group = mz)) +
                 #geom_point()+
                 annotate("rect", xmin = unique(data$apex_rt)-20, 
                          xmax = unique(data$apex_rt)+20,
                          ymin = 0, ymax = max(data$intensity_exp),
                          fill = "blue", 
                          alpha = 0.2)+
                 #guides(color = "none") +
                 labs(x = "Retention time [s]",
                      y = "Intensity",
                      title = i) +
                 geom_vline(xintercept = c(input$Lower_RT_Limit,input$Upper_RT_Limit), color = "black", linewidth = 1)+
                 theme(axis.text.y=element_text(vjust=1,size=20,face = "bold"),
                       axis.text.x=element_text(vjust=1,size=20,face = "bold"))
                 #ggforce::facet_zoom(y = intens < median(result$intensity_exp),zoom.size = 1) +
                 #theme_minimal()
        output$ploteics.8<- renderPlot(p)
      }else{
        output$ploteics.8 <- renderPlot(plot(x=rep(0,60),y=rep(0,60),type="l",lwd = 3,col= "red",
                                             xlab ="rt",ylab="intensity",main= paste0("no this feature !!!"))) 
      }      
    })
    #
    
  })
  ##Filter-Annotation-end
  
  ###################################Bioinformatics tools#########################################################
  observeEvent(input$xcms.start, {
      featuretable.xcms <- pickpeak.XCMS(rawData = input$msRawData.4$datapath,
                                         ppm = input$ppm,
                                         mzdiff = input$mzdiff,
                                         snthresh = input$snthresh,
                                         noise = input$noise,
                                         peakwidth=c(input$peakwidth.1,input$peakwidth.2))
      # Output featuretable
      # Sys.sleep(20),
      output$FeatureTable.xcms <- renderDataTable(featuretable.xcms, options = list(pageLength = 15))
      # download
      output$down.FeatureTable.xcms <- downloadHandler(
        filename = function() {
          paste("FeatureTablexcms-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          openxlsx::write.xlsx(featuretable.xcms, file)
        }
      )
  })
  ##############################Bioinformatics tools-end###########################################################
  
  ###################################construct database#########################################################
  observeEvent(input$construct.start, {
    #EISA-EXPOSOME OR CD
    std-database <- construct_db(MS2.files = input$ms2RawData$datapath,
                                 featuretale_path = input$information,
                                 clean_spe = input$remove_ions,
                                 keep_ion_num = input$num_ions,
                                 type = input$type)
    # Output featuretable
    # Sys.sleep(20),
    output$Std-database <- renderDataTable(std-database, options = list(pageLength = 50))
    # download
    output$Std-database.2 <- downloadHandler(
      filename = function() {
        paste("Std-database-", Sys.Date(), ".xlsx", sep="")
      },
      content = function(file) {
        openxlsx::write.xlsx(std-database, file)
      }
    )
  })
  ##############################construct database-end###########################################################
}

shinyApp(ui, server)
