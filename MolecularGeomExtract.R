#################################################################################################################
###                                                                                                           ###
###                                                                                                           ###
###   Feature extraction from gaussian input file                                                             ###
###   1- Aggregate percentile-based features                                                                  ###
###   2- Coulomb matrix eigen values                                                                          ###
###   3- Electron density surrogate eigen values                                                              ###
###                                                                                                           ###
###                                                                                                           ###
#################################################################################################################

closeAllConnections()
rm(list = ls())
#set the data reproducible
set.seed(2)
setwd("C:\\Working Directory\\")
dir = "C:\\Working Directory\\"

filenames = list.files(path=dir, pattern="\\.com");
filenames
dftemp <- read.csv("ctemp.csv", header=TRUE)
#instead of atemp use btemp and put everything in a single row all calculations and make totalvarfile
# all this will be done in the same for loop below
#dftemp <- NULLFileName_Energy_Filled_complete
Gaussian_energy_results <- read.csv("FileName_Energy_Filled_complete.csv" , header=TRUE, na.strings ="", stringsAsFactors= F)
atomicNumdf <- read.csv("catomnumber.csv", na.strings=c("", "NA"), stringsAsFactors=FALSE, header=TRUE)

for (i in 1:length(filenames)){
  #dftemp = dftemp[-1,]
  dftemp[nrow(dftemp)+1,] <- NA
  dftemp$filenames[i]<- sub(".com","", filenames[i])
  dftemp$filenamesfull[i] <- filenames[i]
  Gaussian_input_data <- read.table(file.path(dir,filenames[i]), colClasses = rep("character", 4), sep = "" , na.strings ="", stringsAsFactors= F, fill = TRUE, comment.char = "")
  dftemp$Charge[i] <- Gaussian_input_data[5,1]
  if (identical(Gaussian_input_data[which(Gaussian_input_data=="b3lyp/aug-cc-pvtz", arr.ind=TRUE)], character(0))==FALSE){
    dftemp$Theory.Basis.set[i] <- "b3lyp/aug-cc-pvtz"
  }else if (identical(Gaussian_input_data[which(Gaussian_input_data=="wb97xd/dgdzvp", arr.ind=TRUE)], character(0))==FALSE){
    dftemp$Theory.Basis.set[i] <-"wb97xd/dgdzvp"
  }else if (identical(Gaussian_input_data[which(Gaussian_input_data=="B1B95/3-21G", arr.ind=TRUE)], character(0))==FALSE){
    dftemp$Theory.Basis.set[i] <- "B1B95/3-21G"
  }else if (identical(Gaussian_input_data[which(Gaussian_input_data=="B3LYP/6-31+G**", arr.ind=TRUE)], character(0))==FALSE){
    dftemp$Theory.Basis.set[i] <- "B3LYP/6-31+G**"
  }
  #which(Gaussian_energy_results==dftemp["filenamesfull"], arr.ind = TRUE)
  if (identical(which(Gaussian_energy_results$filenames==paste0(dftemp$filenames[i],".out"), arr.ind = TRUE),integer(0))){
    dftemp$Energy[i] <- NA
  }else{
  dftemp$Energy[i] <- Gaussian_energy_results$energy[which(Gaussian_energy_results$filenames==paste0(dftemp$filenames[i],".out"), arr.ind = TRUE)]
  }
  atomspresent <- Gaussian_input_data[6:nrow(Gaussian_input_data),1]
  X<-as.numeric(Gaussian_input_data[6:nrow(Gaussian_input_data),2])
  Y<-as.numeric(Gaussian_input_data[6:nrow(Gaussian_input_data),3])
  Z<-as.numeric(Gaussian_input_data[6:nrow(Gaussian_input_data),4])
  numberofprotons = 0
  numberofprotonsarray <- matrix(0, nrow=1, ncol=length(atomspresent))
  for (iter2 in 1:length(atomspresent)){ 
    numberofprotons = numberofprotons + atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter2])]
    numberofprotonsarray[iter2] = atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter2])]
  }
  
  Charge = as.numeric(dftemp$Charge[i])
  dftemp$N.electron[i] = numberofprotons - Charge
  dftemp$M.nuclei[i] <- length(atomspresent)
  dftemp$M.protons[i] <- numberofprotons
  centroidX=mean(X)
  centroidY=mean(Y)
  centroidZ=mean(Z)
  centroidXW=weighted.mean(X, numberofprotonsarray)
  centroidYW=weighted.mean(Y, numberofprotonsarray)
  centroidZW=weighted.mean(Z, numberofprotonsarray)
  
  XYZsum=0
  distantMatrix <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  WeighteddistantMatrix <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  centroiddist <- matrix(0, nrow=1, ncol=length(atomspresent))
  centroidWdist <- matrix(0, nrow=1, ncol=length(atomspresent))
  Elec_Den_Surg1_Mx <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  errorFlagdistance<-0
  for (iter3 in 1:length(atomspresent)){ 
    for (iter4 in 1:length(atomspresent)){ 
      if (iter3<=iter4) {next} 
      Xsumtemp=(X[iter3]-X[iter4])^2
      Ysumtemp=(Y[iter3]-Y[iter4])^2
      Zsumtemp=(Z[iter3]-Z[iter4])^2
      distantMatrix[iter3,iter4] =sqrt(Xsumtemp+Ysumtemp+Zsumtemp)
      if (distantMatrix[iter3,iter4] <= 0.001){
        errorFlagdistance<-1
        errorFlagdistanceindex<- errorFlagdistanceindex + 1
        removerows[errorFlagdistanceindex] <- a
        break
      }
      WeighteddistantMatrix[iter3,iter4] = distantMatrix[iter3,iter4] * (atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter3])]+atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter4])])
      Elec_Den_Surg1_Mx[iter3,iter4] = (atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter3])]+atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter4])])/distantMatrix[iter3,iter4]
    }
    if(errorFlagdistance==1){break}
    centroiddist[iter3]= sqrt(((X[iter3]-centroidX)^2)+((Y[iter3]-centroidY)^2)+((Z[iter3]-centroidZ)^2))
    centroidWdist[iter3]= sqrt(((X[iter3]-centroidXW)^2)+((Y[iter3]-centroidYW)^2)+((Z[iter3]-centroidZW)^2))
  }
  if(errorFlagdistance==1){next}
  sorteddisMx <-NULL
  sorteddisMx=sort(distantMatrix)
  sortedWtdisMx=sort(WeighteddistantMatrix)
  sortedElecDenSurMx=sort(Elec_Den_Surg1_Mx)
  sorteddisMx <- sorteddisMx[sorteddisMx!=0]
  sortedWtdisMx <- sortedWtdisMx[sortedWtdisMx!=0]
  sortedElecDenSurMx <- sortedElecDenSurMx[sortedElecDenSurMx!=0]
  
  HistFitTrial <- function() normalmixEM(sorteddisMx)
  HistFitError <- function() "HistFitError"
  fit <- tryCatch(
    {
      HistFitTrial()
    },
    error = function(e){
      HistFitError()
    })
  
  if(fit!="HistFitError"){
    fit = normalmixEM(sorteddisMx)
    dftemp$Nuclei.Gaussian.Distance.Average[i] = fit$mu[1]
    dftemp$Nuclei.Gaussian.Standard.Deviation[i] = fit$sigma[1]
  }else{
    dftemp$Nuclei.Gaussian.Distance.Average[i] = mean(sorteddisMx)
    dftemp$Nuclei.Gaussian.Standard.Deviation[i] = sd(sorteddisMx)
  }
  #if(is.na(df_geom$Nuclei.Gaussian.Standard.Deviation)){df_geom$Nuclei.Gaussian.Standard.Deviation=0}
  dftemp$Nuclei.Distance.Average[i] = mean(sorteddisMx)
  dftemp$Nuclei.Distance.Average.Proton.Weighted[i] = mean(sortedWtdisMx)
  normalizedsortedWtdisMx= sortedWtdisMx/max(sortedWtdisMx)
  histnormalizedsortedWtdisMx <- hist(normalizedsortedWtdisMx, breaks=seq(0,max(normalizedsortedWtdisMx),l=11), plot = FALSE)
  dftemp$Nuclei.Distance.Proton.Weighted.Max[i] <- max(sortedWtdisMx)
  dftemp$Nuclei.Distance.Proton.Weighted.Min[i] <- min(sortedWtdisMx)

  dftemp$Centroid.average.dist.from.atoms[i] <- mean(centroiddist)
  dftemp$CentroidW.average.dist.from.atoms[i] <- mean(centroidWdist)
  
  dftemp[i,21:31] <- quantile(numberofprotonsarray, prob = seq(0,1,0.1))
  dftemp[i,32:41] <- histnormalizedsortedWtdisMx$counts
  dftemp[i,42:52] <- quantile(sortedWtdisMx, prob = seq(0,1,0.1))
  dftemp[i,53:63] <- quantile(centroiddist, prob = seq(0,1,0.1))
  dftemp[i,64:74] <- quantile(centroidWdist, prob = seq(0,1,0.1))
  dftemp[i,75:85] <- quantile(sortedElecDenSurMx, prob= seq(0,1,0.1))
  dftemp[i,86:96] <- quantile(atomicNumdf$atomic.radius[(atomicNumdf$Atom %in% atomspresent)], prob=seq(0,1,0.1))
  X_highest_p=X[numberofprotonsarray == max(numberofprotonsarray)]
  Y_highest_p=Y[numberofprotonsarray == max(numberofprotonsarray)]
  Z_highest_p=Z[numberofprotonsarray == max(numberofprotonsarray)]
  X_highest_p2=X[numberofprotonsarray==max(numberofprotonsarray[numberofprotonsarray!=max(numberofprotonsarray)])]
  Y_highest_p2=Y[numberofprotonsarray==max(numberofprotonsarray[numberofprotonsarray!=max(numberofprotonsarray)])]
  Z_highest_p2=Z[numberofprotonsarray==max(numberofprotonsarray[numberofprotonsarray!=max(numberofprotonsarray)])]
  X_lowest_p=X[numberofprotonsarray == min(numberofprotonsarray)]
  Y_lowest_p=Y[numberofprotonsarray == min(numberofprotonsarray)]
  Z_lowest_p=Z[numberofprotonsarray == min(numberofprotonsarray)]
  X_lowest_p2=X[numberofprotonsarray==min(numberofprotonsarray[numberofprotonsarray!=min(numberofprotonsarray)])]
  Y_lowest_p2=Y[numberofprotonsarray==min(numberofprotonsarray[numberofprotonsarray!=min(numberofprotonsarray)])]
  Z_lowest_p2=Z[numberofprotonsarray==min(numberofprotonsarray[numberofprotonsarray!=min(numberofprotonsarray)])]
  #min(XX[XX!=min(XX)])
  centroiddistH <- matrix(0, nrow=1, ncol=length(X_highest_p))
  centroiddistL <- matrix(0, nrow=1, ncol=length(X_highest_p))
  centroiddistH2 <- matrix(data=NA, nrow=1, ncol=length(X_highest_p))
  centroiddistL2 <- matrix(data=NA, nrow=1, ncol=length(X_highest_p))
  for (iter5 in 1:length(X_highest_p)){centroiddistH[iter5]=(centroidX-X_highest_p)^2+(centroidY-Y_highest_p)^2+(centroidZ-Z_highest_p)^2}
  for (iter5 in 1:length(X_lowest_p)){centroiddistL[iter5]=(centroidX-X_lowest_p)^2+(centroidY-Y_lowest_p)^2+(centroidZ-Z_lowest_p)^2}
  if (!length(X_lowest_p2)==0){
    for (iter5 in 1:length(X_highest_p2)){centroiddistH2[iter5]=(centroidX-X_highest_p2)^2+(centroidY-Y_highest_p2)^2+(centroidZ-Z_highest_p2)^2}
    for (iter5 in 1:length(X_lowest_p2)){centroiddistL2[iter5]=(centroidX-X_lowest_p2)^2+(centroidY-Y_lowest_p2)^2+(centroidZ-Z_lowest_p2)^2}
  }
  dftemp$Centroid.dist.from.High.e[i] <- mean(centroiddistH)
  dftemp$Centroid.dist.from.low.e[i] <- mean(centroiddistL)
  dftemp$Centroid.dist.from.High.e2[i] <- mean(centroiddistH2)
  dftemp$Centroid.dist.from.low.e2[i] <- mean(centroiddistL2)
  
  dftemp$average.atomic.radius[i] <- mean(atomicNumdf$atomic.radius[(atomicNumdf$Atom %in% atomspresent)])
  dftemp$atomic.rad.H[i] <- max(atomicNumdf$atomic.radius[(atomicNumdf$Atom %in% atomspresent)])
  dftemp$atomic.rad.L[i] <- min(atomicNumdf$atomic.radius[(atomicNumdf$Atom %in% atomspresent)])
  
  
  #constructing diagonal columb matrix
  DistantMatrix <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  ColumbMatrix <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  Elec_Den_Surg2_Mx <- matrix(0, nrow=length(atomspresent), ncol=length(atomspresent))
  
  n_atoms = length(atomspresent)
  n_edges = n_atoms*(n_atoms-1)/2
  for (iter6 in 1:length(atomspresent)){ 
    for (iter7 in 1:length(atomspresent)){ 
      if (iter6<iter7) {next} 
      if (iter6==iter7){
        ColumbMatrix[iter6,iter7] = 0.5*(atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter6])])^2.4
        Elec_Den_Surg2_Mx[iter6,iter7] = atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter6])]-Charge/n_atoms
        next
      }
      Xsumtemp=(X[iter6]-X[iter7])^2
      Ysumtemp=(Y[iter6]-Y[iter7])^2
      Zsumtemp=(Z[iter6]-Z[iter7])^2
      DistantMatrix[iter6,iter7] =sqrt(Xsumtemp+Ysumtemp+Zsumtemp)
      ColumbMatrix[iter6,iter7] = distantMatrix[iter6,iter7] * atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter6])]+atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter7])]
      ColumbMatrix[iter7,iter6] = ColumbMatrix[iter6,iter7]
      Elec_Den_Surg2_Mx[iter6,iter7] = (atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter6])]/(n_atoms-1)+atomicNumdf$atomic.number[(atomicNumdf$Atom %in% atomspresent[iter7])]/(n_atoms-1)-Charge/n_edges)/distantMatrix[iter6,iter7]
      Elec_Den_Surg2_Mx[iter7,iter6] = Elec_Den_Surg2_Mx[iter6,iter7]
    }
  }
  
  #extracting eigenvalues
  ev_ColumbMatrix <- rep(0, 100)
  ev_Elec_Den_Surg2_Mx <- rep(0, 100)
  ev_ColumbMatrix [1:n_atoms]<- eigen(ColumbMatrix)$values
  ev_Elec_Den_Surg2_Mx [1:n_atoms]<- eigen(Elec_Den_Surg2_Mx)$values
  #adding misc vars (eg. number of atoms, charge)
  dftemp[i,97:196] <- ev_ColumbMatrix
  dftemp[i,197:296] <- ev_Elec_Den_Surg2_Mx
  
  #df_geom_total <- rbind(df_geom_total, df_geom)
  
  
} 


mypath <- file.path(dir ,paste("Unoptimizedstructures2_set3.csv", sep = ""))
write.csv(dftemp,file=mypath ,row.names=FALSE, na="")



