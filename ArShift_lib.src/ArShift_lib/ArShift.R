################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
# This is the interfacing script to the ArShift program, which predicts        #
# 1H chemical shifts of the Tyr and Phe aromatic rings.                        #
################################################################################

### ENTER COMMAND (cmd):
 #title        = "first_test" 
 #pdbname      = "a7_b03.pdb"
 #experdata    = NULL #a7_out.shift
 #examine      = 0 #0 for all, 1,2,3 etc... for specific ones
 #pdbtype      = "almost" # "complex", "almost", "medium", "simple"
 #seqshift     = 0
 #rereference  = FALSE #TRUE/FALSE
 #outorder     = "bytype" #"bytype", "byseq", "byshift"
 #outputname   = "out.txt"
 #bias         = 1        # positive number to bias color palette
### EXAMPLE:
#cmd <- 'pdbname="example.pdb", title="Ubiquitin_1", experdata="example.exp", bias=1, examine=1, pdbtype="almost", seqshift=0, rereference=FALSE, outorder="bytype", outputname="out.txt"'
  
# # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS
# # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS

doval  = TRUE
x11lib = TRUE

if(IsThisWebServer==TRUE) {
  if(title=="")    {title="Not specified."}
  if(pdbname=="")  {stop("The PDB file is not uploaded")}
  outputname  =  "out.txt"
  prefix.ws   =  "../../../"
  load("../../../ArShift_lib/ArShift.ini")
} else {
  source("command.cmd")
  load("ArShift_lib/ArShift.ini")
  # generates and reads the commands based on created cmd line
  eval(parse(text=paste("ArShift.cmd(",cmd,")")))
  load("ArShift.cmd")
  file.remove("ArShift.cmd")
  if(pdbname=="")  {stop("The PDB file is not uploaded")}
  # 
  prefix.ws <- NULL 
}
if(experdata==""){experdata=NULL}
# # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS


# START

OUT <- "NOTE: * * * ArShift * * * ArShift * * *  "
OUT[2] <- paste("NOTE: Title -- ", title, sep="")
write(OUT, file="process_info.txt")

# # # # # # # # ####### PDB file examination: ####### # # # # # # # # # #

  isNMR <- isNMRpdb(pdbname)
  OUT <- c(OUT, paste("NOTE: Hydrogen atom presence - ",isNMR$isH,sep=""))
         print(OUT[length(OUT)]) # PRINTOUT
  OUT <- c(OUT, paste("NOTE: Number of structures in the model - ",isNMR$Nmodels,sep=""))
         print(OUT[length(OUT)]) # PRINTOUT
  models2analyse <- 1:isNMR$Nmodels

  if(isNMR$isH!=TRUE){
    OUT <- c(OUT,"ERROR: Please add hydrogens to your molecule!")
           stop(OUT[length(OUT)]) # PRINTOUT
  }
  if(isNMR$Nmodels > 1 & examine==0){
    OUT <- c(OUT, paste("NOTE: All the ",isNMR$Nmodels,
             " models will be analyzed with the shifts averaged!",sep=""))
             print(OUT[length(OUT)]) # PRINTOUT
    models2analyse <- 1:isNMR$Nmodels
  }

  if(isNMR$Nmodels > 1 & examine > 0 & examine <= isNMR$Nmodels){
    OUT <- c(OUT,paste("NOTE: Only the structure ",examine," will be analyzed.",sep=""))
           print(OUT[length(OUT)]) # PRINTOUT
    isNMR$Nmodels <- 1
    models2analyse <- examine
  }
  
  write(OUT[3:length(OUT)], file="process_info.txt", append=TRUE)

# # # # # # # # ######## End of PDB file examination: ####### # # # # # # 



# Loading the parameters:
OUT <- c(OUT, "NOTE: Loading the parameters..."); print(OUT[length(OUT)])
write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
load(paste(prefix.ws,PATH.PAR,"ArShift.par",sep=""))
OUT <- c(OUT, "NOTE: Parameters are loaded."); print(OUT[length(OUT)])
write(OUT[length(OUT)], file="process_info.txt", append=TRUE)


#########################################################################
RESID <- SEQ <- NUCL <- CHAINID <- INTERCHAIN.firstM <- NULL
ChSh.pred <- SE_fit <- SD_train <- list(NULL)
ChSh.pred[[maxnumofentry+1]] <- 
SE_fit[[maxnumofentry+1]]    <- 
SD_train[[maxnumofentry+1]]  <- NA
#########################################################################

count.tried.model <- 1
for(model in models2analyse){

  # Generating the file for a particular model (model.pdb)
  file.copy(pdbname, "model.pdb", overwrite=TRUE)

  #--
  if(model==1){                                                                          ###15OCT12###
    resName.ForQuickCheck <- pdbparser.simple("model.pdb")$resName                       ###15OCT12###
    wrong.naming.ind <- c(which(resName.ForQuickCheck=="LYP"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="CYN"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="CYS2"),                          ###15OCT12###
                          which(resName.ForQuickCheck=="HID"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="HIE"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="HISE"),                          ###15OCT12###
                          which(resName.ForQuickCheck=="HISD"))                          ###15OCT12###
    if(length(wrong.naming.ind)!=0){                                                     ###15OCT12###
      OUT <- c(OUT, "ERROR: The pdb file uses force field related naming for residues,") ###15OCT12###
      OUT <- c(OUT, "ERROR: such as LYP, CYN, CYS2, HID, HIE, HISE, HIDE etc. Please")   ###15OCT12###
      OUT <- c(OUT, "ERROR: change them into normal amino acid naming convention.")      ###15OCT12###
      OUT <- c(OUT, "ERROR: Also, make sure to have delta-protonated histidine in your") ###15OCT12###
      OUT <- c(OUT, "ERROR: structure (still named as HIS) for an extra precision!")     ###15OCT12###
      write(OUT[(length(OUT)-4):length(OUT)], file="process_info.txt", append=TRUE)      ###15OCT12###
      stop(OUT[(length(OUT)-4):length(OUT)])                                             ###15OCT12###
    }                                                                                    ###15OCT12###
  }                                                                                      ###15OCT12###
  #--

  list.entry <- 1
  isNMRpdb("model.pdb", model.num=model)

  # ###### Parsing the PDB file: ######
  if(pdbtype=="complex"){
    OUT <- c(OUT, "NOTE: Parsing with a 'complex' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.complex("model.pdb")$PROTEIN
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="almost"){
    OUT <- c(OUT, "NOTE: Parsing with an 'almost' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.almost("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="medium"){
    OUT <- c(OUT, "NOTE: Parsing with a 'medium' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.medium("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="simple"){  # At present is the same as medium
    OUT <- c(OUT, "NOTE: Parsing with a 'simple' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.simple("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  # Converting the $name convention to the one used in ALMOST, which is the
  #    first record in topology.lib:
  OUT <- c(OUT, "NOTE: Converting the PDB naming convention."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  # assigning "UN" for the chainIDs that are either "" or NA:                    ###18AUG11###
  repair.chain.ind <- which(parsed.pdb$chainID=="" | is.na(parsed.pdb$chainID))  ###18AUG11###
  if(length(repair.chain.ind) > 0){                                              ###18AUG11###
    parsed.pdb$chainID[repair.chain.ind] <- "UN"                                 ###18AUG11###
  }; rm(repair.chain.ind)                                                        ###18AUG11###
  #                                                                              ###18AUG11###
  record.variants <- pdbProc(parsed.pdb, topology=proc.topol)
  parsed.pdb <- convert.topology(parsed.pdb=parsed.pdb, proc.topol=proc.topol,
                                 record.variants=record.variants)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # # # # # # # # # # # # # # # # 
  if(count.tried.model==1){
    # Finding out the Phe/Tyr group presence:
    # Will be done for only the first model, and of course the rest are
    # assumed to be the same:
    Ar.res <- Ar.seq <- Ar.chain <- NULL
    phe.pres <- amac.find(parsed.pdb, AmAc="PHE")
    if(phe.pres$num!=0){
      Ar.res   <- c(Ar.res, rep(phe.pres$resName, phe.pres$num))
      Ar.seq   <- c(Ar.seq, phe.pres$resSeq)
      Ar.chain <- c(Ar.chain, phe.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",phe.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    tyr.pres <- amac.find(parsed.pdb, AmAc="TYR")
    if(tyr.pres$num!=0){
      Ar.res   <- c(Ar.res, rep(tyr.pres$resName, tyr.pres$num))
      Ar.seq   <- c(Ar.seq, tyr.pres$resSeq)
      Ar.chain <- c(Ar.chain, tyr.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",tyr.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }

  } else { #of if(count.tried.model==1)
    OUT <- c(OUT, paste("NOTE: Model ",model," of ",isNMR$Nmodels,".",sep="")); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  }
  # # # # # # # # # # # # # # # # 

  OUT <- c(OUT, "NOTE: Extracting the terms for prediction..."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)


  OUT <- c(OUT, "NOTE: Predicting chemical shifts..."); print(OUT[length(OUT)])
   write(OUT[length(OUT)], file="process_info.txt", append=TRUE)


  for(ar.ind in 1:length(Ar.seq)){

    if(Ar.res[ar.ind]=="PHE"){
        # nucl="HD"
        a <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="PHE", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HD1", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           intrch <- a$interchain
           a <- a$DataToFit
        b <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="PHE", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HD2", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           b <- b$DataToFit
        av.df <- average.pars(DataToFit.m1=a, DataToFit.m2=b, proc.topol=proc.topol)
        av.df <- addCos(Data=av.df, ang="PHI")
        av.df <- addCos(Data=av.df, ang="PSI")
        av.df <- addCos(Data=av.df, ang="CHI1")
        av.df <- addCos(Data=av.df, ang="CHI2")
        av.df <- addMergeDist(Data=av.df)
        av.df <- addROT(Data=av.df)
        #options(warn=-1)
        pr.lm <- predict.lm(object=PHE_HD, newdata=av.df, type='response', se.fit=TRUE)
        rm(a, b, av.df)
        # Filtering anomalous results, by setting them to NA:
        if(pr.lm$fit<PHE_HD_mincs | pr.lm$fit>PHE_HD_maxcs){pr.lm$fit <- NA}
    
       if(count.tried.model==1){
        RESID <- c(RESID, Ar.res[ar.ind])
        SEQ <- c(SEQ, Ar.seq[ar.ind])
        NUCL <- c(NUCL, "HD")
        CHAINID <- c(CHAINID, Ar.chain[ar.ind])
        INTERCHAIN.firstM <- c(INTERCHAIN.firstM, intrch)
       }
        ChSh.pred[[list.entry]] <- c(ChSh.pred[[list.entry]], as.vector(pr.lm$fit))
        SE_fit[[list.entry]] <- c(SE_fit[[list.entry]], as.vector(pr.lm$se.fit))
        SD_train[[list.entry]] <- c(SD_train[[list.entry]], as.vector(pr.lm$residual.scale))
        
        list.entry <- list.entry + 1

        # nucl="HE"
        a <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="PHE", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HE1", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           intrch <- a$interchain
           a <- a$DataToFit
        b <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="PHE", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HE2", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           b <- b$DataToFit
        av.df <- average.pars(DataToFit.m1=a, DataToFit.m2=b, proc.topol=proc.topol)
        av.df <- addCos(Data=av.df, ang="PHI")
        av.df <- addCos(Data=av.df, ang="PSI")
        av.df <- addCos(Data=av.df, ang="CHI1")
        av.df <- addCos(Data=av.df, ang="CHI2")
        av.df <- addMergeDist(Data=av.df)
        av.df <- addROT(Data=av.df)
        #options(warn=-1)
        pr.lm <- predict.lm(object=PHE_HE, newdata=av.df, type='response', se.fit=TRUE)
        rm(a, b, av.df)
        # Filtering anomalous results, by setting them to NA:
        if(pr.lm$fit<PHE_HE_mincs | pr.lm$fit>PHE_HE_maxcs){pr.lm$fit <- NA}

       if(count.tried.model==1){
        RESID <- c(RESID, Ar.res[ar.ind])
        SEQ <- c(SEQ, Ar.seq[ar.ind])
        NUCL <- c(NUCL, "HE")
        CHAINID <- c(CHAINID, Ar.chain[ar.ind])
        INTERCHAIN.firstM <- c(INTERCHAIN.firstM, intrch)
       }
        ChSh.pred[[list.entry]] <- c(ChSh.pred[[list.entry]], as.vector(pr.lm$fit))
        SE_fit[[list.entry]] <- c(SE_fit[[list.entry]], as.vector(pr.lm$se.fit))
        SD_train[[list.entry]] <- c(SD_train[[list.entry]], as.vector(pr.lm$residual.scale))

        list.entry <- list.entry + 1

        # nucl="HZ"
        av.df <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="PHE", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HZ", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           intrch <- av.df$interchain
           av.df <- av.df$DataToFit
        av.df <- addCos(Data=av.df, ang="PHI")
        av.df <- addCos(Data=av.df, ang="PSI")
        av.df <- addCos(Data=av.df, ang="CHI1")
        av.df <- addCos(Data=av.df, ang="CHI2")
        av.df <- addMergeDist(Data=av.df)
        av.df <- addROT(Data=av.df)
        #options(warn=-1)
        pr.lm <- predict.lm(object=PHE_HZ, newdata=av.df, type='response', se.fit=TRUE)
        rm(av.df)
        # Filtering anomalous results, by setting them to NA:
        if(pr.lm$fit<PHE_HZ_mincs | pr.lm$fit>PHE_HZ_maxcs){pr.lm$fit <- NA}

       if(count.tried.model==1){
        RESID <- c(RESID, Ar.res[ar.ind])
        SEQ <- c(SEQ, Ar.seq[ar.ind])
        NUCL <- c(NUCL, "HZ")
        CHAINID <- c(CHAINID, Ar.chain[ar.ind])
        INTERCHAIN.firstM <- c(INTERCHAIN.firstM, intrch)
       }
        ChSh.pred[[list.entry]] <- c(ChSh.pred[[list.entry]], as.vector(pr.lm$fit))
        SE_fit[[list.entry]] <- c(SE_fit[[list.entry]], as.vector(pr.lm$se.fit))
        SD_train[[list.entry]] <- c(SD_train[[list.entry]], as.vector(pr.lm$residual.scale))

        list.entry <- list.entry + 1 
    }

    if(Ar.res[ar.ind]=="TYR"){
        # nucl="HD"
        a <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="TYR", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HD1", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           intrch <- a$interchain
           a <- a$DataToFit
        b <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="TYR", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HD2", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           b <- b$DataToFit
        av.df <- average.pars(DataToFit.m1=a, DataToFit.m2=b, proc.topol=proc.topol)
        av.df <- addCos(Data=av.df, ang="PHI")
        av.df <- addCos(Data=av.df, ang="PSI")
        av.df <- addCos(Data=av.df, ang="CHI1")
        av.df <- addCos(Data=av.df, ang="CHI2")
        av.df <- addMergeDist(Data=av.df)
        av.df <- addROT(Data=av.df)
        #options(warn=-1)
        pr.lm <- predict.lm(object=TYR_HD, newdata=av.df, type='response', se.fit=TRUE)
        rm(a, b, av.df)
        # Filtering anomalous results, by setting them to NA:
        if(pr.lm$fit<TYR_HD_mincs | pr.lm$fit>TYR_HD_maxcs){pr.lm$fit <- NA}
    
       if(count.tried.model==1){
        RESID <- c(RESID, Ar.res[ar.ind])
        SEQ <- c(SEQ, Ar.seq[ar.ind])
        NUCL <- c(NUCL, "HD")
        CHAINID <- c(CHAINID, Ar.chain[ar.ind])
        INTERCHAIN.firstM <- c(INTERCHAIN.firstM, intrch)
       }
        ChSh.pred[[list.entry]] <- c(ChSh.pred[[list.entry]], as.vector(pr.lm$fit))
        SE_fit[[list.entry]] <- c(SE_fit[[list.entry]], as.vector(pr.lm$se.fit))
        SD_train[[list.entry]] <- c(SD_train[[list.entry]], as.vector(pr.lm$residual.scale))
        
        list.entry <- list.entry + 1

        # nucl="HE"
        a <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="TYR", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HE1", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           intrch <- a$interchain
           a <- a$DataToFit
        b <- Ar.get.single.dframe(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc="TYR", resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HE2", readlib.dih=readlib.dih, arealR=6.5, chainID=Ar.chain[ar.ind])
           b <- b$DataToFit
        av.df <- average.pars(DataToFit.m1=a, DataToFit.m2=b, proc.topol=proc.topol)
        av.df <- addCos(Data=av.df, ang="PHI")
        av.df <- addCos(Data=av.df, ang="PSI")
        av.df <- addCos(Data=av.df, ang="CHI1")
        av.df <- addCos(Data=av.df, ang="CHI2")
        av.df <- addMergeDist(Data=av.df)
        av.df <- addROT(Data=av.df)
        #options(warn=-1)
        pr.lm <- predict.lm(object=TYR_HE, newdata=av.df, type='response', se.fit=TRUE)
        rm(a, b, av.df)
        # Filtering anomalous results, by setting them to NA:
        if(pr.lm$fit<TYR_HE_mincs | pr.lm$fit>TYR_HE_maxcs){pr.lm$fit <- NA}

       if(count.tried.model==1){
        RESID <- c(RESID, Ar.res[ar.ind])
        SEQ <- c(SEQ, Ar.seq[ar.ind])
        NUCL <- c(NUCL, "HE")
        CHAINID <- c(CHAINID, Ar.chain[ar.ind])
        INTERCHAIN.firstM <- c(INTERCHAIN.firstM, intrch)
       }
        ChSh.pred[[list.entry]] <- c(ChSh.pred[[list.entry]], as.vector(pr.lm$fit))
        SE_fit[[list.entry]] <- c(SE_fit[[list.entry]], as.vector(pr.lm$se.fit))
        SD_train[[list.entry]] <- c(SD_train[[list.entry]], as.vector(pr.lm$residual.scale))

        list.entry <- list.entry + 1
    }
  } # of   for(ar.ind in 1:length(Ar.seq))

  count.tried.model <- count.tried.model + 1
} # of for(model in models2analyse)

file.remove("model.pdb")

#########################################################################
list.entry <- list.entry-1
#RESID; SEQ; NUCL; CHAINID
ChSh.pred <- ChSh.pred[1:list.entry]
SE_fit    <-    SE_fit[1:list.entry]
SD_train  <-  SD_train[1:list.entry]

# Averaging the data obtained for multiple conformers
Pred.CS               <- sapply(ChSh.pred, FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)
Pred.sefit            <- sapply(SE_fit,    FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)
Pred.SDresidual.train <- sapply(SD_train,  FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)


#########################################################################
## Filtering out the erroneous and violating predictions + reformating data

violating.shift.ind <- which(is.na(Pred.CS))
if(length(violating.shift.ind)>0){
  Pred.CS[violating.shift.ind] <- "xxxxx"
  Pred.CS[-violating.shift.ind] <- format(round(as.numeric(Pred.CS[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.CS <- format(round(Pred.CS,3),nsmall=3)
}

violating.shift.ind <- which(Pred.sefit>999.99)
if(length(violating.shift.ind)>0){
  Pred.sefit[violating.shift.ind] <- "xxxxx"
  Pred.sefit[-violating.shift.ind] <- format(round(as.numeric(Pred.sefit[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.sefit <- format(round(Pred.sefit,3),nsmall=3)
}

violating.shift.ind <- which(Pred.SDresidual.train>999.99)
if(length(violating.shift.ind)>0){
  Pred.SDresidual.train[violating.shift.ind] <- "xxxxx"
  Pred.SDresidual.train[-violating.shift.ind] <- format(round(
                                      as.numeric(Pred.SDresidual.train[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.SDresidual.train <- format(round(as.numeric(Pred.SDresidual.train),3),nsmall=3)
}

# RESID; SEQ; NUCL; Pred.CS; Pred.sefit; Pred.SDresidual.train

#########################################################################
# # # # Reading in the experimental data
if(length(experdata)!=0){
  OUT <- c(OUT, "NOTE: Processing the experimental data."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  exp.data <- readLines(experdata)

  # dummy and blank line exclusion is also added.
  comment.rows <- unique(c(grep("#", exp.data, fixed=TRUE),  
                           which(exp.data==""),
                           which(exp.data==" "),
                           which(exp.data=="  "),
                           which(exp.data=="   "),
                           which(exp.data=="    "),
                           which(exp.data=="     "),
                           which(exp.data=="      ")  ))

  if(length(comment.rows!=0)){ 
    exp.data <- exp.data[-comment.rows] 
  }

  exp.amac <- exp.seq <- exp.nucl <- exp.shift <- exp.chain <- NULL

  # # # # 
  for(i in 1:length(exp.data)){
    line <- linesplit(exp.data[i])
    if(length(line)==4 | length(line)==5){
      amacic    <- as.character(line[1])
      exp.amac  <- c(exp.amac, amacic)
      exp.seq   <- c(exp.seq, as.numeric(line[2]))
      nuclic    <- as.character(line[3])
      # switching to internal convention:
      if(amacic=="PHE"& nuclic=="HD1"){nuclic<-"HD"}     
      if(amacic=="PHE"& nuclic=="HD2"){nuclic<-"HD"}
      if(amacic=="PHE"& nuclic=="HE1"){nuclic<-"HE"}
      if(amacic=="PHE"& nuclic=="HE2"){nuclic<-"HE"}
      if(amacic=="PHE"& nuclic=="HZ") {nuclic<-"HZ"}
      if(amacic=="TYR"& nuclic=="HD1"){nuclic<-"HD"}
      if(amacic=="TYR"& nuclic=="HD2"){nuclic<-"HD"}
      if(amacic=="TYR"& nuclic=="HE1"){nuclic<-"HE"}
      if(amacic=="TYR"& nuclic=="HE2"){nuclic<-"HE"}
      exp.nucl  <- c(exp.nucl, nuclic)
      exp.shift <- c(exp.shift, as.numeric(line[4]))
      chainic   <- as.character(line[5])
      if(is.na(chainic)){chainic <- "UN"}
      exp.chain <- c(exp.chain, chainic)
    }
  };rm(i)
  # # # # 

  Exp.CS <- Exp_Pred.CS <- rep(NA, length(NUCL))

  for(i in 1:length(exp.shift)) {
    position <-  which( RESID   == exp.amac[i]  &
                        NUCL    == exp.nucl[i]  & 
                        SEQ     == exp.seq[i]   &
                        CHAINID == exp.chain[i]   )
    if(length(position)!=0){
      Exp.CS[position] <- round(exp.shift[i], 3)
      prediction <- as.numeric(Pred.CS[position])
      if(!is.na(prediction)){
        Exp_Pred.CS[position] <- exp.shift[i]-prediction
      }
    } else {
      OUT <- c(OUT, "ERROR: Some of the exp. data records do not match the structure,") ###15OCT12###
      OUT <- c(OUT, "ERROR: or the experimental data file is incorrectly formatted.")   ###15OCT12###
      OUT <- c(OUT, paste("ERROR: i=",i,", RESID=",exp.amac[i],", NUCL=",exp.nucl[i],   ###15OCT12###
                          ", SEQ=",exp.seq[i],", CHAINID=",exp.chain[i],sep=""))        ###15OCT12###
      write(OUT[c(length(OUT)-2,length(OUT)-1,length(OUT))],                            ###15OCT12###
           file="process_info.txt", append=TRUE)                                        ###15OCT12###
      stop(OUT[c(length(OUT)-2,length(OUT)-1,length(OUT))])                             ###15OCT12###
    } 
  }

  Exp.CS.nonfilt <- Exp.CS

  ndx <- which(is.na(Exp.CS))
  if(length(ndx)!=0){
   Exp.CS[ndx]  <- "xxxxx"
   Exp.CS[-ndx] <- format(round(as.numeric(Exp.CS[-ndx]),3),nsmall=3)
  } else {
   Exp.CS <- format(round(as.numeric(Exp.CS),3),nsmall=3)  
  }
  ndx <- which(is.na(Exp_Pred.CS))
  if(length(ndx)!=0){
   Exp_Pred.CS[ndx]  <- "xxxxxx"
   Exp_Pred.CS[-ndx] <- format(round(as.numeric(Exp_Pred.CS[-ndx]),3),nsmall=3)
  } else {
   Exp_Pred.CS <- format(round(as.numeric(Exp_Pred.CS),3),nsmall=3)  
  }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # # # REREFERENCING
  if(rereference==TRUE){
    Pred.CS.num <- as.numeric(Pred.CS)
    # rereferencing proton chemical shifts:
    h.reref.ind <- which(!is.na(Pred.CS.num) & !is.na(Exp.CS.nonfilt))
    exp.reref  <- Exp.CS.nonfilt[h.reref.ind]
    pred.reref <- Pred.CS.num[h.reref.ind]
    if(length(exp.reref)>minRerefNum) {
      lm.reref <- lm(formula = exp.reref ~ 1, offset = pred.reref)
      Pred.CS.num <- Pred.CS.num + lm.reref$coefficients[1]
    }
    # recalculating errors
    Exp_Pred.CS <- Exp.CS.nonfilt-Pred.CS.num

    # reformating back
    na.pred.CS <- which(is.na(Pred.CS.num))
    if(length(na.pred.CS)>0){
      Pred.CS.num[na.pred.CS] <- "xxxxx"
      Pred.CS.num[-na.pred.CS] <- format(round(as.numeric(Pred.CS.num[-na.pred.CS]),3),nsmall=3)
    } else {
      Pred.CS.num <- format(round(as.numeric(Pred.CS.num),3),nsmall=3)
    }
    Pred.CS <- Pred.CS.num
    na.err.CS <- which(is.na(Exp_Pred.CS))
    if(length(na.err.CS)>0){
     Exp_Pred.CS[na.err.CS]  <- "xxxxxx" 
     Exp_Pred.CS[-na.err.CS] <- format(round(as.numeric(Exp_Pred.CS[-na.err.CS]),3),nsmall=3)
    } else {
     Exp_Pred.CS <- format(round(as.numeric(Exp_Pred.CS),3),nsmall=3)  
    }
    #
  }
  # # # END OF REREFERENCING

  # RESID; SEQ; NUCL; Pred.CS; Pred.sefit; Pred.SDresidual.train; Exp.CS; Exp_Pred.CS

  # Generating the Probs and Stars data objects that comment on probabilities
  Exp_Pred.CS.num <- as.vector(sapply(Exp_Pred.CS, fun5, simplify=TRUE))
  Probs <- rep("xxxxx", length(Exp_Pred.CS.num))
  Stars <- rep("xxx", length(Exp_Pred.CS.num))
  UNQ.SEQ.CHAINID <- paste(SEQ, CHAINID, sep="_")

  for(SEQ.CH.i in unique(UNQ.SEQ.CHAINID)) {
    is <- which(SEQ==as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]) &
                CHAINID==unlist(strsplit(SEQ.CH.i,"_"))[2] & 
                !is.na(Exp_Pred.CS.num))  
    if(length(is)!=0) {
      P.data <- getP(AMC=ResidueName(unique(RESID[is]), capitalize=FALSE,
                     threeletter=TRUE, oneletter = FALSE),
                     NCL=NUCL[is], De_c=Exp_Pred.CS.num[is], er.data=er.data)
      Probs[is] <- format(round(as.numeric(P.data$Prob),3),nsmall=3)  
      Stars[is] <- paste(rep("*", P.data$Stars), collapse="")
      rm(P.data)
    }
    rm(is)   
  }; rm(SEQ.CH.i)

  # # #
  if(doval==TRUE){

     # Generating chimera_cmd.txt file with UCSF Chimera command to color residues
     # according to their probabilities:
     colfun <- colorRampPalette(c("red","orange","yellow","green",5,"blue"), bias=bias)
     colors <- colfun(500)
     prob.range <- seq(from=0, to=1, by=0.002001)
     # plotting pallete:
     if(x11lib==TRUE){
       jpeg(quality=100, height=500, width=500, filename="protein_prob_palette.jpg")
         plot(x=prob.range,y=rep(0,500), col=colors, cex=5, pch=3, ylim=c(-0.1,0.1),
              main=paste("Bias = ", bias,sep=""))
       dev.off()
     }
     #
     chimera.cmd <- NULL
     for(SEQ.CH.i in unique(UNQ.SEQ.CHAINID)){
       is <- which(SEQ==as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]) & 
                   CHAINID==unlist(strsplit(SEQ.CH.i,"_"))[2] &
                   Probs!="xxxxx")
       if(length(is)!=0) {
         is <- is[1]
         prb <- as.numeric(Probs[is]) 
         colind <- which(abs(prob.range-prb) == min(abs(prob.range-prb)))  
         chimera.cmd <- c(chimera.cmd," color ",colors[colind],
                          " :",as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]),
                           ".",unlist(strsplit(SEQ.CH.i,"_"))[2],";")
         rm(prb, colind)
       }
       rm(is)
     }; rm(SEQ.CH.i, colors, colfun, prob.range)
     
     write(paste(chimera.cmd, collapse=""), file="chimera_cmd.txt")
     rm(chimera.cmd)
  }
  # # #
 
  # Printing out the results:
  if(outorder=="bytype"){ # Prints results odered by type.
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")," ; ",
          format(CHAINID, width=3, justify="centre")," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre"),"; ",
          format(NUCL, width=4, justify="centre")," ; ",
          format(Pred.CS, width=6, justify="centre")," ; ",
          format(Pred.sefit, width=5, justify="centre")," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")," ; ",
          format(Exp.CS, width=6, justify="centre")," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")," ; ",
          format(Probs, width=5, justify="centre")," ; ",
          format(Stars, width=3, justify="centre")
          ,sep=""))
  }
  if(outorder=="byseq"){ # Prints results odered by sequence number.
    byseq.ord <- order(SEQ)
    if(length(unique(CHAINID)) > 1){
      byseq.ord <- byseq.ord[order(CHAINID[byseq.ord])]
    }
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byseq.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byseq.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byseq.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byseq.ord],"; ",
          format(NUCL, width=4, justify="centre")[byseq.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byseq.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byseq.ord]," ; ",
          format(Exp.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Probs, width=5, justify="centre")[byseq.ord]," ; ",
          format(Stars, width=3, justify="centre")[byseq.ord]
          ,sep=""))
  }
  if(outorder=="byshift"){ # Prints results odered by sequence number.
    byshift.ord <- order(as.numeric(Pred.CS))
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byshift.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byshift.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byshift.ord],"; ",
          format(NUCL, width=4, justify="centre")[byshift.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byshift.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byshift.ord]," ; ",
          format(Exp.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Probs, width=5, justify="centre")[byshift.ord]," ; ",
          format(Stars, width=3, justify="centre")[byshift.ord]
          ,sep=""))
  }

  rm(Probs, Stars, Exp_Pred.CS.num, fun5, er.data, getP)

} else {  # of if(length(experdata)!=0) => experimental data does not exist!!!


  if(outorder=="bytype"){ # Prints results odered by type.
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")," ; ",
          format(CHAINID, width=3, justify="centre")," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre"),"; ",
          format(NUCL, width=4, justify="centre")," ; ",
          format(Pred.CS, width=6, justify="centre")," ; ",
          format(Pred.sefit, width=5, justify="centre")," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")
          ,sep=""))
  }
  if(outorder=="byseq"){ # Prints results odered by thype.
    byseq.ord <- order(SEQ)
    if(length(unique(CHAINID)) > 1){
      byseq.ord <- byseq.ord[order(CHAINID[byseq.ord])]
    }
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byseq.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byseq.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byseq.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byseq.ord],"; ",
          format(NUCL, width=4, justify="centre")[byseq.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byseq.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byseq.ord]
          ,sep=""))
  }
  if(outorder=="byshift"){ # Prints results odered by thype.
    byshift.ord <- order(as.numeric(Pred.CS))
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byshift.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byshift.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byshift.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byshift.ord],"; ",
          format(NUCL, width=4, justify="centre")[byshift.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byshift.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byshift.ord]
          ,sep=""))
  }
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

OUT <- c(OUT, RESULTS.rep, 
         "NOTE: ArShift - Good luck!                    ",
         "NOTE: For questions and bug reports contact - aleksahak@cantab.net ",
         "NOTE: N.B. Please, attach all the related files to the e-mail.  ")

write(OUT, file=outputname)
print("NOTE: Output file is written!")
write("NOTE: Output file is written!", file="process_info.txt", append=TRUE)
  
# DONE (P added)
