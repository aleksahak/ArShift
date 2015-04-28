################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

average.pars <- function(DataToFit.m1, DataToFit.m2, proc.topol=proc.topol){

  
    DataToFit <- NULL      
    row <- ( as.vector(as.numeric(as.matrix(DataToFit.m1[1,3:1408]))) +
             as.vector(as.numeric(as.matrix(DataToFit.m2[1,3:1408]))) )/2       
    row <- c(as.vector(as.matrix(DataToFit.m1[1,1:2])), row)
    DataToFit <- rbind(DataToFit, row)
 
    dimnames(DataToFit)[[2]] <- dimnames(DataToFit.m1)[[2]]
    dimnames(DataToFit)[[1]] <- 1:length(DataToFit[,1])

    DataToFit <- data.frame(DataToFit)


  ### REASSEMBLING THE DATA FRAME FOR REFORMATING:

  # Generating the distance cell names, which will hold the further distance data:
  amac.list <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS",
                "ILE", "LEU", "LYS", "MET", "PHE", "VAL", "SER", "THR",
                "TYR", "TRP", "GLY", "PRO") 
  obj <- NULL
  for(i in amac.list) {
    i1 <- i
    if(i == "CYS"){i1 <- "CYN"}
    if(i == "HIS"){i1 <- "HID"}
  
    m <- topol.retrieve(AmAc=i1,proc.topol=proc.topol, nvar=1)
    resid <- m[,"residues"]
    rec <- m[,"record"]
    for(k in 1:length(rec)){
      obj <- c(obj, paste(i, ".", rec[k],sep=""),
                    paste(i, ".", rec[k],".R1",sep=""), 
                    paste(i, ".", rec[k],".R3",sep=""),
                    paste(i, ".", rec[k],".R6",sep="")  )
    }
  }; distance.cell.names <- obj; rm(resid,m,rec,i,i1,k,obj)
  # # # #

  DataToFitB <- DataToFit; rm(DataToFit)
  DataToFit <- data.frame(  AMACID         =as.character(as.vector(as.matrix(DataToFitB$AMACID))),
                            NMRNUCL        =as.character(as.vector(as.matrix(DataToFitB$NMRNUCL))),
                            EF.FARCC       =as.numeric(as.vector(as.matrix(DataToFitB$EF.FARCC))),
                            ANISO.ARG      =as.numeric(as.vector(as.matrix(DataToFitB$ANISO.ARG))),
                            ANISO.PEPT     =as.numeric(as.vector(as.matrix(DataToFitB$ANISO.PEPT))), 
                            ANISO.AMD      =as.numeric(as.vector(as.matrix(DataToFitB$ANISO.AMD))), ###
                            ANISO.ACD      =as.numeric(as.vector(as.matrix(DataToFitB$ANISO.ACD))), ###
                            RING.HIS       =as.numeric(as.vector(as.matrix(DataToFitB$RING.HIS))),
                            RING.TRP5      =as.numeric(as.vector(as.matrix(DataToFitB$RING.TRP5))),
                            RING.TRP6      =as.numeric(as.vector(as.matrix(DataToFitB$RING.TRP6))),
                            RING.TYR       =as.numeric(as.vector(as.matrix(DataToFitB$RING.TYR))),
                            RING.PHE       =as.numeric(as.vector(as.matrix(DataToFitB$RING.PHE))),
                            PHI            =as.numeric(as.vector(as.matrix(DataToFitB$PHI))),
                            PHI2           =as.numeric(as.vector(as.matrix(DataToFitB$PHI2))),
                            PHI3           =as.numeric(as.vector(as.matrix(DataToFitB$PHI3))),
                            PHI4           =as.numeric(as.vector(as.matrix(DataToFitB$PHI4))),
                            PSI            =as.numeric(as.vector(as.matrix(DataToFitB$PSI))),
                            PSI2           =as.numeric(as.vector(as.matrix(DataToFitB$PSI2))),
                            PSI3           =as.numeric(as.vector(as.matrix(DataToFitB$PSI3))),
                            PSI4           =as.numeric(as.vector(as.matrix(DataToFitB$PSI4))),
                            CHI1           =as.numeric(as.vector(as.matrix(DataToFitB$CHI1))), 
                            CHI12          =as.numeric(as.vector(as.matrix(DataToFitB$CHI12))),
                            CHI13          =as.numeric(as.vector(as.matrix(DataToFitB$CHI13))), 
                            CHI14          =as.numeric(as.vector(as.matrix(DataToFitB$CHI14))), 
                            CHI2           =as.numeric(as.vector(as.matrix(DataToFitB$CHI2))),
                            CHI22          =as.numeric(as.vector(as.matrix(DataToFitB$CHI22))),
                            CHI23          =as.numeric(as.vector(as.matrix(DataToFitB$CHI23))),
                            CHI24          =as.numeric(as.vector(as.matrix(DataToFitB$CHI24)))  )
  
  # Adding all the distance information (in Angstrom) to the created data frame,
  # with only [indi] cells:
  text.cmd <- "DataToFit <- data.frame(DataToFit"
  for(i in distance.cell.names) {
    text.cmd <- paste(text.cmd,",",i,"=as.numeric(as.vector(as.matrix(DataToFitB$",i,")))",sep="")
  }; rm(i)
  text.cmd <- paste(text.cmd,")",sep="")
  eval(parse(text=text.cmd)); rm(text.cmd)
  # # # #

  return(DataToFit)
}

################################################################################
