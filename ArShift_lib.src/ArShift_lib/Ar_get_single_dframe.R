################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
# VERSION: 18 August 2011
# VERSION: 28 October 2010, modified to account chainIDs:

Ar.get.single.dframe <- function(parsed.pdb=parsed.pdb, proc.topol=proc.topol, record.variants=record.variants, AmAc=Ar.res[ar.ind], resSeqnum=Ar.seq[ar.ind], whichNMRnucl="HD1", readlib.dih=readlib.dih, arealR=6.5, chainID){

   # Generating the distance cell names, which will hold the further distance data:
   amac.list <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "VAL", "SER", "THR",
                                              "TYR", "TRP", "GLY", "PRO") 
   obj <- NULL
   for(i in amac.list) {
     i1 <- i
     if(i == "CYS"){i1 <- "CYN"}
     if(i == "HIS"){i1 <- "HID"}  #-- will be initialized as HIS but with atom names of HID                    ###18AUG11###
                                  #-- so that additional step is needed to create a dummy objects for "HE2",   ###18AUG11###
                                  #-- in case the users forget to convert all HIS to HID, so that HE2 remains. ###18AUG11###
   
     m <- topol.retrieve(AmAc=i1,proc.topol=proc.topol, nvar=1)
     resid <- m[,"residues"]
     rec <- m[,"record"]
     rec <- c(rec, "NAN")                                                                                      ###18AUG11###
     if(i == "HIS"){ rec <- c(rec,"HE2") }                                                                     ###18AUG11###
   
     for(k in 1:length(rec)){
       obj <- c(obj, paste(i, ".", rec[k],sep=""),
                     paste(i, ".", rec[k],".R1",sep=""), 
                     paste(i, ".", rec[k],".R3",sep=""),
                     paste(i, ".", rec[k],".R6",sep="")  )
     }
   }; distance.cell.names <- obj; rm(resid,m,rec,i,i1,k,obj)
   # # # #

   # Generating the cell names for the other data holders
   dihedral.cell.names <- c("PHI","PSI","CHI1","CHI2")
   ring.cell.names     <- c("RING.HIS","RING.TRP5","RING.TRP6","RING.TYR","RING.PHE")
   aniso.cell.names    <- c("ANISO.ASP","ANISO.ASN","ANISO.GLU","ANISO.GLN","ANISO.ARG",
                            "ANISO.PEPT")
   ef.cell.names       <- c("EF.FARCC")
   info.cell.names     <- c("AMACID", "NMRNUCL")
   # # # #

   # All the cell names, ready for initialization.
   all.cell.names <- c(dihedral.cell.names, ring.cell.names, aniso.cell.names,
                       ef.cell.names, distance.cell.names, info.cell.names)
   # # # #

   # Initializing all the cells as 0:
   for(i in all.cell.names){
     eval(parse(text=paste(i, "<- 0",sep="")))
   }; rm(i)
   # # # #

  
  dihedrals <- all.dihedral(resSeqnum, parsed.pdb, READLIB=readlib.dih, chainID=chainID)
           
  query.xyz <- get.coord(name=whichNMRnucl, resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)

  areal <- getarea(xyz=query.xyz, parsed.pdb=parsed.pdb, r=arealR, charge=TRUE, proc.topol=proc.topol)

  DistDepos <- Ar.Merge.Same.Type.dist(arealt=areal, resSeqnum=resSeqnum, amac.type=AmAc, chainID=chainID)

  # Storing distance information:
  for(i in 1:length(DistDepos$name)) {
   if(DistDepos$name[i]!="OT1" & DistDepos$name[i]!="OT2" &
      DistDepos$name[i]!="HT1" & DistDepos$name[i]!="HT2" &
      DistDepos$name[i]!="HT3") {  
    if(DistDepos$name[i]!="HN2" & DistDepos$resName[i]!="PRO") { 
      eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],"<-",DistDepos$Sdist[i],sep=""))) 
      eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R1<-",DistDepos$SdistR1[i],sep=""))) 
      eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R3<-",DistDepos$SdistR3[i],sep="")))
      eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R6<-",DistDepos$SdistR6[i],sep=""))) 
    }
   }
  }; rm(i) # all the distances are stored

  # Storing identification information
  AMACID  <- ResidueName(AmAc, threeletter=TRUE)
  NMRNUCL <- whichNMRnucl

  # Storing dihedral angle data
  PHI  <- dihedrals[1]
  PSI  <- dihedrals[2]
  CHI1 <- dihedrals[4]
  chi2nosym <- dihedrals[5]
  if(chi2nosym < 0){ chi2nosym <- chi2nosym + 180 } # correcting for symmetry
  CHI2 <- chi2nosym

    # Storing COO and CON anisotropy data
    asp <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                       atom1="CG", atom2="OD1", atom3="OD2", Amino="ASP")$ANISO.SUM
    if(length(asp)!=0){ ANISO.ASP <- asp }
    #
    asn <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                       atom1="CG", atom2="OD1", atom3="ND2", Amino="ASN")$ANISO.SUM
    if(length(asn)!=0){ ANISO.ASN <- asn }
    #
    glu <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                       atom1="CD", atom2="OE1", atom3="OE2", Amino="GLU")$ANISO.SUM
    if(length(glu)!=0){ ANISO.GLU <- glu }
    #
    gln <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                       atom1="CD", atom2="OE1", atom3="NE2", Amino="GLN")$ANISO.SUM
    if(length(gln)!=0){ ANISO.GLN <- gln }
    #
    arg <- Aniso.Arg(query.xyz, areal, parsed.pdb, 
                       atom1="CZ", atom2="NH1", atom3="NH2", Amino="ARG")$ANISO.SUM
    if(length(arg)!=0){ ANISO.ARG <- arg }
    rm(asp,asn,glu,gln,arg)
  
    # Storing peptide moiety anisotropy data
    An.Ppt <- Aniso.Pept(query.xyz, resSeqnum, areal, parsed.pdb, chainID=chainID) 
    if(length(An.Ppt$ANISO.PEPT)!=0) { ANISO.PEPT <- An.Ppt$ANISO.PEPT }
  
    # Storing ring current data
    his.ring <- Ar.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="HIS", excludeSeqID=resSeqnum, chainID=chainID)$RC.SUM 
     if(length(his.ring)!=0){ RING.HIS <- his.ring }
    trp5.ring <- Ar.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TRP", excludeSeqID=resSeqnum, chainID=chainID)$RC1.SUM 
     if(length(trp5.ring)!=0){ RING.TRP5 <- trp5.ring }
    trp6.ring <- Ar.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TRP", excludeSeqID=resSeqnum, chainID=chainID)$RC2.SUM 
     if(length(trp6.ring)!=0){ RING.TRP6 <- trp6.ring }
    tyr.ring <- Ar.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TYR", excludeSeqID=resSeqnum, chainID=chainID)$RC.SUM 
     if(length(tyr.ring)!=0){ RING.TYR <- tyr.ring }
    phe.ring <- Ar.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="PHE", excludeSeqID=resSeqnum, chainID=chainID)$RC.SUM
     if(length(phe.ring)!=0){ RING.PHE <- phe.ring }
    rm(his.ring, trp5.ring, trp6.ring, tyr.ring, phe.ring)
  
  # Storing electric field data
  EF.FARCC <- ARgetEF(resSeqnum, AmAc, query.xyz, whichNMRnucl, areal, parsed.pdb, chainID=chainID)$FARCC

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # NA check for validation (is extra)
  PHI[is.na(PHI)]               <- 0.00
  PSI[is.na(PSI)]               <- 0.00
  CHI1[is.na(CHI1)]             <- 0.00                                                        
  CHI2[is.na(CHI2)]             <- 0.00                                       
  RING.HIS[is.na(RING.HIS)]     <- 0.00   
  RING.TRP5[is.na(RING.TRP5)]   <- 0.00
  RING.TRP6[is.na(RING.TRP6)]   <- 0.00
  RING.TYR[is.na(RING.TYR)]     <- 0.00   
  RING.PHE[is.na(RING.PHE)]     <- 0.00     
  ANISO.ASP[is.na(ANISO.ASP)]   <- 0.00
  ANISO.ASN[is.na(ANISO.ASN)]   <- 0.00
  ANISO.GLU[is.na(ANISO.GLU)]   <- 0.00
  ANISO.GLN[is.na(ANISO.GLN)]   <- 0.00
  ANISO.ARG[is.na(ANISO.ARG)]   <- 0.00
  ANISO.PEPT[is.na(ANISO.PEPT)] <- 0.00
  EF.FARCC[is.na(EF.FARCC)]     <- 0.00             
                                                                
  # Generating the final data frame:
  DataToFit <- data.frame(  AMACID         =as.character(AMACID),
                            NMRNUCL        =as.character(NMRNUCL),
                            EF.FARCC       =as.numeric(EF.FARCC),
                            ANISO.ARG      =as.numeric(ANISO.ARG),
                            ANISO.PEPT     =as.numeric(ANISO.PEPT), 
                            ANISO.AMD      =as.numeric(ANISO.ASN)+as.numeric(ANISO.GLN), ###
                            ANISO.ACD      =as.numeric(ANISO.ASP)+as.numeric(ANISO.GLU), ###
                            RING.HIS       =as.numeric(RING.HIS),
                            RING.TRP5      =as.numeric(RING.TRP5),
                            RING.TRP6      =as.numeric(RING.TRP6),
                            RING.TYR       =as.numeric(RING.TYR),
                            RING.PHE       =as.numeric(RING.PHE),
                            PHI            =(as.numeric(PHI)*pi/180),
                            PHI2           =(as.numeric(PHI)*pi/180)^2,
                            PHI3           =(as.numeric(PHI)*pi/180)^3,
                            PHI4           =(as.numeric(PHI)*pi/180)^4,
                            PSI            =(as.numeric(PSI)*pi/180),
                            PSI2           =(as.numeric(PSI)*pi/180)^2,
                            PSI3           =(as.numeric(PSI)*pi/180)^3,
                            PSI4           =(as.numeric(PSI)*pi/180)^4,
                            CHI1           =as.numeric(CHI1)*pi/180, 
                            CHI12          =(as.numeric(CHI1)*pi/180)^2,
                            CHI13          =(as.numeric(CHI1)*pi/180)^3, 
                            CHI14          =(as.numeric(CHI1)*pi/180)^4, 
                            CHI2           =as.numeric(CHI2)*pi/180,
                            CHI22          =(as.numeric(CHI2)*pi/180)^2,
                            CHI23          =(as.numeric(CHI2)*pi/180)^3,
                            CHI24          =(as.numeric(CHI2)*pi/180)^4  )

  # Adding all the distance information (in Angstrom) to the created data frame.
  text.cmd <- "DataToFit <- data.frame(DataToFit"
  for(i in distance.cell.names) {
    cmd <- paste(i,"[is.na(",i,")] <- 0.00",sep="")
    eval(parse(text=cmd)); rm(cmd)
    text.cmd <- paste(text.cmd,",",i,"=",i,sep="")
  }; rm(i)
  text.cmd <- paste(text.cmd,")",sep="")
  eval(parse(text=text.cmd)); rm(text.cmd)
  # # # #


  DATA.TO.RETURN <- NULL
  DATA.TO.RETURN$DataToFit <- DataToFit
  DATA.TO.RETURN$interchain <- areal$interchain
  return(DATA.TO.RETURN)

}

################################################################################
