################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
## This function determines whether the given pdb file contains an NMR        ##
## structure or an Xray one. If the structure is an NMR one, the func-        ##
## tion returns an information on the number of models inside, and if         ##
## second numeric argument is not 0, it modifies the file, so that it         ##
## contains only the specified model (1st, 2nd etc.).                         ##
## Now, takes into account Nuclei acids as well!!                             ##
################################################################################
# VERSION: 23 Sep,2010 - an important bug is fixed: some pdb files include
#    the word MODEL before the actual pdb records, hence the comment part
#    should be trimmed from the pdb file at first.
#
# VERSION: 1 Feb, 2010 + an important bug is fixed
# ARGUMENTS:
# pdb.filename - name of the pdb filename
#
# model.num - index of the model to be outputted for file modification,
#             no modification if model.num=0, or the file is an Xray one.
# OUTPUTS:
# $type - "NMR" or "Xray".
# $isH  - TRUE/FALSE whether contains hydrogens (non-mainchain) or not. 
# $Nmodels - number of models in the file for NMR structures.
################################################################################

isNMRpdb <- function(pdb.filename, model.num=0){

  RESULTS <- NULL
  PDB <- readLines(pdb.filename)
  
    if( length(grep(" R VALUE ", PDB, fixed=TRUE)) > 0 |
        length(grep(" CRYSTALLOGRAPHIC SYMMETRY ", PDB, fixed=TRUE))>0) {
      RESULTS$type <- "Xray"    
    } else {
      RESULTS$type <- "NMR"
    }
    #####################################################################
    if(  length(grep(" HB ", PDB, fixed=TRUE))  > 0 |
         length(grep(" HB1 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" HB2 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" HG2 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" HD2 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" HE2 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" HZ2 ", PDB, fixed=TRUE)) > 0 |
         length(grep(" H4' ", PDB, fixed=TRUE)) > 0 |
         length(grep(" H4* ", PDB, fixed=TRUE)) > 0 |
         length(grep(" H6 ", PDB, fixed=TRUE))  > 0 |
         length(grep(" H8 ", PDB, fixed=TRUE))  > 0 |
         length(grep(" H3' ", PDB, fixed=TRUE)) > 0 |
         length(grep(" H3* ", PDB, fixed=TRUE)) > 0   ) {
      RESULTS$isH <- TRUE
    } else {
      RESULTS$isH <- FALSE
    }
    #####################################################################
    # At first removing some of the rows which may contain the word ATOM
    remove <- grep("REMARK", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("CONECT", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("SEQRES", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("COMPND", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }
  
    remove <- grep("TITLE", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }
  
    remove <- grep("HEADER", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("KEYWDS", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("EXPDTA", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("REVDAT", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("JRNL", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("SOURCE", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("MDLTYP", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("AUTHOR", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("DBREF", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("MODRES", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("HETNAM", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("FORMUL", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("LINK   ", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("HET  ", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("CRYST1", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("SCALE", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }

    remove <- grep("ORIGX", PDB, fixed=TRUE)
    if(length(c(1,remove)) > 1) { PDB <- PDB[-remove] }
    #####################################################################
    MD <- grep("MODEL  ", PDB, fixed=TRUE)
    if( (length(MD) == 0) | (RESULTS$type == "Xray") ) { 
      RESULTS$Nmodels <- 1 
    } else {
      false.MD <- grep(" MODEL  ", PDB, fixed=TRUE)
      real.MD <- MD[is.na(match(MD, false.MD))]
      RESULTS$Nmodels <- length(real.MD)
      if( model.num > 0 & model.num != RESULTS$Nmodels ) {
        PDBtrim <- PDB[(real.MD[model.num]+1):(real.MD[model.num+1]-1)]
        write(PDBtrim, file=pdb.filename)
      }
      if( model.num > 0 & model.num == RESULTS$Nmodels ) {
        PDBtrim <- PDB[(real.MD[model.num]+1):(length(PDB)-1)]
        write(PDBtrim, file=pdb.filename)
      }
    }
    #####################################################################
  return(RESULTS)
} 

################################################################################
