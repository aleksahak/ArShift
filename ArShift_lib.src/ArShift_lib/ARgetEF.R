################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
# MODIFIED for TYR and PHE.
################################################################################
# This function given the index
# of the methyl group in the usual output of the "methyl.name()" function,
# the output of the "methyl.name()" function (ch3.names) and the object
# from the "getarea()" function (areal), returns the electric field pro-
# jection FACTOR along the axis specified in ch3.names. Also requires
# the parsed pdb object.
################################################################################
# VERSION: 28 October 2010, modified to account chain IDs.
# REQUIRES: "get.coord" SOURCE: "get_coord.R"
# REQUIRES: "Xangle" SOURCE: "Xangle.R"  -- the updated version
# REQUIRES: "Xnormcent" SOURCE: "Xnormcent.R"
# REQUIRES: "Vector.translate" SOURCE: "Vector_translate.R"

ARgetEF <- function(resSeqnum, AmAc, query.xyz, whichNMRnucl, areal, parsed.pdb, chainID) {



  ef.app <- query.xyz  ## xyz coordinates of the EF application site

  ## xyz coordinates of the EF destination
  if(whichNMRnucl=="CD1"){
    ef.dir <- get.coord(name="HD1", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="CD2"){
    ef.dir <- get.coord(name="HD2", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="CE1"){
    ef.dir <- get.coord(name="HE1", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="CE2"){
    ef.dir <- get.coord(name="HE2", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="CZ" & AmAc=="TYR"){
    ef.dir <- get.coord(name="OH", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="CZ" & AmAc=="PHE"){
    ef.dir <- get.coord(name="HZ", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
  }
  if(whichNMRnucl=="HD1"){
    ef.dir <- get.coord(name="CD1", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
    ef.dir <- Vector.translate(Vec1B=ef.dir, Vec1E=ef.app, NewB=ef.app) 
  }
  if(whichNMRnucl=="HD2"){
    ef.dir <- get.coord(name="CD2", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
    ef.dir <- Vector.translate(Vec1B=ef.dir, Vec1E=ef.app, NewB=ef.app) 
  }
  if(whichNMRnucl=="HE1"){
    ef.dir <- get.coord(name="CE1", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
    ef.dir <- Vector.translate(Vec1B=ef.dir, Vec1E=ef.app, NewB=ef.app) 
  }
  if(whichNMRnucl=="HE2"){
    ef.dir <- get.coord(name="CE2", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
    ef.dir <- Vector.translate(Vec1B=ef.dir, Vec1E=ef.app, NewB=ef.app) 
  }
  if(whichNMRnucl=="HZ"){
    ef.dir <- get.coord(name="CZ", resSeq=resSeqnum, parsed.pdb=parsed.pdb, chainID=chainID)
    ef.dir <- Vector.translate(Vec1B=ef.dir, Vec1E=ef.app, NewB=ef.app) 
  }

  geo.factor <- NULL
  for(i in 1:length(areal$dist)) {
    q.xyz <- get.coord(areal$name[i], areal$resSeq[i], parsed.pdb, chainID=areal$chain[i])
    cos.angle <- Xangle(Xnormcent(ef.app, q.xyz), Xnormcent(ef.app, ef.dir), cosine=TRUE)
    geo.factor[i] <- (areal$charge[i]*cos.angle)/(areal$dist[i])^2
  }
  EF <- NULL
  EF$CC <- sum(geo.factor)

  # discarding charges of the own ring:
  excl <- which(areal$resName==AmAc     &  # excluding everything own, except the non-ring atoms
                areal$resSeq==resSeqnum & 
                areal$chain==chainID    &
                areal$name!="N"   &
                areal$name!="HN"  &
                areal$name!="CA"  &
                areal$name!="HA"  &
                areal$name!="C"   &
                areal$name!="O"   &
                areal$name!="CB"  &
                areal$name!="HB1" &
                areal$name!="HB2"  )


  EF$FARCC <- sum(geo.factor[-excl])  

  return(EF)
}

################################################################################
