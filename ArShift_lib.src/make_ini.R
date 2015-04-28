################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

compiled = TRUE
if(compiled==TRUE){ suppressPackageStartupMessages(library(compiler)) }
###############################################
# Making sure the existing folders are removed
cur.folders <- dir("..")
delete <- which(cur.folders=="ArShift_exe")
if(length(delete)!=0){ system("rm -r ../ArShift_exe") }
delete <- which(cur.folders=="ArShift_srv")
if(length(delete)!=0){ system("rm -r ../ArShift_srv") }
###############################################
# reading the ArShift.R script:
ArShift.R <- readLines("ArShift_lib/ArShift.R")
system("mkdir ../ArShift_exe")
system("mkdir ../ArShift_exe/ArShift_lib")
write(c("IsThisWebServer = FALSE", ArShift.R), file="../ArShift_exe/ArShift.R")
system("cp ArShift_lib/example.pdb ../ArShift_exe/example.pdb")
system("cp ArShift_lib/example.exp ../ArShift_exe/example.exp")
system("cp ArShift_lib/command.cmd ../ArShift_exe/command.cmd")
###############################################
# Preparing the joint parameter object:
load(paste("ArShift_par/Phe_HD_all.methcs",sep=""))
PHE_HD = lmfit
load(paste("ArShift_par/Phe_HE_all.methcs",sep=""))
PHE_HE = lmfit
load(paste("ArShift_par/Phe_HZ_all.methcs",sep=""))
PHE_HZ = lmfit
load(paste("ArShift_par/Tyr_HD_all.methcs",sep=""))
TYR_HD = lmfit
load(paste("ArShift_par/Tyr_HE_all.methcs",sep=""))
TYR_HE = lmfit
save(PHE_HD, PHE_HE, PHE_HZ, TYR_HD, TYR_HE, 
     file="../ArShift_exe/ArShift_lib/ArShift.par")
rm( list=ls()[which(ls()!="compiled")] )

############################# Internal Flags:
maxnumofentry   = 1000
PATH.LIB        = "ArShift_lib/"
PATH.PAR        = "ArShift_lib/"
# # # # # # 
PHE_HD_mincs    = 5.7; PHE_HD_maxcs    = 8.0
PHE_HE_mincs    = 5.7; PHE_HE_maxcs    = 8.0
PHE_HZ_mincs    = 5.7; PHE_HZ_maxcs    = 8.0
TYR_HD_mincs    = 5.7; TYR_HD_maxcs    = 8.0
TYR_HE_mincs    = 5.7; TYR_HE_maxcs    = 8.0
# # # # # # 
minRerefNum     = 10
###################### End of Internal Flags:

# READING FUNCTIONS:
source("ArShift_lib/ArShift_cmd.R"            )
source("ArShift_lib/isNMRpdb.R"               )
source("ArShift_lib/pdbparser_almost.R"       )
source("ArShift_lib/pdbparser_complex.R"      )
source("ArShift_lib/pdbparser_medium.R"       )
source("ArShift_lib/pdbparser_simple.R"       )
source("ArShift_lib/linesplit.R"              ) 
source("ArShift_lib/get_amac_types.R"         ) 
source("ArShift_lib/pdbProc.R"                ) 
source("ArShift_lib/convert_topology.R"       ) 
source("ArShift_lib/is_integer0.R"            )
source("ArShift_lib/amac_find.R"              )
source("ArShift_lib/get_torsion.R"            )
source("ArShift_lib/get_coord.R"              )
source("ArShift_lib/all_dihedral.R"           )
source("ArShift_lib/topol_retrieve.R"         )
source("ArShift_lib/getarea.R"                )
source("ArShift_lib/Xdist.R"                  )
source("ArShift_lib/Ar_Merge_Same_Type_dist.R")
source("ArShift_lib/ResidueName.R"            )
source("ArShift_lib/normal_vec.R"             )
source("ArShift_lib/Xangle.R"                 )
source("ArShift_lib/Xnormcent.R"              )
source("ArShift_lib/PlaneCurrent.R"           )
source("ArShift_lib/Aniso_Plane.R"            )
source("ArShift_lib/PlaneCurrent_ARG.R"       )
source("ArShift_lib/Aniso_Arg.R"              )
source("ArShift_lib/Aniso_Pept.R"             )
source("ArShift_lib/Xnormcent_length.R"       )
source("ArShift_lib/Vector_translate.R"       )
source("ArShift_lib/vec2_aver.R"              )
source("ArShift_lib/RingCurrent.R"            )
source("ArShift_lib/Ar_Ring_AmAc.R"           )
source("ArShift_lib/ARgetEF.R"                )
source("ArShift_lib/Ar_get_single_dframe.R"   )
source("ArShift_lib/average_pars.R"           )
source("ArShift_lib/addCos.R"                 )
source("ArShift_lib/addMergeDist.R"           )
source("ArShift_lib/addROT.R"                 )
source("ArShift_lib/getP.R"                   )
source("ArShift_lib/uniqizer.R"               )
# BIT-COMPILING THE FUNCTIONS:
if(compiled==TRUE){
  is.integer0        <- cmpfun(is.integer0, options=list(suppressUndefined=TRUE))
  linesplit          <- cmpfun(linesplit,   options=list(suppressUndefined=TRUE))
  Xangle             <- cmpfun(Xangle,      options=list(suppressUndefined=TRUE))
  Xnormcent          <- cmpfun(Xnormcent,   options=list(suppressUndefined=TRUE))
  Vector.translate   <- cmpfun(Vector.translate, options=list(suppressUndefined=TRUE))
  Xdist              <- cmpfun(Xdist,       options=list(suppressUndefined=TRUE))
  uniqizer           <- cmpfun(uniqizer, options=list(suppressUndefined=TRUE))
  get.coord          <- cmpfun(get.coord,   options=list(suppressUndefined=TRUE))
  get.torsion        <- cmpfun(get.torsion, options=list(suppressUndefined=TRUE))
  convert.topology   <- cmpfun(convert.topology, options=list(suppressUndefined=TRUE))
  topol.retrieve     <- cmpfun(topol.retrieve, options=list(suppressUndefined=TRUE))
  Ar.Merge.Same.Type.dist <- cmpfun(Ar.Merge.Same.Type.dist, options=list(suppressUndefined=TRUE))
  Ar.Ring.AmAc       <- cmpfun(Ar.Ring.AmAc, options=list(suppressUndefined=TRUE))
  ARgetEF            <- cmpfun(ARgetEF, options=list(suppressUndefined=TRUE))
  getarea            <- cmpfun(getarea, options=list(suppressUndefined=TRUE))
  pdbparser.complex  <- cmpfun(pdbparser.complex, options=list(suppressUndefined=TRUE))
  pdbparser.almost   <- cmpfun(pdbparser.almost, options=list(suppressUndefined=TRUE))
  pdbparser.medium   <- cmpfun(pdbparser.medium, options=list(suppressUndefined=TRUE))
  pdbparser.simple   <- cmpfun(pdbparser.simple, options=list(suppressUndefined=TRUE))
  pdbProc            <- cmpfun(pdbProc, options=list(suppressUndefined=TRUE))
  get.amac.types     <- cmpfun(get.amac.types, options=list(suppressUndefined=TRUE))
  normal.vec         <- cmpfun(normal.vec, options=list(suppressUndefined=TRUE))
  vec2.aver          <- cmpfun(vec2.aver, options=list(suppressUndefined=TRUE))
  Xnormcent.length   <- cmpfun(Xnormcent.length, options=list(suppressUndefined=TRUE))
  addCos             <- cmpfun(addCos, options=list(suppressUndefined=TRUE))
  addMergeDist       <- cmpfun(addMergeDist, options=list(suppressUndefined=TRUE))
  average.pars       <- cmpfun(average.pars, options=list(suppressUndefined=TRUE))
  addROT             <- cmpfun(addROT, options=list(suppressUndefined=TRUE))
  getP               <- cmpfun(getP, options=list(suppressUndefined=TRUE))
  amac.find          <- cmpfun(amac.find, options=list(suppressUndefined=TRUE))
  all.dihedral       <- cmpfun(all.dihedral, options=list(suppressUndefined=TRUE))
  RingCurrent        <- cmpfun(RingCurrent, options=list(suppressUndefined=TRUE))
  Aniso.Plane        <- cmpfun(Aniso.Plane, options=list(suppressUndefined=TRUE))
  Aniso.Pept         <- cmpfun(Aniso.Pept, options=list(suppressUndefined=TRUE))
  Aniso.Arg          <- cmpfun(Aniso.Arg, options=list(suppressUndefined=TRUE))
  PlaneCurrent.ARG   <- cmpfun(PlaneCurrent.ARG, options=list(suppressUndefined=TRUE))
  PlaneCurrent       <- cmpfun(PlaneCurrent, options=list(suppressUndefined=TRUE))
  ResidueName        <- cmpfun(ResidueName, options=list(suppressUndefined=TRUE))
  isNMRpdb           <- cmpfun(isNMRpdb, options=list(suppressUndefined=TRUE))
  Ar.get.single.dframe <- cmpfun(Ar.get.single.dframe, options=list(suppressUndefined=TRUE))
  ArShift.cmd        <- cmpfun(ArShift.cmd, options=list(suppressUndefined=TRUE))
}

load("ArShift_par/er_data.Rdata")             # object for P calculations
readlib.dih <- readLines("ArShift_lib/dihedral.lib")

# Processing the topology file for further accelerated usage:
proc.topol <- readLines( "ArShift_lib/topology.lib" )
proc.topol <- proc.topol[-grep("#", proc.topol)]                
proc.topol <- lapply(proc.topol, linesplit)                    
# # #

# Utility function for P calculations
fun5 <- function(i){
  if(i=="xxxxxx"){
    return(NA)
  } else {
    return(as.numeric(i))
  }
}
# # # # # # # # # # # # # # # # # # # 

save(list=ls(), file="../ArShift_exe/ArShift_lib/ArShift.ini")
system("cp ArShift_manual/ArShift_manual.pdf ../ArShift_exe/README.pdf")

startdir <- getwd()
setwd("..")
system("tar czvf ArShift_exe.tar.gz ArShift_exe")
setwd(startdir)

# # # # # Generating the Web Server version of ArShift:
# reading the ArShift.R script:
ArShift.R <- readLines("ArShift_lib/ArShift.R")
system("mkdir ../ArShift_srv")
system("cp -r ArShift_srv_lib/ArShift ../ArShift_srv")
system("cp -r ../ArShift_exe/ArShift_lib ../ArShift_srv/ArShift/WEB-INF/ArShift_lib")
system("cp -r ArShift_manual/ArShift_manual.pdf ../ArShift_srv/ArShift/ArShift_manual.pdf")
system("mv  ../ArShift_exe.tar.gz ../ArShift_srv/ArShift/ArShift_exe.tar.gz")
write(c("IsThisWebServer = TRUE", ArShift.R), file="../ArShift_srv/ArShift/WEB-INF/ArShift.R")
setwd("../ArShift_srv/ArShift")
system("jar cf ArShift.war . ")
system("mv ArShift.war ..")
setwd("..")
system("rm -r ArShift")
setwd("../ArShift_lib.src")

print("DONE")
