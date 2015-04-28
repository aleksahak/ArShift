################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

ArShift.cmd <- function(pdbname="",
                        title="Not specified.",
                        experdata=NULL,
                        examine=0,     #0 for all models, 1,2,3 etc... for specific ones
                        pdbtype="complex",     # "complex", "almost", "medium", "simple"
                        seqshift=0,
                        rereference=FALSE,
                        outorder="bytype",
                        outputname="out.txt",
                        bias = 2   ) {


 
     save(pdbname, title, experdata, examine, pdbtype, seqshift,
          rereference, outorder, outputname, bias, file="ArShift.cmd")
}

################################################################################