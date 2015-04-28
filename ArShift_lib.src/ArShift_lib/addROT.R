################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

addROT <- function(Data) {

ROT1 <- ROT2 <- ROT3 <- NULL
 
 for(i in 1:length(Data[,1])) {

   amacid <- as.character(as.vector((Data[i,"AMACID"])))
   if(amacid!="Ala") {
     rt1.ind <- which(Data[i,"CHI1"]>(-120*pi/180) & Data[i,"CHI1"]<=(0*pi/180))
   
     rt2.ind <- which(Data[i,"CHI1"] > (0*pi/180) & Data[i,"CHI1"] <= (120*pi/180))
    
     rt3.ind <- c( which(Data[i,"CHI1"] > (120*pi/180) & Data[i,"CHI1"] <= (180*pi/180)),
                  which(Data[i,"CHI1"] >= (-180*pi/180) & Data[i,"CHI1"] <= (-120*pi/180)) )

     if(length(rt1.ind)==1) { ROT1[i]<-1; ROT2[i]<-0; ROT3[i]<-0 }
     if(length(rt2.ind)==1) { ROT1[i]<-0; ROT2[i]<-1; ROT3[i]<-0 }   
     if(length(rt3.ind)==1) { ROT1[i]<-0; ROT2[i]<-0; ROT3[i]<-1 }     
     if(length(rt1.ind)==0 & length(rt2.ind)==0 & length(rt3.ind)==0){
       print("Warning: Indefinite Chi1 angle!!! from addROT.R !!!")
     }
   } else {
     ROT1[i]<-0; ROT2[i]<-0; ROT3[i]<-0
   } 

 }; rm(i)

 Data <- data.frame(Data, ROT1=ROT1, ROT2=ROT2, ROT3=ROT3)
 return(Data)
}

################################################################################