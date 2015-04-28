################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

addCos <- function(Data, ang="PHI") {

 COS1 <- COS2 <- COS3 <- COS4 <- COS5 <- COS6 <- COS7 <- COS8 <- COS9 <- COS10 <- NULL
  
 for(i in 1:length(Data[,1])) {
  if(Data[i, ang]!=0) {  
   COS1[i]  <- cos( Data[i, ang] )                         
   COS2[i]  <- cos( 5*Data[i, ang] )                       
   COS3[i]  <- cos( 3*Data[i, ang] )                       
   COS4[i]  <- (cos( Data[i, ang] ))^2                     
   COS5[i]  <- cos( 2*Data[i, ang] )^2                    
   COS6[i]  <- cos( 3*Data[i, ang] )^2      
   COS7[i]  <- cos( Data[i, ang] + pi/2)     
   COS8[i]  <- cos( Data[i, ang] + pi/2)^3       
   COS9[i]  <- sin( Data[i, ang] )*cos( Data[i, ang] ) 
   COS10[i] <- cos( 2*(Data[i, ang] + (30*pi/180)) ) 
  } else {
    COS1[i] <- COS2[i] <- COS3[i] <- COS4[i] <- COS5[i] <- COS6[i] <- COS7[i] <- COS8[i] <- COS9[i] <- COS10[i] <- 0
  }  
 }; rm(i)

 text <- paste("Data <- data.frame(Data, ",ang,"COS1=COS1, ",ang,"COS2=COS2, ",ang,"COS3=COS3, ",
                                           ang,"COS4=COS4, ",ang,"COS5=COS5, ",ang,"COS6=COS6, ",
                                           ang,"COS7=COS7, ",ang,"COS8=COS8, ",ang,"COS9=COS9, ",
                                           ang,"COS10=COS10)",sep="")
  eval(parse(text=text))
  return(Data)
}

################################################################################