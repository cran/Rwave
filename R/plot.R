
#########################################################################
#    (c) Copyright 1997
#               by
#   Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, Andrea Wang
#              Princeton University
#              All right reserved
########################################################################


plotwt <- function(original, psi, phi, maxresoln, scale=F, yaxtype="s")
#*********************************************************************#
# plotwt
# ------
# Plot wavelet transform
# Output the original signal, wavelet transform Wf starting resoln 1,
# and Sf at the last resoln  
# Used by mw and dw
#*********************************************************************#
{
  par.orig <- par()  # Save the original plot parameters

  par(mfrow = c((maxresoln+2), 1))

  if ( !scale )
  {
    plot.ts(original,bty="n",yaxt=yaxtype)
    for (j in 1:maxresoln)
      plot.ts(psi[,j],bty="n",yaxt=yaxtype)
    plot.ts(phi, bty="n",yaxt=yaxtype)
  }
  else
  {
    limit <- c(min(original, psi, phi), max(original, psi, phi))
    plot.ts(original, ylim=limit, bty="n",yaxt=yaxtype)
    for (j in 1:maxresoln)
      plot.ts(psi[,j], ylim=limit, bty="n",yaxt=yaxtype)
    plot.ts(phi, ylim=limit, bty="n",yaxt=yaxtype)
  }

  par(par.orig)  
}


plot.result <- function(result, original, maxresoln, scale=F, yaxtype="s")
#*********************************************************************#
# Plot.result
# -----------
#    function for the output of the following S functions:
# 	dnpval, dbpval, region, dntrim
#	mnpval,	ext, mntrim, localvar
#*********************************************************************#

{
  par.orig <- par()  # Save the original plot parameters
  par(mar=c(2.0, 2.0, 0.5, 0.5))
  par(oma=c(1,0,0,0)) # Change 1 to a larger number to allow a larger
		      # bottom margin 

  par(mfrow = c(maxresoln+1, 1))

  if ( !scale )
  {
    plot.ts(original, bty="n",yaxt=yaxtype)
    for (j in 1:maxresoln)
      plot.ts(result[,j], bty="n",yaxt=yaxtype)
  }
  else
  {
    ymin <- min(result)
    ymax <- max(result)
    plot.ts(original, bty="n")
    for (j in 1:maxresoln)
      plot.ts(result[,j], ylim=c(ymin, ymax), bty="n",yaxt=yaxtype)
  }
 
  par(par.orig)
  par(mfrow=c(1,1))
}


wpl <- function(dwtrans)
#*********************************************************************#
# wpl
# ---
# Plot dyadic wavelet transform (output of mw).
#*********************************************************************#
{
  maxresoln <- dwtrans$maxresoln
  plotwt(dwtrans$original, dwtrans$Wf, dwtrans$Sf,maxresoln)
  cat("")
}


epl <- function(dwext)
#*********************************************************************#
# epl
# ---
# Plot wavelet transform extrema (output of ext).
#*********************************************************************#
{
  maxresoln <- dwext$maxresoln
  plot.result(dwext$extrema, dwext$original, maxresoln)
  cat("")
}



