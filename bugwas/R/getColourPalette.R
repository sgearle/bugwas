# File KmerAnalysisMain.R
# Authors: Earle, S. G., Wu, C.-H. and Wilson, D. J.
#
# Copyright (C) 2015 University of Oxford
#
# This file is part of the bugwas R package.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

# Get a random sample of colours to use for all plots, equal to the number of significant PCs
getColourPalette = function(colourCount = 20, p.pca.bwt, signifCutOff, pc.lim, palette = NULL){
	randomCOL=rep("grey50",colourCount) #maybe fix the 20 colours
	
	if(is.null(palette)){
		palette = getPalette20()
		
	}
	
	if(length(palette) < length(pc.lim)){
		stop("Not enough colours provided in the palette (", length(palette),")for colouring significant PCs (", length(pc.lim),").")
	}
  	
  	if(length(which(p.pca.bwt >= signifCutOff))>0){
    	randomCOL[pc.lim] = palette[1:length(pc.lim)]
  	}
  	
  	return(randomCOL)
	
}