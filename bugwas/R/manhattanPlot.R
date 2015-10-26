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

# Mahattan plot
.manhattanPlot = function(prefix = NULL,
                         snpPos = NULL,
                         snpType = NULL,
                         pValues = NULL,
                         main = NULL,
                         col = NULL,
                         xlab = "Position in reference genome"){ #c(COL,COL2)
  
  .pl(prefix,{
    
    plot(x = snpPos,
         y = pValues,
         col = col, 
         pch = snpType,
         xlab=xlab, ylab="-log10(p-value)",
         main = main)
  })
  
}