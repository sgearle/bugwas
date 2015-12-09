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

# Get tree patterns
.tree2patterns = function(tree,tiporder=tree$tip.label) {
  n = length(tree$tip.label)
  mtp = matrix(0,n,n+tree$Nnode)
  mtp[,1:n] = diag(n)
  for(i in 1:tree$Nnode) {
    wh = match(extract.clade(tree,n+i)$tip.label,tree$tip.label)
    mtp[wh,n+i] = 1
  }
  mtp = mtp[match(tiporder,tree$tip.label),]
  mtp.f = apply(mtp,2,mean)
  mtp[,mtp.f>0.5] = 1-mtp[,mtp.f>0.5]
  mtp.f = apply(mtp,2,mean)
  rownames(mtp) = tiporder
  colnames(mtp) = paste("Node",1:ncol(mtp),sep="_")
  edge = match(1:ncol(mtp),tree$edge[,2])
  edge.length = tree$edge.length[edge]
  return(list("pat"=mtp,"ancestral_edge"=edge,"ancestral_edge.length"=edge.length))
}