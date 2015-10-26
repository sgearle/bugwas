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

###################################################################
## Plot manhattan plot ordered on the x-axis by PCs
## @o: PCs in order of significance
## @which.pc: For each pattern, which PC is it most correlated to
## @pattern: For each variant, what pattern it is
## @p.pca.bwt: Bayesian Wald test results for each PC
## @pc.lim: Which PCs are significant by the BWT
## @lmm: LMM results matrix
## @pat.weight: How many variants does each pattern represent
## @prefix: output prefix
## Outputs:
## Manhattan plot with the x-axis ordered by PCs
###################################################################
.plot_pc_manhattan <- function(o = NULL,
                              which.pc = NULL,
                              pattern = NULL,
                              p.pca.bwt = NULL,
                              pc.lim = NULL,
                              negLog10 = NULL,
                              pat.weight = NULL,
                              prefix = NULL,
                              colourPalette = NULL,
                              npcs = NULL){
  
  # Find which PCs have a variant that is most correlated to it
  m <- match(o,unique(which.pc))
  o.pats <- o[which(is.na(m)==FALSE)]
  # Subset the Bayesian Wald Test -log10(p) to those PCs which have kmers most correlated to them
  p.pca.bwt.pats <- p.pca.bwt[o.pats]
  
  
  num.variants <- sapply(o.pats, function(x, pat.weight=NULL, which.pc = NULL)
    sum(pat.weight[which(which.pc==x)]),pat.weight = pat.weight,
    which.pc = which.pc, USE.NAMES=FALSE)
  
  cols.pcs <- rep_len(c("#5a5a5a","#c6c6c6"), length.out = length(o.pats))
  if(!is.null(pc.lim)){
    cols.pcs[1:length(pc.lim)] <- colourPalette[1:length(pc.lim)]
  }
  m <- match(which.pc,o.pats)
  cols <- cols.pcs[m[pattern]]
  #cols[which(max.cor.pc[which(which.pc==o.pats[i])]<0.3)]="grey50"
  
  
  pos.gap <- sapply(num.variants, function(x) 10000/x, USE.NAMES=FALSE)
  pos.gap[num.variants > 10000] <- 1
  if(!is.null(pc.lim)){
    pos.gap[is.na(match(o.pats, o[pc.lim]))] <- pos.gap[is.na(match(o.pats, o[pc.lim]))]/ (npcs/10)
  } else {
    pos.gap[21:length(pos.gap)] <- pos.gap[21:length(pos.gap)] / (npcs/10)
  }
  
  pos <- rep(0,length(negLog10))
  max.pos <- 0
  plot.lines <- matrix(rep(0, length(o.pats)*2), ncol=2)
  
  which.pc <- which.pc[pattern]
  
  for(i in 1:length(o.pats)){
    
    s = seq(from = (max.pos+pos.gap[i]), by = pos.gap[i], length.out = num.variants[i])
    if(num.variants[i]>1){
      pos[which(which.pc==o.pats[i])] <- sample(s, num.variants[i], replace=FALSE)
    } else {
      pos[which(which.pc==o.pats[i])] <- s
    }
    max.pos <- max(s)
    plot.lines[i,] <- c((max.pos-pos.gap[i]*num.variants[i]+((num.variants[i]*0.15)*pos.gap[i])), (max.pos-((num.variants[i]*0.15)*pos.gap[i])))
    
  }
  
  ymax = max(c(p.pca.bwt.pats, negLog10))
  
  .pl2(paste0(prefix,"_PC_manhattan"),{
    plot(x = pos, y = negLog10, col = cols, xlab = "", ylab = "LMM -log10(p)", ylim = c(0, ymax), xaxt = "n")
    
    for(i in 1:length(o.pats)){
      lines(x = c(plot.lines[i, ]), y = c(p.pca.bwt.pats[i], p.pca.bwt.pats[i]),
            type = "l", col = cols.pcs[i], lwd=2)
      if(i<=20){
        text(x = c(plot.lines[i, 1]+((plot.lines[i, 2]-plot.lines[i, 1])/2)),
             y = c(p.pca.bwt.pats[i]+(max(negLog10)/90)), labels = o.pats[i],
             col = cols.pcs[i], font = 2, cex = 0.9)
      }
    }
  })
  
}