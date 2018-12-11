#  plotPriPost.R
#' Plot 1 marginal prior and posterior pdf.
#' @param pri prior sample
#' @param pos posterior sample
#' @param tag name of the variable
#' @param xlim an interval
#' @param gPars a list of graphical parameters and colors
#' @param kl a logical to compute and display KL-divergence
#' @return Produces a plot. The KL divergence in bits is
#'         displayed in the title if `kl == TRUE`.
#' @author Pascal PERNOT
#' @export

plotPriPost      <- function(pri,pos,tag,xlim=range(c(pri,pos)),gPars,
                             kl = TRUE) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  # Compute KL-divergence in bits
  if(kl)
    kld = 1.443 * mean(FNN::KL.divergence(pos, pri))
  
  # Plot overlapped densities for 2 samples
  d = density(pri)
  d$y = d$y/max(d$y)
  plot(d, type = 'l', col = cols[4],
       main = ifelse(kl,
                     paste0(tag,' / KL: ',signif(kld,2)),
                     tag),
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(0,1.1), yaxs = 'i')
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
  d = density(pos)
  d$y = d$y/max(d$y)
  lines(d$x,d$y,col=cols[6])
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[6],border=NA)
}
