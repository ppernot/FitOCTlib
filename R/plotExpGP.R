#  plotExpGP.R
#' Plot outputs from \code{fitExpGP}.
#' @param x a numeric vector
#' @param y a numeric vector of responses/data
#' @param uy a numeric vector of uncertainties
#' @param dataType an numeric (1 or 2) defining the type of data 
#' @param ySmooth a numeric vector of smoothed data
#' @param out a list, output of fitEpxGP
#' @param modScale a real defining the plotting scale for yGP
#' @param nMC an integer number of spaghetti lines
#' @param gPars a list of graphical parameters and colors
#' @param br a list output from \code{printBr()}
#' @return Produces a plot. If \code{gPars$graphTable} is TRUE and
#'         out$prior_PD != 0, the best values of the fit are shown
#'         in tabular form instead of a QQ-plot of the residuals.
#'           
#' @author Pascal PERNOT
#' @export

plotExpGP       <- function(x, y, uy, ySmooth, out,
                            modScale=0.3, nMC=100, gPars,
                            dataType = 2, br=NULL) {
  # Extract graphical parameters
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  
  fit      = out$fit
  method   = out$method
  xGP      = out$xGP
  prior_PD = out$prior_PD
  # lasso    = out$lasso
  
  par(mfrow=c(2-prior_PD,2),pty=pty, mar=mar, mgp=mgp,
      tcl=tcl, lwd=lwd, cex=cex)
  
  if(dataType==1){
    ylabel = "mean amplitude (a.u.)" 
    A0 ='A[0]~(a.u.) '
  }
  if(dataType==2){
    ylabel = "mean intensity (a.u.)" 
    A0='I[0]~(a.u.) '
  }
  
  if(method == 'sample') {
    
    theta   = rstan::extract(fit,'theta')[[1]]
    yGP     = rstan::extract(fit,'yGP')[[1]]
    sigma   = mean(rstan::extract(fit,'sigma')[[1]])
    if(prior_PD == 0) {
      resid   = rstan::extract(fit,'resid')[[1]]
      mod     = rstan::extract(fit,'m')[[1]]
      dL      = rstan::extract(fit,'dL')[[1]]
      lp      = rstan::extract(fit,'lp__')[[1]]
      map     = which.max(lp)
      y_map   = mod[map,]
    }
    
    iMC = sample.int(nrow(theta),nMC)
    
    # Fit
    plot(x,y,pch=20,col=cols[6],
         main= plot_title,
         xlab= xlabel,
         ylab= ylabel)
    grid()
    
    if(prior_PD == 0) {
      # Posterior sample
      
      if(nMC >0) {
        for (i in 1:nMC)
          lines(x, mod[iMC[i],], col=col_tr[4])
      }
      
      # Calculate AVerage Exponential Decay
      mExp = x*0
      for (i in 1:nrow(theta))
        mExp = mExp + expDecayModel(x,theta[i,1:3],dataType)
      mExp = mExp/nrow(theta)
      lines(x,mExp,col=cols[7])
      
      legend('topright', bty='n', cex=cex.leg,
             title = '', title.adj = 1,
             legend=c('data','mean exp. fit','post. sample'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6], cols[7], col_tr2[4])
      )
      legend('topright', bty='n', cex=cex.leg,
             legend=c('','','','','',
                      as.expression(
                        bquote(
                          "br    " == .(formatC(br$br,digits=3)))
                      ),
                      as.expression(
                        bquote(
                          "CI95" == .(paste0(signif(br$CI95,2),collapse='-')))
                      )
             )
      )
    } else {
      # Prior sample
      
      if(nMC >0)
        for (i in 1:nMC)
          lines(x,expDecayModel(x,theta[i,1:3],dataType),col=col_tr[7])
      
      legend('topright', bty='n', cex=cex.leg,
             title = ' ', title.adj = 1,
             legend=c('data','prior sample'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6], col_tr2[7])
      )
    }
    
    
    box()

    # Residuals
    if(prior_PD == 0) {
      res = colMeans(resid)
      ylim=1.2*max(abs(res))*c(-1,1)
      plot(x,res,type='n',
           main ='Residuals',
           xlab = xlabel,
           ylim = ylim, 
           ylab ='residuals (a.u.)' )
      grid(); abline(h=0)
      polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
      points(x,res,pch=20,cex=0.75,col=cols[6])
      lines(x, ySmooth-y_map, col=cols[7])
      legend('topright', bty='n', cex=cex.leg,
             legend=c('mean resid.','data 95% uncert.','smooth - fit'),
             pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
             col=c(cols[6],col_tr2[4],cols[7])
      )
      box()
    }
    
    # Local deviations
    if(prior_PD == 0) {
      plot(x, dL[map,], type = 'n',
           main = expression('Deviation from mean '*L[s]),
           xlab = xlabel,
           ylim = modScale * c(-1,1),
           ylab = 'relative deviation')
      abline(h=0)
      grid()

      if(nMC >0)
        for (i in 1:nMC)
          lines(x, dL[iMC[i],], col=col_tr[4])
      
    } else {
      plot(x, x, type = 'n',
           ylim = modScale*c(-1,1),
           main = expression('Deviation from mean '*L[s]),
           xlab = xlabel,
           ylab = 'relative deviation')
      abline(h=0)
      grid()
    }
    
    Q = t(apply(yGP,2,
                function(x)
                  quantile(x,probs = c(0.025,0.25,0.75,0.975))
    )
    )
    segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])           # 95 %
    segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=3*lwd) # 50 %
    
    legend('topright', bty='n', cex=cex.leg,
           legend=c('50% CI','95% CI','post. sample'),
           pch=NA ,lty=c(1,1,1),lwd=c(3*lwd,lwd,lwd),
           col=c(cols[6],cols[7],col_tr2[4])
    )
    box()
    
  } else {
    
    theta = fit$par$theta
    mod   = fit$par$m
    resid = fit$par$resid
    dL    = fit$par$dL
    yGP   = fit$par$yGP
    sigma = fit$par$sigma
    
    plot(x,y,pch=20,cex=0.5,col=cols[6],
         main= plot_title,
         xlab= xlabel,
         ylab= ylabel)
    grid()
    lines(x,expDecayModel(x,theta,dataType),col=cols[7])
    lines(x,mod, col=cols[4])
    
    legend('topright', bty='n', cex=cex.leg,
           legend=c('data','expo. best fit','best fit'),
           pch=c(20,NA,NA),lty=c(-1,1,1),
           col=c(cols[6],cols[7], col_tr2[4])
    )
    legend(
      'topright',
      bty = 'n',
      cex = cex.leg,
      legend = c(
        '',
        '',
        '',
        '',
        '',
        as.expression(bquote("br    " == .(
          formatC(br$br, digits = 3)
        ))),
        as.expression(bquote("CI95" == .(
          paste0(signif(br$CI95, 2), collapse = '-')
        )))
      )
    )
    box()
    
    # Residus
    ylim=1.2*max(abs(resid))*c(-1,1)
    res = resid
    plot(x,res,type='n',
         main = 'Residuals',
         xlab = xlabel,
         ylim = ylim, 
         ylab = 'residuals (a.u.)')
    grid()
    abline(h=0)
    polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
    points(x,res,pch=20,cex=0.75,col=cols[6])
    lines(x, ySmooth - mod, col=cols[7])
    legend('topright', bty='n', cex=cex.leg,
           legend=c('mean resid.','data 95% uncert.','smooth - fit'),
           pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
           col=c(cols[6],col_tr2[4],cols[7])
    )
    box()
    
    # Local deviations
    plot(x, dL, type = 'l',
         ylim = modScale*c(-1,1),
         col  = cols[4],
         main = expression('Deviation from mean '*L[s]),
         xlab = xlabel,
         ylab = 'relative deviation')
    grid()
    abline(h=0)
    
    if(!is.null(fit$theta_tilde)) {
      S = fit$theta_tilde
      iMC = sample.int(nrow(S),nMC)
      
      c = which(grepl(pattern = 'dL\\[', x=colnames(S)))
      for (i in 1:nMC)
        lines(x, S[iMC[i],c], col=col_tr[4])
      
      c = which(grepl(pattern = 'yGP\\[', x=colnames(S)))
      yGP = S[iMC,c]
      Q = t(apply(yGP,2,
                  function(x)
                    quantile(x,probs = c(0.025,0.25,0.75,0.975))
      )
      )
      segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])       # 95 %
      segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=2*lwd) # 50 %
      legend('topright', bty='n', cex=cex.leg,
             legend=c('50% CI','95% CI','post. sample'),
             pch=NA ,lty=c(1,1,1),lwd=c(2*lwd,lwd,lwd),
             col=c(cols[6],cols[7],col_tr2[4])
      )
      
    } else {
      points(xGP,yGP,pch=19,col=cols[7])
      segments(xGP,yGP,xGP,0*yGP,col=cols[7])
      legend('topright', bty='n', cex=cex.leg,
             legend=c('ctrl points','modulation'),
             pch=c(19,NA) ,lty=c(1,1),lwd=c(-1,lwd),
             col=c(cols[7],cols[4])
      )
    }
    box()
  }
  
  # if(prior_PD == 0) {
  #   # Plot true modulation for synthetic signals
  #   fName = paste0('./Modulation.csv')
  #   if(file.exists(fName)) {
  #     M = read.csv(fName)
  #     lines(M[,1],M[,2],lty=2)
  #   }
  # }
  
  sum = NULL
  if(prior_PD == 0) {
    
    if(method == 'sample') {
      sum  = rstan::summary(fit,
                            use_cache=FALSE,
                            probs=c()
      )$summary[1:3,c(1,3)]
      # fit_summary = fit@.MISC[["summary"]][["msd"]][1:3,1:2]
      # mean = formatC(fit_summary[1:3,1],digits=3, format="g", decimal.mark=".")
      # std  = formatC(fit_summary[1:3,2],digits=1, format="g", decimal.mark=".")
    } else {
      
      pars = c('theta')
      opt = list()
      for (par in pars)
        opt[[par]] = fit$par[[par]]
      opt = unlist(opt,use.names = TRUE)
      
      se = rep(NA, length(opt))
      if(!is.null(fit$hessian)) {
        H = fit$hessian
        tags = colnames(H)
        tags = gsub('\\.','',tags)
        colnames(H) = rownames(H) = tags
        se = list()
        for (par in names(opt))
          se[[par]]  = sqrt(-1/H[par,par])
        se = unlist(se)
      }
      sum = data.frame(mean = opt, sd = se)
    }
    
    if(graphTable) {
      frame()
      # vps <- gridBase::baseViewports() 
      grid::pushViewport(grid::viewport(x=0.74,y=0.28))
      
      parametres = data.frame(row.names=c('C~(a.u.) ',A0,'L[s]~(µm) '),
                              signif(sum[,1],3),signif(sum[,2],2))
      names(parametres) <- c("mean","std")
      
      tt = gridExtra::ttheme_minimal(
        base_size = 56, 
        base_colour = "black", 
        base_family = "", 
        parse = FALSE,
        padding = grid::unit(c(30, 20), "mm"), 
        rowhead = list(fg_params = list(parse=TRUE))
      )
      tbl = gridExtra::tableGrob(parametres, theme=tt)
      separator1 = grid::segmentsGrob(
        x0 = grid::unit(0.1,"npc"), 
        y0 = grid::unit(0,"npc"), 
        x1 = grid::unit(1,"npc"), 
        y1 = grid::unit(0,"npc"),
        gp = grid::gpar(lwd = 4.0)
      )
      separator2 = grid::segmentsGrob(
        x0 = grid::unit(0,"npc"), 
        y0 = grid::unit(0,"npc"), 
        x1 = grid::unit(0,"npc"), 
        y1 = grid::unit(1,"npc"),
        gp = grid::gpar(lwd = 4.0)
      )
      tbl = gtable::gtable_add_grob(
        tbl, 
        grobs = separator1, 
        t = 1, 
        b = 1, 
        l = 1, 
        r = ncol(tbl)
      )
      tbl = gtable::gtable_add_grob(
        tbl, 
        grobs = separator2, 
        t = 1, 
        b = nrow(tbl), 
        l = 2, 
        r = 3
      )
      grid::grid.draw(tbl)
      grid::popViewport()
      
    } else {

      # QQ-plot of weighted residuals
      resw = res/(uy*sigma)
      xlim = range(resw)
      qqnorm(resw, xlim=xlim,ylim=xlim,
             main='Norm. Q-Q plot of wghtd resid.',
             col=cols[6], pch=20)
      abline(a=0,b=1,lty=2)
      grid();box()
      
    }
    
  }
  invisible(
    list(sum = sum)
  )
}
