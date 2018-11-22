#  estimateExpPrior.R
#' Define/estimate normal multivariate prior pdf for exponential decay parameters.
#' @param x         numeric vector of depths
#' @param uy        numeric vector of uncertainties
#' @param dataType  integer defining the type of data (1:intensity or 
#'                  2:amplitude) )
#' @param priorType string defining the type of prior ('mono' or 'abc')
#' @param out       output list from \code{fitMonoExp}
#' @param ru_theta  optional positive real defining the relative uncertainty 
#'                  on parameters (priorType='mono')
#' @param eps       tolerance parameter for the moments matching method 
#'                  (priorType='abc')
#' @param nb_chains number of MCMC chains (priorType='abc')
#' @param nb_warmup number of warmup steps (priorType='abc')
#' @param nb_iter   number of steps (priorType='abc')
#'  
#' @return A list containing: the center, covariance matrix 
#'         of the prior pdf and a list with the constraint 
#'         and realized statistics.  
#' 
#' @details Provides two ways to buil a normal multivariate prior for 
#'          the exponential decay parameters of the \code{ExpGP} model.
#'          The mean value is in both cases the MAP issued by fitMonoExp.
#'          For the covariance matrix one has two options: 
#'          \describe{
#'             \item{priorType='mono'}{
#'             a covariance matrix is built from the correlation matrix 
#'             estimated by the \code{fitMonoExp} model and a relative 
#'             uncertainty on the parameters (\code{ru_theta}). 
#'             This accounts for the fact that the parameters 
#'             uncertainties provided by fitMonoExp are not reliable,
#'             as they are issued from an invalid model.
#'             }
#'             \item{priorType='abc'}{
#'             the parameters uncertainties/variances are optimized
#'             by a moments matching strategy: (1) the 2-sigma prediction 
#'             uncertainty has to match the 95-th quantile of the absolute
#'             errors of the fitMonoExp model (this statistics is weighted 
#'             by \code{uy}); and (2) the standard deviation of the prediction 
#'             uncertainty has to be as small as possible.
#'             This prior assumes a diagonal covariance matrix.
#'             }
#' } 
#' 
#' @author Pascal PERNOT
#' 
#' @export

estimateExpPrior <- function(x, uy, dataType, priorType = 'mono',
                             out, ru_theta = 0.05, eps = 1e-3,    
                             nb_chains = 4, nb_warmup = 800,
                             nb_iter = nb_warmup + 200) {

  statsObs = function(x){
    # Stats for residuals are Q95 and 0.
    c(quantile(abs(x),probs=c(0.95)), 0.)
  }
  
  theta0    = out$best.theta   # Used by next stage
  
  if(priorType == 'mono') {
    # Scale Monoexp covariance matrix by ru_theta
    cor_theta = out$cor.theta
    u_theta   = ru_theta * theta0
    rList = NULL
    
  } else {
    # ABC nocorr approx.
    resid = out$fit$par$resid
    Sobs  = statsObs(resid/uy)
    
    stanData = list(
      N        = length(x), 
      x        = x, 
      uy       = uy, 
      dataType = dataType,
      Sobs     = Sobs,
      Np       = 3,
      theta    = theta0,
      eps      = eps
      )
    init = list(
      u_theta = ru_theta * theta0
    )
    
    # Sample
    fit = rstan::sampling(
      stanmodels$modUQExp,
      data      = stanData,
      init      = function() {init},
      control   = list(adapt_delta=0.995),
      iter      = nb_iter,
      chains    = nb_chains,
      warmup    = nb_warmup,
      verbose   = FALSE
    )
    
    lp   = rstan::extract(fit,'lp__')[[1]]
    map  = which.max(lp)
    u_theta   = rstan::extract(fit,'u_theta')[[1]][map,]
    cor_theta = diag(c(1,1,1)) # Hyp. nocorr
    rList =     list(
      Sobs   = Sobs,
      Ssim   = rstan::extract(fit,'Ssim')[[1]][map,]
    )
    
  }
  
  # Cov. matrix
  Sigma0 = diag(u_theta) %*% cor_theta %*% diag(u_theta)
  
  return(
    list(
      theta0 = theta0,
      Sigma0 = Sigma0,
      rList  = rList
    )
  )
}
