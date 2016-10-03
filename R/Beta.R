####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/

#library(roxygen2)
#library(roxygen) # not working
#roxygenize("SIBER")

library(devtools)
build('SIBER')
install('SIBER')
  
load_all('SIBER')
 

##
detach("package:SIBER", unload=TRUE)
library(SIBER)

    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/SIBER/R/Beta.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/SIBER/R/siberRaw2.R'))

	
}
###################
#
# this is based on direct optimization and thus might be too slow; further improvement would be using EM algorithm
#
###################

# log: whether parameter is parameterized at the log scale
BetaMixutreLogLik<-function(parm=parm,y=y, transform=TRUE){
	a1=parm[1]
	b1=parm[2]
	a2=parm[3]
	b2=parm[4]
	pi=parm[5]
	if(transform){
        a1 <- exp(a1)
        b1 <- exp(b1)
        a2 <- exp(a2)
        b2 <- exp(b2)
        pi <- inv_logit(pi) # logit transform
    }
    #return((-1)*sum(log(dbeta(y,a1,b1)*pi+dbeta(y,a2,b2)*(1-pi))))
	return(sum(log(dbeta(y,a1,b1)*pi+dbeta(y,a2,b2)*(1-pi))))
}

#' fit a 2-component beta mixture
#'
#' can impose the two component to have equal variance (model='E') or unequal variance (model='V')
#'
#' @param y ration values between 0 and 1
#' @param inits Initial value to fit the mixture model. A vector with elements alpha1, beta1, alpha2, beta2 and pi1. 
#' @param model Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. 
#' @return A vector consisting parameter estimates of mu1, mu2, phi1, phi2, pi1, logLik and BIC. For 0-inflated model, mu1=phi1=0. 
#' @export
alphaBeta_to_meanVar <- function(alpha, beta){
    m <- alpha/(alpha+beta)
    v <- alpha*beta/((alpha+beta+1)*(alpha+beta)^2)
    res <- c(m, v)
    names(res) <- c('mean', 'variance')
    res
}
#' convert mean/variance parameterization to Beta(alpha, beta) parameterization for a beta distribution
meanVar_to_alphaBeta <- function(m, v){
    alpha <- (m^2 - m^3 - m*v)/v
    beta <- (m - 2*m^2 + m^3 - v + m*v)/v
    res <- c(alpha, beta)
    names(res) < c('alpha', 'beta')
    res
}

transformAlphaBetaToMuSigma <- function(est){
    alpha1 <- est[1]
    beta1 <- est[2]
    alpha2 <- est[3]
    beta2 <- est[4]
    tt1 <- alphaBeta_to_meanVar(alpha1, beta1)
    tt2 <- alphaBeta_to_meanVar(alpha2, beta2)
    res <- c(tt1[1], tt2[1], sqrt(tt1[2]), sqrt(tt2[2]), est[5])
    names(res) <- c('mean1',  'mean2', 'sigma1','sigma2', 'pi1')
    res
}

transformMuSigmaToAlphaBeta <- function(mu, sigma){
    #browser()
    sigma2 <- sigma^2
    alpha <- (1-mu)*mu*mu/sigma2-mu
    beta <- (1-mu)/mu*alpha
    c(alpha=alpha, beta=beta)
}
#transformMuSigmaToAlphaBeta(mu=bi_thym_beta[ss, 'mu1'], sigma=bi_thym_beta[ss, 'sigma1'])

logit <- function(x) {
    log(x/(1-x))
}
inv_logit <- function(x) {
    exp(x)/(1+exp(x))
}    
#' Fit Beta Mixture Model
#' The function fits a two-component Beta mixture model. 
#' This function directly maximize the log likelihood function through optimization. 
#' @param y A vector representing the RNAseq raw count.
#' @param inits Initial value to fit the mixture model. A vector with elements alpha1, beta1, alpha2, beta2 and pi1. 
#' @param model Character specifying E or V model. E model fits the mixture model with equal dispersion phi while V model doesn't put any constraint. Currently only V model is implemented. 
#' @return A vector consisting parameter estimates of alpha1, beta1, alpha2, beta2 and pi1, logLik and BIC. 
#' @export
#' @examples
#' set.seed(100)
#' y=c(rbeta(100,1,4),rbeta(200,4,1))
#' fitBeta(y=y)
#' SIBER(y, model='Beta')
fitBeta <- function(y,inits=NULL, model='V'){
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
    # data check
    y <- y[!is.na(y)] # remove NA
    if(any(y<0)) stop('Values in beta-mixture cannot be negative!\n')
    if(any(y>1)) stop('Values in beta-mixture cannot be larger than 1!\n')
    if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
    if(model=='E') stop('E model for Beta mixture has not been implemented yet!\n')
    res <- rep(NA, 7)	# res=c(alpha1, beta1, alpha2, beta2 and pi1, logLik, BIC)
    names(res) <- c('alpha1', 'beta1', 'alpha2', 'beta2', 'pi1', 'logLik', 'BIC')
    if(length(y)<10) return(res) # refuse computing when N<10
    nPar <- ifelse(model=='V', 5, 4) # number of parameters
    #browser()
    if(model=='V'){ # V model
		# initial value
		initials <- rep(NA, 5)
		if(is.null(inits)){
			m <- mean(y, na.rm=TRUE)
			v <- var(y, na.rm=TRUE)
			tt <- meanVar_to_alphaBeta(m, v)
            initials[1] <- log(tt[1]*0.8)
            initials[2] <- log(tt[1]*1.2)
            initials[3] <- log(tt[2]*0.8)
            initials[4] <- log(tt[2]*1.2)
            initials[5] <- logit(0.5)
		} else {
			initials <- inits
		}
		#constrL <- c(1e-5, 1e-5, 1e-5, 1e-5, 1e-3) # box constraint
        #constrU <- c(5e2, 5e2, 5e2, 5e2, 1-1e-3)
	}
    #browser()
    optimRes <- try(optim(par=initials, BetaMixutreLogLik, y=y, transform=TRUE, 
            control=list(fnscale=-1, pgtol=1e-8, maxit=1000), method="BFGS"), silent=TRUE) #factr=1e3, 
   if(class(optimRes)!='try-error'){
        tt <- c(optimRes$par, optimRes$value)
        res[1:4] <- exp(tt[1:4])
        res[5] <- inv_logit(tt[5])
        res[6] <- tt[6]
        res[7] <- getBIC(logLik=optimRes$value, nPar=nPar, nObs=length(y)) 
        #names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
    }
    res 
}

#example
#set.seed(100)
#y=c(rbeta(800,1,3),rbeta(200,2,4))
#BetaMixutreLogLik(parm=c(1,2,1,3,0.5),y=y)
#init=c(1,2,1,3,0.5)
#fitBeta(y=y)