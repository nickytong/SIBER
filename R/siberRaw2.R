#######################################
### raw code to build SIBER package
#######################################
###
### NEW in 0.9.3: update fitLN such that 
###

### (1) Explicitely deal with 0-inflation in NB and GP models
###     Rather than using BIC, we detect 0-inflated genes by empirical approach based on observed percentage of 0s.

### (2) actually 3 models are fitted: E, V and 0-inflated. When zeroPercentThr is not achieved, 0-inflated model is inactive. However,
###     when this is achieved, E and V are deactivated while 0-inflated model is chosen.

### (3) For output, all results return c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') or 
###     c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
###     For 0-inflated models, mu1=phi1=0 or mu1=sigma1=0

### (4) note for the parameter inits:
###     This parameter is assumed to be a vector of length 5 corresponding to c('mu1', 'mu2', 'phi1', 'phi2', 'pi1'). Note LN model needs no initial values.
###     For E model, phi1=phi2
###     For 0-inflated model, only mu2, phi2, pi1 is used to initialize. The other elements can be NA or anything else which are not used.

## compute BIC
getBIC <- function(logLik, nPar, nObs){
	-2*logLik+nPar*log(nObs)
}
# convert mean mu and variance Var to overdispersion phi by: Var=mu+phi*mu^2 (NB) or Var=mu phi (GP)
# modified on 07/31/12: add scenarios where dist=LN,Normal or Gaussian: just return sqrt(Var) since we only focus on parameterization on log scale in LNmix or normalMix
muVarToPhi <- function(mu, Var, dist='NB'){
	if(dist=='NB') {
	res <- (Var-mu)/mu^2
	res[which(res<0)] <- 0 # under-dispersion, really rare
	
	}
	if(dist=='GP'){
	res <- Var/mu
	}
	if(dist=='LN'|dist=='Gaussian'|dist=='Normal'){ # no action on normal or LN since we've only focus on parameterization on log scale in LNmix or normalMix: just return SD
	res <- sqrt(Var) # returned phi is sd itself
	}
	
	res
}
# given mu and phi, calculates variance. i.e. in NB(mu, phi), compute its variance as mu+phi*mu^2; for GP, var=phi*mu
muPhiToVar <- function(mu, phi, dist='NB'){
	if(dist=='NB') {
	res <- mu+phi*mu^2
	} 
	if(dist=='GP'){
	res <- phi*mu
	}
	res
}

# given mu and phi, calculates variance. i.e. in NB(mu, phi), compute its variance as mu+phi*mu^2; for GP, var=phi*mu
muPhiToVar <- function(mu, phi, dist='NB'){
	if(dist=='NB') {
	res <- mu+phi*mu^2
	} 
	if(dist=='GP'){
	res <- phi*mu
	}
	res
}

#######################################################################################################
############################################# NB mix ##################################################
#######################################################################################################

logLikNegBinOneCompByOverdisp <- function(mu=100, phi=10, y, d){
#browser()
	lgamma(1/phi+y)-lgamma(y+1)-lgamma(1/phi)-
		1/phi*log(phi*mu*d+1)+y*log(phi*mu*d)-y*log(phi*mu*d+1)
}
##### check:
## logLikNegBinOneCompByOverdisp(mu=200, phi=100, datasetsGP_H1[[1]][1, ], d=d)

logLikNegBinMixByOverdisp <- function(par, y=y, d, model='V'){
	if(model=='V'){
	mu1 <- par[1]; mu2 <- par[2]
	phi1 <- par[3]; phi2 <- par[4]
	pi1 <- par[5]
	res <- sum(log(pi1*exp(logLikNegBinOneCompByOverdisp(mu1, phi1, y, d))+
				(1-pi1)*exp(logLikNegBinOneCompByOverdisp(mu2, phi2, y, d))))
	
	} else { # for the E model, only 4 parameters
	mu1 <- par[1]; mu2 <- par[2]
	phi <- par[3];
	pi1 <- par[4]	
	res <- sum(log(pi1*exp(logLikNegBinOneCompByOverdisp(mu1, phi, y, d))+
				(1-pi1)*exp(logLikNegBinOneCompByOverdisp(mu2, phi, y, d))))
	}
#browser()	
	if(is.infinite(res)){
		res <- -1e100 # very important
	}
	res
}

## logLik for 0-inflated model: 3 parameters: mu2, phi2, pi1
logLikNegBin0inflByOverdisp <- function(par, y=y, d){
	mu2 <- par[1]
	phi2 <- par[2]
	pi1 <- par[3]
	res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*exp(logLikNegBinOneCompByOverdisp(mu2, phi2, y, d))))
#browser()	
	if(is.infinite(res)){
		res <- -1e100 # very important
	}
	res
}

#' Fit Negative Binomial Mixture Model
#' The function fits a two-component Negative Binomial mixture model. 
#'
#' This function directly maximize the log likelihood function through optimization. 
#' With this function, three models can be fitted: (1) negative binomial mixture with equal dispersion (E model); 
#' (2) negative binomial mixture with unequal dispersion (V model); (3) 0-inflated negative binomial model.
#' The 0-inflated negative binomial has the following density function:

#' \eqn{P(Y=y)=\pi D(y) + (1-\pi)NB(\mu, \phi)} where D is the point mass at 0 while \eqn{NB(\mu, \phi)} is the density
#' of negative binomial distribution with mean \eqn{\mu} and dispersion \eqn{\phi}. The variance is \eqn{\mu+\phi \mu^2}.

#' The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
#' rule overrides the model argument when observed percentae of zero count exceeds the threshold.
#' 
#' @param y A vector representing the RNAseq raw count.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data. For the LN model, the original data would be devided by this vector. 
#' @param inits Initial value to fit the mixture model. A vector with elements mu1, mu2, phi1, phi2 and pi1. For 0-inflated model,
#' only mu2, phi2, pi1 are used while the other elements can be arbitrary.
#' @param model Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. 
#' @param zeroPercentThr A scalar specifying the minimum percent of zero counts needed when fitting a zero-inflated Negative Binomial model. This
#' parameter is used to deal with zero-inflation in RNAseq count data. When the percent of zero exceeds this threshold, 
#' rather than fitting a 2-component negative binomial mixture, a mixture of point mass at 0 and negative binomial is fitted. 
#' @return A vector consisting parameter estimates of mu1, mu2, sigma1, sigma2, pi1, logLik and BIC. For 0-inflated model, mu1=sigma1=0.
#' @examples
#' set.seed(1000)
#' dat <- rnbinom(100, mu=1000, size=1/0.2)
#' fitNB(y=dat)
#' @export
fitNB <- function(y, d=NULL, inits=NULL, model='V', zeroPercentThr=0.2){
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
	percent0 <- mean(y==0, na.rm=TRUE) 
	res <- rep(NA, 7)	# res=c(mu1, mu2, phi1, phi2, pi1, logLik, BIC)
	names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
	nPar <- ifelse(model=='V', 5, 4) # number of parameters
	
	if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
	
	# Case I: no severe 0-inflation, fit 2-comp E or V models
	if(percent0<=zeroPercentThr){ 
		if(model=='V'){ # V model
			# initial value
			initials <- rep(NA, 5)
			if(is.null(inits)){
				temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
				initials[1] <- temp[1]
				initials[2] <- temp[2]
				overallPhi <- muVarToPhi(mean(y), var(y), dist='NB')
				initials[3] <- overallPhi*0.5
				initials[4] <- overallPhi*2
				initials[5] <- 0.5
			} else {
				initials <- inits
			}
			constrL <- c(1e-5, 1e-5, 1e-5, 1e-5, 0) # box constraint. phi not set to exactly 0, mu not 0. otherwise, pdf not defined
			constrU <- c(5e7, 5e7, 1e3, 1e3, 1)
		} else { # E model
			# initial value
			initials <- rep(NA, 4)
			if(is.null(inits)){
				temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
				initials[1] <- temp[1]
				initials[2] <- temp[2]
				overallPhi <- muVarToPhi(mean(y), var(y), dist='NB')
				initials[3] <- overallPhi
				initials[4] <- 0.5
			} else {
				if(length(inits)==5) { # in case mu1, mu2, phi, phi, pi1 is specified, we only need 4 elements
					initials <- inits[c(1:3, 5)]
				} else {
				initials <- inits # assumed 4 elements
				}
			}
			constrL <- c(0, 0, 1e-5, 0) # box constraint
			constrU <- c(5e7, 5e7, 1e3, 1)
		}
		#browser()	
		optimRes <- try(optim(par=initials, logLikNegBinMixByOverdisp, y=y, d=d, model=model, lower=constrL, 
			upper=constrU, control=list(fnscale=-1, pgtol=1e-16, factr=1e3, maxit=3000), method="L-BFGS-B"), silent=TRUE)
		if(class(optimRes)!='try-error'){
			if(model=='V'){
				res <- c(optimRes$par, optimRes$value)
			} else { # E model, formulate result such that it is like mu1, mu2, phi1, phi2, pi1, logLik (phi1=phi2)
				temp <- optimRes$par
				res <- c(temp[1:3], temp[3], temp[4], optimRes$value)
			}
			res[7] <- getBIC(logLik=optimRes$value, nPar=nPar, nObs=length(y)) 
		}
	} else {
	# Case II: 0-inflation, override E or V models	
		# initial value
		initials <- rep(NA, 3)
		if(is.null(inits)){
			initials[1] <- mean(y[y!=0], na.rm=TRUE)
			initials[2] <- muVarToPhi(mean(y[y!=0]), var(y[y!=0]), dist='NB')
			initials[3] <- 0.5
		} else {
			initials <- inits[, c(2, 4, 5)] # pick mu2, phi2, pi1 as initials
		}
		constrL <- c(1e-5, 1e-5, 0) 
		constrU <- c(5e7, 1e3, 1)
		optimRes <- try(optim(par=initials, logLikNegBin0inflByOverdisp, y=y, d=d, lower=constrL, 
			upper=constrU, control=list(fnscale=-1, pgtol=1e-16, factr=1e3, maxit=3000), method="L-BFGS-B"), silent=TRUE)
		if(class(optimRes)!='try-error'){
			temp <- optimRes$par
			res[1:6] <- c(0, temp[1], 0, temp[2], temp[3], optimRes$value)
			res[7] <- getBIC(logLik=optimRes$value, nPar=3, nObs=length(y)) 
		}
	}	
	res
}
## check
## fitNB(y=dat_2comp, model='V')
## fitNB(y=dat_2comp, model='E')





#######################################################################################################
############################################# GP mix ##################################################
#######################################################################################################
## logLik for GenPoi parameterized by overdispersion and mean. this is a 
## temporary function since it is only for one component
# phi: needs phi>1
# d: normalization constant, a vector of length(y)
# returns a vector of logLik of length(y)

logLikGenPoiOneCompByOverdisp <- function(mu=100, phi=10, y, d){
#browser()
	log(mu*d/sqrt(phi))+(y-1)*log(mu*d/sqrt(phi)+(1-1/sqrt(phi))*y)-mu*d/sqrt(phi)-(1-1/sqrt(phi))*y-lfactorial(y)
}
##### check:
## logLikGenPoiOneCompByOverdisp(mu=200, phi=100, datasetsGP_H1[[1]][1, ], d)
## dnbinom(dat_2comp, mu=4000, size=4000/9, log=TRUE)


## logLik for GenPoi mixture distribution parameterized by overdispersion and V model
# par: in the order of (mu1, mu2, phi1, phi2, pi1)
# returns a scalar
logLikGenPoiMixByOverdisp <- function(par, y=y, d, model='V'){
	if(model=='V'){
	mu1 <- par[1]; mu2 <- par[2]
	phi1 <- par[3]; phi2 <- par[4]
	pi1 <- par[5]
	res <- sum(log(pi1*exp(logLikGenPoiOneCompByOverdisp(mu1, phi1, y, d))+
				(1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi2, y, d))))
	
	} else { # for the E model, only 4 parameters
	mu1 <- par[1]; mu2 <- par[2]
	phi <- par[3];
	pi1 <- par[4]	
	res <- sum(log(pi1*exp(logLikGenPoiOneCompByOverdisp(mu1, phi, y, d))+
				(1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi, y, d))))
	}
#browser()	
	if(is.infinite(res)){
		res <- -1e100 # very important
	}
	res
}
### check
#par <- c(1000, 4000, 10, 10, 0.5); logLikGenPoiMixByOverdisp(par, dat_2comp, model='V')

logLikGenPoi0inflByOverdisp <- function(par, y=y, d){
	mu2 <- par[1]
	phi2 <- par[2]
	pi1 <- par[3]
	res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*exp(logLikGenPoiOneCompByOverdisp(mu2, phi2, y, d))))
#browser()	
	if(is.infinite(res)){
		res <- -1e100 # very important
	}
	res
}


###################
### fitting GP mixture of 2-component by direct optimization. E and V model are available parameterized by over-dispersion
###################
# suppressWarn: whether to suppress warning that comes from optim() convergence criterion
#' Fit Generalized Poisson Mixture Model
#' The function fits a two-component Generalized Poisson mixture model. 
#'
#' This function directly maximize the log likelihood function through optimization. 
#' With this function, three models can be fitted: (1) Generalized Poisson mixture with equal dispersion (E model); 
#' (2) Generalized Poisson mixture with unequal dispersion (V model); (3) 0-inflated Generalized Poisson model.
#' The 0-inflated Generalized Poisson has the following density function:
#' 
#' \eqn{P(Y=y)=\pi D(y) + (1-\pi)GP(\mu, \phi)} where D is the point mass at 0 while \eqn{GP(\mu, \phi)} is the density
#' of Generalized Poisson distribution with mean \eqn{\mu} and dispersion \eqn{\phi}. The variance is \eqn{\phi \mu}.
#' 
#' The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
#' rule overrides the model argument when observed percentae of zero count exceeds the threshold.
#' @param y A vector representing the RNAseq raw count.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data.
#' @param inits Initial value to fit the mixture model. A vector with elements mu1, mu2, phi1, phi2 and pi1.
#' @param model Character specifying E or V model. E model fits the mixture model with equal dispersion phi while V model doesn't put any constraint.
#' @param zeroPercentThr}{A scalar specifying the minimum percent of zero counts needed when fitting a zero-inflated 
#' Generalized Poisson model. This parameter is used to deal with zero-inflation in RNAseq count data. When the percent of zero exceeds this threshold, 
#' rather than fitting a 2-component Generalized Poisson mixture, a mixture of point mass at 0 
#' and Generalized Poisson is fitted. 
#' @return A vector consisting parameter estimates of mu1, mu2, phi1, phi2, pi1, logLik and BIC. 
#' For 0-inflated model, mu1=phi1=0.
#' @export
fitGP <- function(y, d=NULL, inits=NULL, model='V', zeroPercentThr=0.2){
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
	percent0 <- mean(y==0, na.rm=TRUE) 
	
	res <- rep(NA, 7)	# res=c(mu1, mu2, phi1, phi2, pi1, logLik, BIC)
	names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
	nPar <- ifelse(model=='V', 5, 4) # number of parameters
	
	if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
	
	# Case I: no severe 0-inflation, fit 2-comp E or V models
	if(percent0<=zeroPercentThr){ 
		if(model=='V'){ # V model
			# initial value
			initials <- rep(NA, 5)
			if(is.null(inits)){
				temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
				initials[1] <- temp[1]
				initials[2] <- temp[2]
				overallPhi <- muVarToPhi(mean(y), var(y), dist='GP')
				initials[3] <- overallPhi*0.5
				initials[4] <- overallPhi*2
				initials[5] <- 0.5
			} else {
				initials <- inits
			}
			constrL <- c(1e-5, 1e-5, 1, 1, 0) # box constraint
			constrU <- c(5e7, 5e7, 1e10, 1e10, 1)
		} else { # E model
			# initial value
			initials <- rep(NA, 4)
			if(is.null(inits)){
				temp <- quantile(y, prob=c(1/8, 7/8), na.rm=T)
				initials[1] <- temp[1]
				initials[2] <- temp[2]
				overallPhi <- muVarToPhi(mean(y), var(y), dist='GP')
				initials[3] <- overallPhi
				initials[4] <- 0.5
			} else {
				if(length(inits)==5) { # in case mu1, mu2, phi, phi, pi1 is specified, we only need 4 elements
				initials <- inits[c(1:3, 4)]
				} else {
				initials <- inits
				}
			}
			constrL <- c(1e-5, 1e-5, 1, 0) # box constraint
			constrU <- c(5e7, 5e7, 1e10, 1)
		}
#browser()	
		optimRes <- try(optim(par=initials, logLikGenPoiMixByOverdisp, y=y, d=d, model=model, lower=constrL, 
			upper=constrU, control=list(fnscale=-1, pgtol=1e-16, factr=1e3, maxit=3000), method="L-BFGS-B"), silent=TRUE)
		if(class(optimRes)!='try-error'){
			if(model=='V'){
				res <- c(optimRes$par, optimRes$value)
			} else { # E model, formulate result such that it is like mu1, mu2, phi1, phi2, pi1, logLik (phi1=phi2)
				temp <- optimRes$par
				res <- c(temp[1:3], temp[3], temp[4], optimRes$value)
			}
			
			#res$conv <- optimRes$convergence
			res[7] <- getBIC(logLik=optimRes$value, nPar=nPar, nObs=length(y)) 
			names(res) <- c('mu1', 'mu2', 'phi1', 'phi2', 'pi1', 'logLik', 'BIC')
		}
	} else {
	# Case II: 0-inflation, override E or V models	
		# initial value
		initials <- rep(NA, 3)
		if(is.null(inits)){
			initials[1] <- mean(y[y!=0], na.rm=TRUE)
			initials[2] <- muVarToPhi(mean(y[y!=0]), var(y[y!=0]), dist='GP')
			initials[3] <- 0.5
		} else {
			initials <- inits[, c(2, 4, 5)] # pick mu2, phi2, pi1 as initials
		}
		constrL <- c(1e-5, 1, 0) 
		constrU <- c(5e7, 1e10, 1)
		optimRes <- try(optim(par=initials, logLikGenPoi0inflByOverdisp, y=y, d=d, lower=constrL, 
			upper=constrU, control=list(fnscale=-1, pgtol=1e-16, factr=1e3, maxit=3000), method="L-BFGS-B"), silent=TRUE)
		if(class(optimRes)!='try-error'){
			temp <- optimRes$par
			res[1:6] <- c(0, temp[1], 0, temp[2], temp[3], optimRes$value)
			res[7] <- getBIC(logLik=optimRes$value, nPar=3, nObs=length(y)) 
		}
	}
	res
}
## check
## fitGP(y=dat_2comp, model='V')
## fitGP(y=dat_2comp, model='E')



#######################################################################################################
############################################# LN mix ##################################################
#######################################################################################################
# d: normalization constant. The definition is different from TMM or RLE. Our d is the scale factor from
#    true expression level to observed count. Hence, Cs~d*mu_{c(s)}. However, TMM is the scale factor from
#    observed count to true expression level, which is the opposite direction: Cs*TMM~mu_{c(s)}. Thus, d=1/TMM
#    For LN model, we need fit normal mixture on log((Cs+eps)/d)=log(Cs+eps)-log(d)
##### capable dealing with 0-inflation

#' Fit Log Normal Mixture Model
#' The function fits a two-component log normal mixture model. 
#'
#' The parameter estimates from log normal mixture is obtained by taking logarithm and fit normal mixture. We use
#' mclust package to obtain parameter estimates of normal mixture model. In particular, \eqn{log_{base}(\frac{y+eps}{d})} is used to
#' fit to normal mixture model. 
#' 
#' With this function, three models can be fitted: (1) log normal mixture with equal variance (E model); 
#' (2) Generalized Poisson mixture with unequal variance (V model); (3) 0-inflated log normal model.
#' The 0-inflated log normal has the following density function:
#' 
#' \eqn{P(Y=y)=\pi D(y) + (1-\pi)LN(\mu, \sigma)} where D is the point mass at 0 while \eqn{LN(\mu, \sigma)} is the density
#' of log normal distribution with mean \eqn{\mu} and variance \eqn{\sigma^2}. 
#' 
#' The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
#' rule overrides the model argument (E or V) when observed percentae of zero count exceeds the threshold.
#' 
#' @param y A vector representing the RNAseq raw count.
#' @param base The logarithm base defining the parameter estimates in the logarithm scale. This is also the base of log transformation applied to the data. 
#' @param eps A scalar to be added to the count data to avoid taking logarithm of zero.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data. For the LN model, the original data would be devided by this vector. 
#' @param model Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. 
#' @param zeroPercentThr}{A scalar specifying the minimum percent of zero counts needed when fitting a zero-inflated 
#' log normal model. This
#' parameter is used to deal with zero-inflation in RNAseq count data. When the percent of zero exceeds this threshold, 
#' 1-comp mixture LN model is used to estimate mu and sigma from nonzero count.
#' @param logLikToLN logical indicating if the log likelihod is defined on the transformed value or the orginal value from log noral distribution. 
#' @return A vector consisting parameter estimates of mu1, mu2, sigma1, sigma2, pi1, logLik and BIC. For 0-inflated model, mu1=sigma1=0.
#' @export
fitLN <- function(y, base=10, eps=10, d=NULL, model='E', zeroPercentThr=0.2, logLikToLN=TRUE){
	if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
	percent0 <- mean(y==0, na.rm=TRUE)
	res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
	names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
#browser()
	if(percent0<=zeroPercentThr){ # no severe 0-inflation, fit 2-comp
		Dat <- (log(y+eps)-log(d))/log(base) # transform data; change of base formula: log_b(x)=log_e(x)/log_e(b)
		mc <- try(Mclust(Dat, G = 2, modelNames = model), silent=TRUE)
		res[1:7] <- extractMclustPar(mc, modelName=model, logLikToLN=logLikToLN, dat=Dat)
	} else { # severe 0-inflation. 1-comp model. need control mu and sd especially for RPKM
		nonzero <- y[y!=0]
		#if(mean(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001, force it at least 1 to avoid FP
		#### median is more robust: not susceptible to single outlier expression
		#### bug here: 10/23/2012: when 0-inflated with small RPKM, i.e. 0.002, the estimated mu2=eps(1) due to this artifact
		#if(median(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001 (mu2=-5),  force it at least 1 to avoid FP
		Dat <- (log(nonzero)-log(d[y!=0]))/log(base) ###### no need to add eps since no 0 now!
		#mc_1comp <- try(Mclust(Dat, G = 1), silent=TRUE)
		#Modified to MLE such that c(00000, 11111) wouldn't fail
		mu2 <- mean(Dat)
		if(mu2<0 & median(nonzero)<1) mu2 <- 0 # mu2 needs to be larger than mu1 when 0-inflated
		sigma2 <- max(sd(Dat), 1e-2) ## in case sd=0.0004
		pi1 <- percent0 # 
		res[1:5] <- c(0, mu2, 0, sigma2, pi1) # mu1=0, sigma1=0, logLik=BIC=NA
		res[6] <- logLik0inflatedLN(y, mu=mu2, sigma=sigma2, pi1=pi1, logLikToLN=logLikToLN)
		res[7] <- getBIC(logLik=res[6], nPar=3, nObs=length(y)) 
	}
	res
}

# calculate log-lik for 0-inflated normal or log normal model. log normal model is achieved by setting logLikToLN=TRUE
logLik0inflatedLN <- function(y, mu, sigma, pi1, logLikToLN=TRUE){
	if(logLikToLN==TRUE){
		res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*dlnorm(y, meanlog=mu, sdlog=sigma)))
	} else {
		res <- sum(log(pi1*as.numeric(y==0)+(1-pi1)*dnorm(y, mean=mu, sd=sigma)))
	}
	res
}
#######################################################################################################
############################################# truncated normal mix ##################################################
#######################################################################################################

fitNL_trunc <- function(y, model='E'){
	if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
	percent0 <- mean(y==0, na.rm=TRUE)
	res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
	names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
#browser()
	if(percent0<=zeroPercentThr){ # no severe 0-inflation, fit 2-comp
		Dat <- (log(y+eps)-log(d))/log(base) # transform data; change of base formula: log_b(x)=log_e(x)/log_e(b)
		mc <- try(Mclust(Dat, G = 2, modelNames = model), silent=TRUE)
		res[1:7] <- extractMclustPar(mc, modelName=model, logLikToLN=logLikToLN, dat=Dat)
	} else { # severe 0-inflation. 1-comp model. need control mu and sd especially for RPKM
		nonzero <- y[y!=0]
		#if(mean(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001, force it at least 1 to avoid FP
		#### median is more robust: not susceptible to single outlier expression
		#### bug here: 10/23/2012: when 0-inflated with small RPKM, i.e. 0.002, the estimated mu2=eps(1) due to this artifact
		#if(median(nonzero)<1) nonzero <- nonzero+max(c(1, eps)) ## when count<1, i.e. in RPKM, usually 0.00001 (mu2=-5),  force it at least 1 to avoid FP
		Dat <- (log(nonzero)-log(d[y!=0]))/log(base) ###### no need to add eps since no 0 now!
		#mc_1comp <- try(Mclust(Dat, G = 1), silent=TRUE)
		#Modified to MLE such that c(00000, 11111) wouldn't fail
		mu2 <- mean(Dat)
		if(mu2<0 & median(nonzero)<1) mu2 <- 0 # mu2 needs to be larger than mu1 when 0-inflated
		sigma2 <- max(sd(Dat), 1e-2) ## in case sd=0.0004
		pi1 <- percent0 # 
		res[1:5] <- c(0, mu2, 0, sigma2, pi1) # mu1=0, sigma1=0, logLik=BIC=NA
		res[6] <- logLik0inflatedLN(y, mu=mu2, sigma=sigma2, pi1=pi1, logLikToLN=logLikToLN)
		res[7] <- getBIC(logLik=res[6], nPar=3, nObs=length(y)) 
	}
	res
}


# a function to extract parameters estimated by mclust with G=2, compatible to try-error of Mclust() function
## modified on 07/22/12: add option logLikToLN and dat. by default, it's disabled so that it is compatible with previous version
## add option logLikToLN on 7/22/2012 so that logLik and BIC can be defined on the exponent LN scale. 
## This is needed for pseudo-LN mixture fitting		
# dat: the data exactly used to fit the normal mixture. it is defined on the log scale and hence can be negative		
extractMclustPar <- function(mc, modelName='E', logLikToLN=FALSE, dat=NA){
res <- rep(NA, 7)
nPar <- ifelse(modelName=='V', 5, 4) # number of parameters
if(class(mc)!="try-error"){
	# extract mu1, mu2
	res[1:2] <- mc$parameters$mean
	# extract sigma1, sigma2
	temp <- sqrt(mc$parameters$variance$sigmasq)
	if(length(temp)==1){ # E model, 1 sigma
		res[3:4] <- rep(temp, 2)
	} else {
		res[3:4] <- temp
	}
	# extract p1
	res[5] <- mc$parameters$pro[1]
	if(logLikToLN){ # modify logLik and BIC to the exponent LN scale.
	# extract logLik
	res[6] <- logLikLN(y=exp(dat), theta=res[1:5]) # data is now transformed to the exponent scale
	# extract BIC
	res[7] <- getBIC(logLik=res[6], nPar=nPar, nObs=length(dat)) ## confirmed by BIC function in R
	} else {
	# extract logLik
	res[6] <- mc$loglik
	# extract BIC
	res[7] <- ifelse(modelName=="V", -bic(modelName="V", loglik=mc$loglik, n=mc$n, d=1, G=2), 
									 -bic(modelName="E", loglik=mc$loglik, n=mc$n, d=1, G=2))
	}								 
}
res
}



## logLik for LN data. used to compute BIC of LN since we fit normalMixture whose logLik is defined on the log-transformed data
# theta: a vector. if length is 2, theta=c(meanlog , sdlog),logLik of 1-comp data is obtained. 
#		 if length=5, logLik of 2-comp data is obtained. theta=c(meanlog1, meanlog2, sdlog1, sdlog2, pi1)
## bug fixed at 07/27/12:  dlnorm(y, meanlog=theta[2], sdlog=theta[5], log=FALSE) --->  dlnorm(y, meanlog=theta[2], sdlog=theta[4], log=FALSE)
logLikLN <- function(y, theta){
	if(length(theta)==2){ # 1-comp logLik
		res <- sum(dlnorm(y, meanlog=theta[1], sdlog=theta[2], log=TRUE))
	} else { # 2-comp logLik
		res <- sum(log(dlnorm(y, meanlog=theta[1], sdlog=theta[3], log=FALSE)*theta[5]+
				   dlnorm(y, meanlog=theta[2], sdlog=theta[4], log=FALSE)*(1-theta[5])))
	}
	res
}




#######################################################################################################
############################################# NL mix ##################################################
#######################################################################################################
# for completeness, we also incorporate normal mixture for microarray data. We use NL to denote normal.
#' Fit Normal Mixture Model
#' The function fits a two-component normal mixture model. 
#'
#' The parameter estimates from log normal mixture is obtained by taking logarithm and fit normal mixture. We use
#' mclust package to obtain parameter estimates of normal mixture model. In particular, \eqn{log_{base}(\frac{y+eps}{d})} is used to
#' fit to normal mixture model. 
#' 
#' With this function, three models can be fitted: (1) log normal mixture with equal variance (E model); 
#' (2) Generalized Poisson mixture with unequal variance (V model); (3) 0-inflated log normal model.
#' The 0-inflated log normal has the following density function:
#' 
#' \eqn{P(Y=y)=\pi D(y) + (1-\pi)LN(\mu, \sigma)} where D is the point mass at 0 while \eqn{LN(\mu, \sigma)} is the density
#' of log normal distribution with mean \eqn{\mu} and variance \eqn{\sigma^2}. 
#' 
#' The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
#' rule overrides the model argument (E or V) when observed percentae of zero count exceeds the threshold.
#' 
#' @param y A vector representing the RNAseq raw count.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data. For the LN model, the original data would be devided by this vector. 
#' @param model Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. 
#' @return A vector consisting parameter estimates of mu1, mu2, sigma1, sigma2, pi1, logLik and BIC. 
#' @references Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes Using RNAseq Data. Bioinformatics, submitted.
#' @export
fitNL <- function(y, d=NULL, model='E') {
	if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
	model <- try(match.arg(model, c('E', 'V'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model E or V can be specified!\n')
	res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
	names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
	Dat <- y/d # normalization
	mc <- try(Mclust(Dat, G = 2, modelNames = model), silent=TRUE)
	res[1:7] <- extractMclustPar(mc, modelName=model, logLikToLN=FALSE, dat=Dat)
	res
}

# given mu1, mu2, phi1, phi2, p1, logLik (optional), BIC (optional), transform phi to sigma
# transformNBfit() is its special case
transformPhiToSigmaInMixFit <- function(est, dist='NB'){
	sigma1 <- sqrt(muPhiToVar(mu=est[1], phi=est[3], dist=dist))
	sigma2 <- sqrt(muPhiToVar(mu=est[2], phi=est[4], dist=dist))
	res <- est
	res[3:4] <- c(sigma1, sigma2)
	res
}
# phi2Var: when mu1, mu2, phi1, phi2, pi1 is specified, we need to set this so that var can be computed
# modified on 08/03/12: deal with 0-inflation: input=c(0, mu2, 0, sigma2, pi1, ...), delta=mu2/sigma2 
# modified on 08/03/12: !any(is.na(est))----->!any(is.na(est[1:5])). necessary since 0-inflation always get logLik and BIC as NA
# modified on 08/04/12: for 0-inflation case, general BI2 formula also works. instead sqrt((1-pi))*mu2/sigma2 is wrong: no detection of such genes!
parToBI <- function(est, phi2Var=FALSE, dist='NB'){
	res <- c(NA, NA)
	names(res) <- c('delta', 'BI')
	if(!any(is.na(est[1:5]))){
		if(phi2Var==TRUE){
			est[3] <- sqrt(muPhiToVar(est[1], est[3], dist=dist))
			est[4] <- sqrt(muPhiToVar(est[2], est[4], dist=dist))
		}
		#est[5] <- restrictPi1(est[5]) # deal with numeric issue from optim fit (GP). i.e. pi= -2.775558e-17
		# BI2 automatically deal with 0-inflation 
		# transformAlphaBetaToMuSigma: names(res) <- c('mean1',  'mean2', 'sigma1','sigma2', 'pi1')
		res[1] <- abs(diff(est[1:2]))/sqrt((1-est[5])*est[3]^2+est[5]*est[4]^2)
		res[2] <- sqrt(est[5]*(1-est[5]))*res[1]
	}
	if(!is.na(est[5]) & (est[5]==0 | est[5]==1)){
		res[1] <- res[2] <- 0 # 1-component data, BI=0
	}
	res
}

noNA <- function (dat, returnIndex = FALSE) 
{
    sel <- complete.cases(dat)
    if (returnIndex) 
        return(sel)
    if (is.null(dim(dat))) {
        res <- dat[sel]
    }
    else {
        res <- dat[sel, ]
    }
    res
}

#######################################################################################################
############################################# SIBER ###################################################
#######################################################################################################
# the final product presented to the user
# parameters (zeroPercentThr=0.1, base=10, eps=10) are only relevant to LN model.
#' Fit Mixture Model on The RNAseq Data and Calculates Bimodality Index
#' The function fits a two-component mixture model and calculate BI from the parameter estimates. 
#'
#' SIBER proceeds in two steps. The first step fits a two-component mixture model. 
#' The second step calculates the Bimodality Index corresponding to the assumed mixture distribution.
#' Four types of mixture models are implemented: log normal (LN), Negative Binomial (NB),  Generalized Poisson (GP), Beta (Beta) and normal mixture (NL). The normal mixture model was developed to identify bimodal genes from microarray data in Wang et al. It is incorporated here
#' in case the user has already transformed the RNAseq data. The Beta mixture model can be applied to methylation data where the observed values are between 0 and 1 representing metylation rate. 
#' Behind the scene, SIBER calls the fitNB, fitGP, fitLN and fitNL function with model=E depending on which
#' distribution model is specified. When the observed percentage of count exceeds the user specified threshold
#' zeroPercentThr, the 0-inflated model overrides the E model and will be fitted. 
#' Type vignette('SIBER') in the R console to pull out the user manual in pdf format. 
#' @param y A vector representing the RNAseq raw count or the transformed values if model=NL.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data.
#' @param model Character string specifying the mixture model type. It can be any of LN, NB, GP, Beta and NL.
#' @param zeroPercentThr A scalar specifying the minimum percent of zero to detect using log normal mixture. This
#' parameter is used to deal with zero-inflation in RNAseq count data. When the percent of zero exceeds this threshold, 
#' 1-comp mixture LN model is used to estimate mu and sigma from nonzero count. This parameter is relevant only if model='LN'.
#' @param base The logarithm base defining the parameter estimates in the logarithm scale from LN model . It is relevant only if model='LN'.
#' @param eps A scalar to be added to the count data when model='LN'. This parameter is relevant only when model='LN'. 
#' @return A vector consisting estimates of mu1, mu2, sigma1, sigma2, p1, delta and BI. 
#' @references Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes Using RNAseq Data. Bioinformatics, submitted.
#' @export
#' @examples
#' set.seed(100)
#' y=c(rbeta(100,1,4),rbeta(200,4,1))
#' fitBeta(y=y)
#' SIBER(y, model='Beta')
SIBER <- function(y, d=NULL, model=c('LN', 'NB', 'GP', 'Beta', 'NL', 'BetaReg'), zeroPercentThr=0.2, base=exp(1), eps=10){
	model <- try(match.arg(model, c('LN', 'NB', 'GP', 'Beta', 'NL', 'BetaReg'), several.ok=FALSE), silent=TRUE)
	# stop if model not recognizable.
	if(class(model)=='try-error') stop('Only model LN, NB, GP, Beta or NL can be specified!\n')
	res <- rep(NA, 7) 
	names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'delta', 'BI') 
	# if y has NA, all result will be NA; thus need to remove NA beforehand
	y <- noNA(y)
	if(model=='LN'){
		fit <- fitLN(y, d=d, model='E', zeroPercentThr=zeroPercentThr, base=base, eps=eps)[1:5]
	} else if(model=='NB') {
		fitPhiScale <- fitNB(y, d=d, model='E')[1:5]
		fit <- transformPhiToSigmaInMixFit(fitPhiScale, dist='NB')
	} else if(model=='GP') {
		fitPhiScale <- fitGP(y, d=d, model='E')[1:5]
		fit <- transformPhiToSigmaInMixFit(fitPhiScale, dist='GP')
	} else if(model=='Beta') {
        fitAlphaBetaScale <- fitBeta(y, model='V')
        fit <- transformAlphaBetaToMuSigma(fitAlphaBetaScale)
    } else if(model=='BetaReg') {
        fitAlphaBetaScale <- fitBetaReg(y)
        #browser()
        fit <- transformAlphaBetaToMuSigma(fitAlphaBetaScale)
    } else {
		#browser()
		fit <- fitNL(y, d=d, model='E')[1:5]
	}
	BIinfo <- parToBI(fit, phi2Var=FALSE)
	res[1:5] <- fit
	res[6:7] <- BIinfo
	#browser()
	res
}

#' fit SIBER on a matrix, one row at a time
#'
#' this function enables parallel through plyr package
#' 
#' @param mat matrix of expression, methylation and so on
#' @param model model as in SIBER
#' @param prune logical, whether to prune sigma1, sigma2 for with a specified percentile; this is to obtain stable BI estimate due to extremely small
#' sigma estimate, especially in the V model
#' @param parallel logical
#' @param core number of cores to register for parallel computing
#' @return a matrix
#' @references Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes Using RNAseq Data. Bioinformatics, submitted.
#' @export
fitSIBERonMat <- function(mat, model='NL', prune=TRUE, q=0.03, parallel=TRUE, core=10){
    require(SIBER)
    if(parallel){
		# required variables (including functions) and packages are exported here
		cl <- createCluster(core=core, logfile = "/dev/null", export='SIBER2', lib = c('SIBER'))
		on.exit(stopCluster(cl))
	}
	#browser()
	tmpRes <- adply(mat, 1, SIBER2, model=model, .parallel=parallel)
    res <- moveColumnToRowName(tmpRes)
	#browser()
	if(prune){
		# truncate sigma
		# res[, 'sigma1'] <- truncByQuantile(res[, 'sigma1'], q1=q, q2=1)
		# res[, 'sigma2'] <- truncByQuantile(res[, 'sigma2'], q1=q, q2=1)
		tpmin <- quantile(c(res[, 'sigma1'], res[, 'sigma2']), prob=q, na.rm=T)
		res[, 'sigma1'] <- truncByLimit(res[, 'sigma1'], Lower=tpmin, Upper=Inf)
		res[, 'sigma2'] <- truncByLimit(res[, 'sigma2'], Lower=tpmin, Upper=Inf)
		#browser()
		# update delta and BI
		res <- updateBImat(res)
	}
    res
}

#' update the result from BI data frame. 
#' @param dat a data frame with results fit from SIBER
#' @export
updateBImat <- function(dat){
	# 1 comp data: pi1=1 or pi1=0: no update
	is1comp <- !is.na(dat[, 'pi1'])  & (dat[, 'pi1']==0 | dat[, 'pi1']==1) 
	is2comp <- !is1comp
	# only update if it is 2-comp data
	t1 <- abs(dat[, 'mu2']-dat[, 'mu1'])/sqrt((1-dat[, 'pi1'])*dat[, 'sigma1']^2+dat[, 'pi1']*dat[, 'sigma2']^2)
	t2 <- sqrt(dat[, 'pi1']*(1-dat[, 'pi1']))*dat[, 'delta']
	dat[is2comp, 'delta'] <- t1[is2comp]
	dat[is2comp, 'BI'] <- t2[is2comp]
	dat
}
#'  Simulated Data From 2-component Mixture Models
#'
#' Data from 2-component mixture models (NB, GP and LN) is simulated
#' with the true parameters given for testing and illustration purpose.
#'
#' The data frame contains the following data objects:
#' \item{parList}{A list of true parameters. There are three named elements 
#'	(NB, GP and LN) corresponding to
#'	the parameters used to simulate gene expression data from NB, GP and LN mixture models. Each
#'	element is a 6 by 5 matrix giving the true parameters generating the simulated data. }
#'    \item{dataList}{A list of matrices for simulated gene expression data. 
#'	There are three named elements (NB, GP and LN) corresponding to
#'	the simulate gene expression data from NB, GP and LN mixture models. Each
#'	element is a 6 by 200 matrix. That is, 6 genes (rows) are simulated with 200 samples (columns). 
#'	The first 3 genes in each matrix are from 2-component mixture model while the last 3 genes
#'	are from 0-inflated models. }
#'
#' @references Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes Using RNAseq Data. Bioinformatics, submitted.
#' @name simDat
NULL