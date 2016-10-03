require('betareg')
#inv_logit <- SIBER:::inv_logit # betamix returns logit(mean); this transforms back
getAlpha <-function(mu,phi) inv_logit(mu)*exp(phi) # from logit(mu), log(precision=a+b) to alpha
getBeta <-function(mu,phi) (1-inv_logit(mu))*exp(phi) # from logit(mu), log(precision=a+b) to beta
getMean <-function(mu,phi) inv_logit(mu)
getPrecision <-function(mu,phi) exp(phi)
getVar <-function(mu,phi) getMean(mu,phi)^2*getBeta(mu, phi)/getAlpha(mu, phi)/(getPrecision(mu, phi)+1)

#getVar(mu1, phi1)
fitBetaReg <- function(y){
	res <- rep(NA, 7)	# res=c(alpha1, beta1, alpha2, beta2 and pi1, logLik, BIC)
    names(res) <- c('alpha1', 'beta1', 'alpha2', 'beta2', 'pi1', 'logLik', 'BIC')
    require('betareg')
	y <- noNA(y)
	#browser()
	rs_mix <- try(betamix(y ~ 1, link="logit", link.phi="log", k = 2), silent=TRUE)
	if(class(rs_mix)!='try-error'){
		out=try(summary(rs_mix), silent=TRUE)
		if(class(out)!='try-error'){
			if(Getter(out, 'k')==2){ # sometimes only 1 component is fitted, thus cannot give all the result
				#capture.output( out=summary(rs_mix), file='NUL')
				logitMu1 = out@components[[1]][[1]]$mean[1]
				logPhi1 = out@components[[1]][[1]]$precision[1]
				logitMu2 = out@components[[1]][[2]]$mean[1]
				logPhi2 = out@components[[1]][[2]]$precision[1]
				# temp result: alpha1 <--> alpha2 may be switched if mu2<mu1
				alpha1_ = getAlpha(logitMu1, logPhi1)
				beta1_ = getBeta(logitMu1, logPhi1)
				pi1_ = 1-inv_logit(out@coef[5])
				alpha2_ = getAlpha(logitMu2, logPhi2)
				beta2_ = getBeta(logitMu2, logPhi2)
				res['logLik'] <- logLik <- Getter(rs_mix$flexmix, 'logLik')
				df = Getter(rs_mix$flexmix, 'df')
				res['BIC'] = SIBER:::getBIC(logLik=logLik, nPar=df, nObs=length(y)) 
				# ensure mu1 < mu2
				if(logitMu1 > logitMu2){
					res['alpha1'] <- alpha2_
					res['beta1'] <- beta2_
					res['alpha2'] <- alpha1_
					res['beta2'] <- beta1_
					res['pi1'] <- 1-pi1_
				} else {
					res['alpha1'] <- alpha1_
					res['beta1'] <- beta1_
					res['alpha2'] <- alpha2_
					res['beta2'] <- beta2_
					res['pi1'] <- pi1_
				}
			}
			if(Getter(out, 'k')==1){ # sometimes only 1 component is fitted, thus set pi1=1, alpha2=beta2=NA, BI=0
				#browser()
				#capture.output( out=summary(rs_mix), file='NUL')
				logitMu1 = out@components[[1]][[1]]$mean[1]
				logPhi1 = out@components[[1]][[1]]$precision[1]
				# temp result: alpha1 <--> alpha2 may be switched if mu2<mu1
				alpha1_ = getAlpha(logitMu1, logPhi1)
				beta1_ = getBeta(logitMu1, logPhi1)
				res['logLik'] <- logLik <- Getter(rs_mix$flexmix, 'logLik')
				df = Getter(rs_mix$flexmix, 'df')
				res['BIC'] = SIBER:::getBIC(logLik=logLik, nPar=df, nObs=length(y)) 
				# ensure mu1 < mu2
				res['alpha1'] <- alpha1_
				res['beta1'] <- beta1_
				res['alpha2'] <- NA
				res['beta2'] <- NA
				res['pi1'] <- 1
			}
		}
		
	}
	res	
	#browser()
}
#fit <- fitBetaReg(y=ReadingSkills$accuracy)
#SIBER(y=ReadingSkills$accuracy, model='BetaReg')