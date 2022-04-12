#' Calculate predictions using scenario accounting for estimation uncertainty
#'
#' This function generates predictions from model and inputted scenario
#' @param modObj model object (glm or lm), it's recommended to just supply values to beta and varcov
#' @param beta vector of regression estimates, usually obtained from coef of the model object
#' @param varcov variance covariance matrix
#' @param scen scenario object generated from scenBuild
#' @param link link function, options include: binary-logit, binary-probit, count, normal
#'
#'
#' @return a data.frame object that includes columns for each of the treatment values
#' and then a final column with the predictions for that treatment. the number of
#' predictions generated per treatment correspond to the number of simulations (sims)
#' run.
#' @keywords predictions
#'
#' @export getPreds

getPreds <- function(
  modObj=NULL,
  beta=NULL,
  varcov=NULL,
  scen,
  link='binary-logit',
  seed=6886, sims=100
){

  # set seed and obtain draws from model object
  set.seed(seed)
  if(is.null(beta) & is.null(varcov)){
    draws = mvtnorm::rmvnorm(sims, coef(modObj), vcov(modObj)) }
  if(!is.null(beta) & !is.null(varcov)){
    draws = mvtnorm::rmvnorm(sims, beta, varcov) }

  # get predictions and apply link function
  preds = lapply(scen$scen, function(s){

		# process into pred for non-ordinal families
		if(!grepl('ordinal', link)){

	    # rearrange cols to be in the same order
	    if( 'intercept' %in% colnames(s) ){
	      s[,2:ncol(s)] = s[,colnames(draws)[-1]]
	    } else { s = s[,colnames(draws)] }

	    # calc xbeta
	    xBeta = draws %*% t(s)

	    # apply link
	    if(link=='binary-logit'){
	      preds = 1/(1+exp(-xBeta)) }
	    if(link=='binary-probit'){
	      preds = pnorm(xBeta) }
	    if(link=='count'){
	      preds = exp(xBeta) }
	    if(link=='normal'){
	      preds = xBeta }   }

		# get preds for ordinal families
		if(grepl('ordinal', link)){

			# rearrange cols to be in the same order
			s = s[,intersect(colnames(draws), colnames(s)),drop=FALSE]

			# organize results from draws
			betas = draws[,1:ncol(s)]
			taus = draws[,(ncol(s)+1):ncol(draws)]

			# calculate prob of lowest cat
			preds = lapply(1:ncol(taus), function(ii){

				# get prob of lowest cat
				if(ii == 1) {
					pred = taus[,ii] - betas%*%t(s)
					if(link=='ordinal-logit'){  pred = plogis(pred)}
					if(link=='ordinal-probit'){  pred = pnorm(pred)}
					return(pred) }

				# get prob of middle cats
				if(ii != 1) {
					pred0 = taus[,ii-1]-betas%*%t(s)
					pred1 = taus[,ii]-betas%*%t(s)
					if(link == 'ordinal-logit'){ pred0=plogis(pred0) ; pred1=plogis(pred1) }
					if(link == 'ordinal-probit'){ pred0=pnorm(pred0) ; pred1=pnorm(pred1) }
					pred = pred1 - pred0
					return(pred) } })

			# if observed vals then avg here
			if(scen$scenType=='observed'){
				preds = lapply(preds, function(predCat){
					matrix(apply(predCat, 1, mean), ncol=1) }) }

			# organize and get preds for final cat
			preds = do.call('cbind', preds)
			preds = cbind(preds, apply(preds, 1, function(x){ 1-sum(x) }) )
			colnames(preds) = paste0('y', 1:ncol(preds))
			return(preds) }

    #
    return(preds) })

  # if observed value approach used get average of predictions
  if(scen$scenType=='observed' & !grepl('ordinal', link)){
    preds = lapply(preds, function(pred){
      avgPred = matrix(apply(pred, 1, mean), ncol=1)
      return(avgPred) }) }

  # add treatment variable info
  preds = lapply(1:nrow(scen$treatCombo), function(ii){
    predT = preds[[ii]]
    tVals = scen$treatCombo[ii,,drop=FALSE]
    rownames(tVals) = NULL
    predT = cbind(tVals, pred = predT)
    return(predT) })

  # reorg
  preds = do.call('rbind', preds)
	if(!grepl('ordinal', link)){ names(preds)[ncol(preds)] = 'pred'	}

  #
  return(preds)
}

# #####
# library(simHelper)
# data(aNigeria)
# #####
#
# #####
# mod = glm(
#   conflict ~
#     lagDV + lagRecip + govActor + postBoko +
#     elecYear + groupSpread + ngbrConfCount,
#   data=aNigeria,
#   family=binomial(link='logit')
# )
#
# bkHypScen = scenBuild(
#   mData=aNigeria,
#   ivs = c('lagDV', 'lagRecip', 'govActor', 'postBoko', 'elecYear', 'groupSpread', 'ngbrConfCount'),
#   scenType='hypothetical',
#   treatVar = 'postBoko',
#   treatCategorical = TRUE)
#
# bkHypPreds = getPreds(
#   modObj = mod,
#   scen=bkHypScen,
#   sims=1000)
#
# bkObsScen = scenBuild(
#   mData=aNigeria,
#   ivs = c('lagDV', 'lagRecip', 'govActor', 'postBoko', 'elecYear', 'groupSpread', 'ngbrConfCount'),
#   scenType='observed',
#   treatVar = 'postBoko',
#   treatCategorical = TRUE)
#
# bkObsPreds = getPreds(
#   modObj = mod,
#   scen=bkObsScen,
#   sims=1000)
#
# head(bkObsPreds)
# #####
#
# #####
# modInt = glm(
#   conflict ~
#     lagDV + lagRecip + govActor + postBoko +
#     elecYear + groupSpread + ngbrConfCount + postBoko*ngbrConfCount,
#   data=aNigeria,
#   family=binomial(link='logit'))
#
# bkNgbrHypScen = scenBuild(
#   mData=aNigeria,
#   ivs = c('lagDV', 'lagRecip', 'govActor', 'postBoko', 'elecYear', 'groupSpread', 'ngbrConfCount'),
#   scenType='hypothetical',
#   treatVar = c('postBoko', 'ngbrConfCount'),
#   treatCategorical = c(TRUE, FALSE))
#
# bkNgbrHypPreds = getPreds(
#   modObj = modInt,
#   scen=bkNgbrHypScen,
#   sims=1000)
#
# bkNgbrObsScen = scenBuild(
#   mData=aNigeria,
#   ivs = c('lagDV', 'lagRecip', 'govActor', 'postBoko', 'elecYear', 'groupSpread', 'ngbrConfCount'),
#   scenType='observed',
#   treatVar = c('postBoko', 'ngbrConfCount'),
#   treatCategorical = c(TRUE, FALSE))
#
# bkNgbrObsPreds = getPreds(
#   modObj = modInt,
#   scen=bkNgbrObsScen,
#   sims=1000)
#
#
# library(dplyr)
# library(ggplot2)
# preds = bkNgbrObsPreds
# preds = preds %>%
#   group_by(postBoko, ngbrConfCount) %>%
#   summarize(
#     mu = mean(pred),
#     hi = quantile(pred, .975),
#     lo = quantile(pred, .025)
#   )
# preds$postBoko = factor(preds$postBoko)
# ggplot(preds, aes(x=ngbrConfCount, y=mu, color=postBoko, fill=postBoko, group=postBoko)) +
#   geom_line() +
#   geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.3)
# #####

# ############
# library(simHelper)
# data(anes)
# library(MASS)
#
# set.seed(6886)
#
# # make sure R reads dv as ordered factor
# anes$obamaLR2 = as.ordered(anes$obamaLR2)
#
# # set up model spec
# dv = 'obamaLR2'
# ivs = c('pid', 'age', 'education', 'income', 'white')
# form = paste0(dv, '~', paste(ivs, collapse='+'))
#
# # run model
# ordMod <- polr( form, data=anes, method="logistic", Hess=T)
#
# scen = scenBuild(
# 	anes,
# 	ivs = ivs,
# 	scenType = 'observed',
# 	# ivValues = list(age=30, education=4, white=1),
# 	treatVar = 'pid',
# 	treatCategorical=FALSE,
# 	intercept=FALSE
# )
#
# ordPreds = getPreds(
# 	beta=c(coef(ordMod), ordMod$zeta),
# 	varcov=solve(ordMod$Hessian),
# 	scen=scen,
# 	link='ordinal-logit',
# 	seed=6886, sims=100
# )
#
# head(ordPreds)
#
# scen2 = scenBuild(
# 	anes,
# 	ivs = ivs,
# 	scenType = 'hypothetical',
# 	# ivValues = list(age=30, education=4, white=1),
# 	treatVar = 'pid',
# 	treatCategorical=FALSE,
# 	intercept=FALSE
# )
#
# ordPreds2 = getPreds(
# 	beta=c(coef(ordMod), ordMod$zeta),
# 	varcov=solve(ordMod$Hessian),
# 	scen=scen2,
# 	link='ordinal-logit',
# 	seed=6886, sims=100
# )
#
# dim(ordPreds)
# dim(ordPreds2)
