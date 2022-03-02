#' Calculate predictions using scenario accounting for estimation uncertainty
#'
#' This function generates predictions from model and inputted scenario
#' @param modObj model object (glm or lm)
#' @param scen scenario object generated from scenBuild
#' @param link link function, options include: logit, probit (we'll build more in as we go)
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
  beta,
  varcov,
  scen,
  link='logit',
  seed=6886, sims=100
){

  # set seed and obtain draws from model object
  set.seed(seed)
  draws = mvtnorm::rmvnorm(sims, beta, varcov)

  # get predictions and apply link function
  preds = lapply(scen$scen, function(s){

    # calc xbeta
    xBeta = draws %*% t(s)

    # apply link
    if(link=='logit'){
      preds = 1/(1+exp(-xBeta)) }
    if(link=='probit'){
      preds = pnorm(xBeta) }
    if(link=='normal'){
      preds = xBeta }

    #
    return(preds) })

  # if observed value approach used get average of predictions
  if(scen$scenType=='observed'){
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
  names(preds)[ncol(preds)] = 'pred'

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
