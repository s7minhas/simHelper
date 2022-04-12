#' Helper to build scenarios
#'
#' This function helps create scenarios for simulation analysis
#' @param mData data.frame object that contains at least the variables to construct the scenario
#' @param ivs character vector with all independent variables to be used
#' @param scenType type of scenario to construct: observed value approach, hypothetical, or counterfactual. if using counterfactual then you must supply values for each iv to ivValues in a named list
#' @param ivStats if hypothetical chosen then you have the option of specifying the summary stats to use
#' @param ivValues a named list where the names should exactly match variables in the mData input and values should be those that you wish to generate a prediction for
#' @param treatVar character vector indicating the treatment or treatments to be used
#' @param treatCategorical logical vector indicating whether each of the treatments specified is categorical or not
#' @param treatVals can supply a set of treatment values to test, if not provided then they will be calculated using the range of the data, default is NULL
#' @param intercept a logical indicating whether model was estimated with an intercept, default is TRUE
#' @param treatVals treatment values to test that should be wrapped in a list
#'
#'
#' @return a list, the first element being a list of scenario matrices, one for each
#' treatment value. second being a data.frame of the treatment values corresponding
#' to each scenario. third is a character indicating the type of scenario (hypothetical
#' or observed. last being a logical indicating whether or not the treatment variables
#' are categorical.
#' @keywords scenario
#'
#' @export scenBuild

# library(simHelper)
# data(saleyhan_etal_2011)
# mData = saleyhan_etal_2011
# ivs = c(
#   'strong', 'cl', 'tk',
#   'gs', 'riv', 'pol6', 'cinc', 'log_rgdpch' )
# ivStats=rep('median',length(ivs))
# scenType = 'counterfactual'
#
# # treatVar = c('cinc')
# # treatCategorical=c(FALSE)
#
# treatVar = c('tk')
# treatCategorical=c(TRUE)
#
# ivValues = list(
# 	strong =1, cl =1 , tk =1,
# 	gs = 0, riv = 1, pol6 = 0, cinc = mean(mData$cinc)
# )
#
# treatVals = NULL
# intercept=TRUE
#
# library(simHelper)
# data(anes)
# 	mData=anes
# 	ivs = c('pid', 'age', 'education', 'income', 'white')
# 	scenType = 'counterfactual'
# 	ivValues = list(age=30, education=4, income=6, white=1)
# 	treatVar = 'pid'
# 	treatCategorical=FALSE
# 	intercept=FALSE
# 	treatVals=NULL
# 	rep('median',length(ivs))

scenBuild <- function(
  mData,
  ivs,
  scenType='hypothetical',
  ivStats=rep('median',length(ivs)),
  ivValues=NULL,
  treatVar,
  treatCategorical=TRUE,
  treatVals=NULL,
  intercept=TRUE
){

  # convert to df and remove NAs from data
  mData = na.omit(data.frame(mData[,ivs], stringsAsFactors=F))

  # make sure all values in mData are numeric
  mData = apply(mData, 2, function(x){ as.numeric(as.character(x)) })

  # get treatment values if not provided
  treatValList = lapply(1:length(treatVar), function(ii){

    # get characs of treatval
    tV = treatVar[ii]
    tvCat = treatCategorical[ii]
    tVals = treatVals[[ii]]

    # get treatment values if they are not provided
    # and treatment is categorical
    if(is.null(tVals) & tvCat){
      tVals = sort(unique(mData[,tV])) }

    # get treatment values if they are not provided
    # and treatment is continuous
    if(is.null(tVals) & !tvCat){
      minVal = min( mData[,tV] )
      maxVal = max( mData[,tV] )
      tVals = seq(minVal, maxVal, length.out=20) }

    #
    return(tVals) })
  names(treatValList) = treatVar
  treatCombo = do.call('expand.grid', treatValList)

  # create base scenario if hypothetical chosen
  if(scenType=='hypothetical'){

    # calculate central tendency of vars in model
    bScen <- lapply( 1:length(ivs), function(ii){
      baseStat = eval(call(ivStats[ii], mData[,ivs[ii]]))
      return(baseStat) } )
    bScen = do.call("cbind", bScen)
    colnames(bScen) = ivs

    # add in values from treatCombo
    scen = lapply(1:nrow(treatCombo), function(ii){
      tVals = data.matrix(treatCombo[ii,,drop=FALSE])
      out = bScen
      out[,colnames(tVals)] = tVals
      out = matrix(out,
        nrow=nrow(bScen), ncol=ncol(bScen),
        dimnames=list(rownames(bScen), colnames(bScen)) )
      return(out) })
  } # end hypothetical scenario construction

  # create scenario based on inputted values
  if(scenType=='counterfactual'){

		#
		if(is.null(ivValues)){
			warning('
			when creating counterfactual scenario,
			supply values for all ivs in named list
			') }

    # find vars to add
    toAdd = setdiff(ivs, treatVar)

		# for variables where no iv value
		# is assigned, use ivStat
		missVars = setdiff(toAdd, names(ivValues))
		if( length(missVars) > 0 ){
		for(ii in 1:length(missVars)){
			ivValues[[missVars[ii]]] = eval(
				call(
					ivStats[ match(missVars[ii], ivs) ],
					mData[,missVars[ii]] ) ) } }

    # create scen base
    scen = treatCombo

    # add values from ivValues
    for(add in toAdd){
      scen = cbind(scen, ivValues[[add]])
      names(scen)[ncol(scen)] = add
    }

    # convert to list
    scen = lapply(1:nrow(scen), function(ii){ scen[ii,,drop=FALSE] })
  }

  # create scenario based on observed value approach
  if(scenType=='observed'){

    # get unique of observed data
    bScen = unique(mData)

    # repeat for every value of treatment
    scen = lapply(1:nrow(treatCombo), function(ii){
      tVals = data.matrix(treatCombo[ii,,drop=FALSE])
      tVals = matrix(
        tVals, nrow=nrow(bScen), ncol=ncol(tVals), byrow=TRUE,
        dimnames=list(NULL, colnames(tVals)) )
      bScen[,colnames(treatCombo)] = tVals
      return(bScen) })
  } # end observed value scenario construction

	# add intercept if relevant
	if(intercept){ scen = lapply(scen, function(s){ cbind(intercept=1, s) }) }

  # if more than one treatVar create interaction columns
  if(length(treatVar)>1){
    intF = formula(paste0('~',paste(treatVar, collapse='*')))
		scen = lapply(scen, function(s){
			toAdd = model.matrix(intF, data.frame(s))
	    toAdd = toAdd[,-c( 1:( length(treatVar)+1 ) ), drop=FALSE]
	    return( cbind(s, toAdd) ) } ) }

  #
	out = list(
		scen=scen, treatCombo=treatCombo,
		scenType=scenType,
		treatCategorical=treatCategorical
	)
  return(out)
}
