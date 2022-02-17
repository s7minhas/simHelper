# miscellaneous functions for use in package

#' char
#'
#' This function allows users to input country and time
#' information to obtain an undirected dyadic dataset
#' @param x vector


char <- function(x){
	as.character(x)
}


#' num
#'
#' This function allows users to input country and time
#' information to obtain an undirected dyadic dataset
#' @param x vector

num <- function(x){
	as.numeric(char(x))
}


#' replScenVal
#'
#' This function helps users to change the values of
#' specific columns in a matrix
#' @param scen a matrix object with column names
#' @param var a character argument that corresponds to a column of the matrix
#' @param newval a numeric argument giving the new value for the scen matrix
#' @export

replScenVal <- function(scen,var,newval){
  scen[,var]=newval
  return(scen) }
