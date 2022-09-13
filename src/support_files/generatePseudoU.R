#' Pseudo-observation generator
#'
#' generatePseudoU() takes as input a matrix or dataframe of observations
#' and returns a matrix of ranked pseudo-observations
#'  where d is the number of columns in the input.
#'
#'  The observations are normalized such that the value 1 is never
#' returned, as this could generate some numerical difficulties
#' later on.
#'
#' @param observedSample Matrix or dataframe of observations
#' @keywords
#' @export
#' @examples
#' generatePseudoU()
#'

generatePseudoU <- function (observedSample){
  udata = matrix(0,dim(observedSample)[1],dim(observedSample)[2])

  n<-length(observedSample[,1])
  for(i in 1:dim(observedSample)[2] ){
    udata[,i]<-rank(observedSample[,i])/(n+1)
  }
  return(udata)
}
