#' Univaraite margin fitting
#'
#'
#' @param observedMargin ECDF values
#' @param listOfModels Univariate margin families to evaluate
#' @param quantilesToEval Quantiles to calculate values for
#' @keywords
#' @export
#' @examples
#' fitMarginFunction()
#'
fitMarginFunction <- function (observedMargin,listOfModels,quantilesToEval){
  fitLc = list()
  aicScore = matrix()
  for (i in 1:length(listOfModels)){




    #if (listOfModels[i] == "weibull3"){
      #  fitLc[[i]] = fitdist(observedMargin,listOfModels[i],start = list(scale=1.92,shape=11.32,thres=0))
    #} else {
      fitLc[[i]] = fitdist(observedMargin,listOfModels[i])
      #}



    # Fit each, collect
    aicScore[i] = fitLc[[i]]$aic
    # Collect Goodness of fit score
  }
  bestFit = which.min(aicScore)

  marginGenerator = fitLc[[bestFit]]

  outputs = list()
  outputs$marginGenerator = marginGenerator
  outputs$marginDraw = drawFromMargin(marginGenerator,quantilesToEval)
  print(summary(marginGenerator))
  return(outputs)
}
