#' Margin sampler
#'
#' This function allows you to express your love of cats.
#' @param marginGenerator Contains fitted univariate margin parameters
#' @param quantilesToEval Contains quantiles to invert to obtain value
#' @keywords cats
#' @export
#' @examples
#' drawFromMargin()
#'
drawFromMargin <- function(marginGenerator,quantilesToEval){
  if (marginGenerator$distname == "weibull"){
    marginDraw = qweibull(quantilesToEval,shape = marginGenerator$estimate[1], scale = marginGenerator$estimate[2])
  } else if (marginGenerator$distname == "gamma"){
    marginDraw = qgamma(quantilesToEval, shape = marginGenerator$estimate[1], rate = marginGenerator$estimate[2])
  } else if (marginGenerator$distname == "lnorm"){
    marginDraw = qlnorm(quantilesToEval, meanlog = marginGenerator$estimate[1], sdlog = marginGenerator$estimate[2])
  } else if (marginGenerator$distname == "norm"){
    marginDraw = qnorm(quantilesToEval, mean = marginGenerator$estimate[1], sd = marginGenerator$estimate[2])
  } else if (marginGenerator$distname == "exp"){
    marginDraw = qexp(quantilesToEval, rate = marginGenerator$estimate[1])
  } else if (marginGenerator$distname == "weibull3"){
    marginDraw = qweibull3(quantilesToEval,marginGenerator$estimate[1],marginGenerator$estimate[2],thres=marginGenerator$estimate[3])
  }
  return(marginDraw)
}
