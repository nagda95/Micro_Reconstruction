#' Conditional sampling
#'
#' conditionalSamplingPulp() takes as input a dataframe
#'
#' @param saveDirName Path to the place where previous characterization was saved
#' @param pulp_cond New data which does not contain the wall thickness or external fibrillation
#' @param pulp_raw Data used as reference
#' @param outputName Where to write the results
#' @keywords
#' @export
#' @examples
#' conditionalSamplingPulp()
#'
conditionalSamplingPulp <- function(saveDirName,pulp_cond,pulp_raw,outputName){

  # Method 1
  #start_time <- Sys.time()
  idxToRetain = length(pulp_cond[,1])
            #conditioningCombined = do.call(rbind, list(pulp_cond, pulp_raw))
  #conditioningCombined = matrix(nrow = idxToRetain,ncol = 3)
  #for(i in 1:idxToRetain ){
  #  for (j in 1:dim(pulp_cond)[2]){
  #   conditioningCombined[i,j] = rank(rbind(as.matrix(pulp_cond[i,j]), as.matrix(pulp_raw[,j])))[1]/(length(pulp_raw[,1])+1+1)
  #  }
  #}
  #end_time <- Sys.time()
  #end_time - start_time


  # Method 2
  denomNormalization = length(pulp_raw[,1])+1+1
  CC = matrix(nrow = idxToRetain,ncol = 3)
  #start_time <- Sys.time()
  for (j in 1:dim(pulp_cond)[2]){
    ap = pulp_raw[,j]
    ap3 = sort(ap)
    for(i in 1:idxToRetain ){
        CC[i,j] = min(min(which(ap3 > pulp_cond[i,j] )),length(pulp_raw[,1])+1)/denomNormalization
    }
  }
  #end_time <- Sys.time()
  #end_time - start_time




  xTemp = generatePseudoU(pulp_cond)


  conditioningData = CC#conditioningCombined
  #conditioningData = generatePseudoU(conditioningCombined)

  #conditioningData = conditioningData[1:idxToRetain,]

  print("Resume file:")
  print(saveDirName)
  load( file = paste(saveDirName,".Rdata",sep = ""))
  # Load previously determined data

  d=dim(RVM$Matrix)[1]

  condition = conditioningData[, cbind(RVM$Matrix[(d+1)-1,(d+1)-1] , RVM$Matrix[(d+1)-2,(d+1)-2] , RVM$Matrix[(d+1)-3,(d+1)-3])-2]
  Sim <- CDVineCondSim(RVM,condition)

  listOfModels = c("weibull","gamma","lnorm","norm","exp")
  print("Margin fit: Lc")
  fitLc      = fitMarginFunction(    pulp_cond$Lc, listOfModels, xTemp[,1])
  print("Margin fit: Width")
  fitWidth   = fitMarginFunction( pulp_cond$Width, listOfModels, xTemp[,2])
  print("Margin fit: Curl")
  fitCurl    = fitMarginFunction(  pulp_cond$Curl, listOfModels, xTemp[,3])

  print("Conditional margin draw: WallTkn")
  fitWallTkn$marginDraw = drawFromMargin(fitWallTkn$marginGenerator,Sim[,1])
  print("Conditional margin draw: Fibrillation")
  fitFibril$marginDraw  = drawFromMargin(fitFibril$marginGenerator,Sim[,2])


  pulp_fit = matrix(nrow=length(conditioningData[,1]),ncol=5)
  pulp_fit[,3] = fitLc$marginDraw
  pulp_fit[,4] = fitWidth$marginDraw
  pulp_fit[,5] = fitCurl$marginDraw
  pulp_fit[,1] = fitWallTkn$marginDraw
  pulp_fit[,2] = fitFibril$marginDraw


  outputNetwork = matrix(nrow = length(pulp_cond[,1]), ncol = 16)
  outputNetwork[, 1] = 0
  outputNetwork[, 2] = pulp_fit[, 3]/(1 + pulp_fit[, 5])
  outputNetwork[, 3] = pulp_fit[, 3]
  outputNetwork[, 4] = pulp_fit[, 4]
  outputNetwork[, 5] = pulp_fit[, 1]
  outputNetwork[, 6] = pulp_fit[, 5] * 100
  outputNetwork[, c(7, 8, 9, 10, 11, 12, 13, 14, 15)] = 0
  outputNetwork[, 16] = pulp_fit[, 2]

  write.table(format(outputNetwork, digits = 3), file = outputName,
              quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE,
              col.names = FALSE)






  pulp_fit = as.data.frame(pulp_fit)
  names(pulp_fit)<- c("Wall","Fibril","Lc","Width","Curl")








  characterizedPulp = list()
  characterizedPulp$fitData = pulp_fit
  return(characterizedPulp)
}
