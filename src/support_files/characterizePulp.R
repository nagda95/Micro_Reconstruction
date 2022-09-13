#' Model selection for pulp characterization using C-Vine copulas
#'
#' characterizePulp() takes as input a dataframe
#'
#' @param pulpPath Path to the fiber characterization file
#' @param saveDirName Path to the place where files are saved
#' @keywords
#' @export
#' @examples
#' characterizePulp()
#'characterizePulp
characterizePulp <- function(pulpPath,saveDirName){
  pulpPath = gsub("[\\]","/",pulpPath)
  saveDirName = gsub("[\\]","/",saveDirName)

  pulp_raw = readAndFormatPulpFile(pulpPath)

  udata = generatePseudoU(pulp_raw)
  # Generate pseudo-observations



  ##########################################################################################
  condVars = 3
  copType = "CVine"
  selMethod = "AIC"
  checkIndep = TRUE
  #fitSet = c(0,1)#NA
  #fitSet = c(0,1,3,4,5,104,204)#NA
  fitSet = c(0,1,3,4,5)#NA

  print("CDVineCondFit input options:")
  print(paste("Nx =",condVars))
  print(paste("type =",copType))
  print(paste("selectioncrit =",selMethod))
  print(paste("indeptest =",checkIndep))
  print(paste("familyset =",fitSet))

  RVM <- CDVineCondFit(udata,Nx=condVars, type=copType, selectioncrit=selMethod,
                       indeptest=checkIndep, level=0.05, familyset = fitSet, rotations = TRUE)
  # Fit Copula
  #RVM <- RVineStructureSelect(udata, type=copType, selectioncrit=selMethod,
  #                     indeptest=checkIndep, level=0.05)

  Sim <- RVineSim(1*dim(pulp_raw)[1],RVM)

  print("******************** Margin fits ********************")
  listOfModels = c("weibull","gamma","lnorm","norm","exp")
  #listOfModels = c("weibull3","gamma","lnorm","norm","exp")
  print("Margin fit: Lc")
  fitLc      = fitMarginFunction(    pulp_raw$Lc, listOfModels, Sim[,3])
  print("Margin fit: Width")
  fitWidth   = fitMarginFunction( pulp_raw$Width, listOfModels, Sim[,4])
  print("Margin fit: Curl")
  fitCurl    = fitMarginFunction(  pulp_raw$Curl, listOfModels, Sim[,5])
  print("Margin fit: WallTkn")
  fitWallTkn = fitMarginFunction(  pulp_raw$Wall, listOfModels, Sim[,1])
  print("Margin fit: Fibrillation")
  fitFibril  = fitMarginFunction(pulp_raw$Fibril, listOfModels, Sim[,2])



  #save(list = "RVM",file = paste(saveDirName,".Rdata",sep = ""))
  # Save the RVM for future use.

  save(RVM,fitLc,fitWidth,fitCurl,fitWallTkn,fitFibril, file = paste(saveDirName,".Rdata",sep = ""))
  # Here we should allow for the files to be saved as well.

  characterizedPulp = list()
  characterizedPulp$rawData = pulp_raw
  return(characterizedPulp)
}




