#' Modified version of the \link[minfi]{preprocessIllumnia} normalization function
#'
#' \code{MNPpreprocessIllumina} works with single samples. Normalize to 10,000 instead of average normalization control intensity of first/reference sample
#' 
#' @param RGset object of class RGset 
#' @param bg.correct background correction 
#' @param normalize either "controls" or "no"
#' @param ref reference intensity 10,000
#' @details See \link[minfi]{preprocessIllumnia} for more details.
#' @export
#' @import minfi

MNPpreprocessIllumina <- function (rgSet, bg.correct = TRUE, normalize = c("controls", "no"),ref=10000) {
  minfi:::.isRGOrStop(rgSet)
  normalize <- match.arg(normalize)
  
  if(normalize == "controls") {
    rgSet <- MNPnormalize.illumina.control(rgSet,ref)
  }
  if(bg.correct) {
    rgSet <- minfi::bgcorrect.illumina(rgSet)
  }
 
  out <- preprocessRaw(rgSet)
  preprocess <- sprintf("Illumina_mnp, bg.correct = %s, normalize = %s, refIntensity = %d",
                        bg.correct, normalize, ref)

  ## TODO:
  ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.          
  ## packageVersion expects a string
  out@preprocessMethod <- c(rg.norm = preprocess,
                            minfi = as.character(packageVersion("minfi")),
                            manifest = as.character(packageVersion(minfi:::.getManifestString(rgSet@annotation))))
  #packageVersion(getManifest(rgSet)))
  out
}

#' @describeIn MNPpreprocessIllumina Modified version of \code{normalize.illumina.control} called by\code{\link{MNPpreprocessIllumina}}
#' 
#' \code{MNPnormalize.illumina.control} works with single samples. Normalize to 10,000 instead of average normalization control intensity of first/reference sample
#' Called by \code{\link{MNPpreprocessIllumina}}
#' 
#' @param rgSet an object of class RGset
#' @param ref reference intensity 10,000
#' @export
#' @import minfi

MNPnormalize.illumina.control <- function (rgSet, ref=10000) {
  ## This function returns an rgset, not a methylset
  ## code duplication
  Green <- getGreen(rgSet)
  Red <- getRed(rgSet)    
  
  if (minfi:::.is450k(rgSet) || minfi:::.isEPIC(rgSet)) {
    AT.controls <- getControlAddress(rgSet, controlType = c("NORM_A", "NORM_T"))
    CG.controls <- getControlAddress(rgSet, controlType = c("NORM_C", "NORM_G"))
  }
  if (minfi:::.is27k(rgSet)) {
    AT.controls <- getControlAddress(rgSet, controlType = "Normalization-Red")
    CG.controls <- getControlAddress(rgSet, controlType = "Normalization-Green")
  }
  Green.avg <- colMeans(Green[CG.controls, , drop = FALSE])
  Red.avg <- colMeans(Red[AT.controls, , drop = FALSE])
  Green.factor <- ref/Green.avg
  Red.factor <- ref/Red.avg
  Green <- sweep(Green, 2, FUN = "*", Green.factor)
  Red <- sweep(Red, 2, FUN = "*", Red.factor)
  assay(rgSet, "Green") <- Green
  assay(rgSet, "Red") <- Red
  rgSet
}



