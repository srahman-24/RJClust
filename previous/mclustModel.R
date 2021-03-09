mclustModel <- function(data, BICvalues, G=NULL, modelNames=NULL, ...)
{
  mc = match.call(expand.dots = FALSE)       # I don't know what is this fucking match.call()
  if (is.null(attr(BICvalues,"initialization")$noise)) {
    mc[[1]] = as.name("summaryMclustBIC")
  }
  else {
    mc[[1]] = as.name("summaryMclustBICn")
  }
  nm = names(mc)
  mc[1:3] = mc[c(1,3,2)]
  nm[1:3] = nm[c(1,3,2)]
  nm[nm == "BICvalues"] <- "object" 
  names(mc) <- nm
  ans <- eval(mc, parent.frame())
  ans$classification <- ans$uncertainty <- NULL
  attr( ans, "bestBICvalues") <- NULL
  attr( ans, "prior") <- NULL
  attr( ans, "control") <- NULL
  attr( ans, "initialization") <- NULL
  oldClass(ans) <- "mclustModel"
  ans
}

mclustModelNames <- function(model)
{
  type <- switch(EXPR = as.character(model),
                 "E" = "univariate, equal variance",
                 "V" = "univariate, unequal variance",
                 "EII" = "spherical, equal volume",
                 "VII" = "spherical, varying volume",
                 "EEI" = "diagonal, equal volume and shape",
                 "VEI" = "diagonal, equal shape",
                 "EVI" = "diagonal, equal volume, varying shape",
                 "VVI" = "diagonal, varying volume and shape",
                 "GMAT"= "G Matrix Covariance",                   # I inlcuded this.
                 "EEE" = "ellipsoidal, equal volume, shape and orientation",
                 "EVE" = "ellipsoidal, equal volume and orientation",
                 "VEE" = "ellipsoidal, equal shape and orientation",
                 "VVE" = "ellipsoidal, equal orientation",
                 "EEV" = "ellipsoidal, equal volume and shape",
                 "VEV" = "ellipsoidal, equal shape",
                 "EVV" = "ellipsoidal, equal volume",
                 "VVV" = "ellipsoidal, varying volume, shape, and orientation",
                 "X"   = "univariate normal",
                 "XII" = "spherical multivariate normal",
                 "XXI" = "diagonal multivariate normal",
                 "XXX" = "ellipsoidal multivariate normal",
                 warning("invalid model"))
  return(list(model = model, type = type))
}