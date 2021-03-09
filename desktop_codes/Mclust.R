Mclust = function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), 
          initialization = NULL, warn = mclust.options("warn"), x = NULL, 
          verbose = interactive(), ...) 
{
  call <- match.call()
  data <- data.matrix(data)
  if (!is.null(x)) 
    if (!inherits(x, "mclustBIC")) 
      stop("If provided, argument x must be an object of class 'mclustBIC'.")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustBIC")
  mc[[2]] <- data
  Bic <- eval(mc, parent.frame())
  G <- attr(Bic, "G")
  modelNames <- attr(Bic, "modelNames")
  Sumry <- summary(Bic, data, G = G, modelNames = modelNames)
  if (length(Sumry) == 0) 
    return()
  if (!(length(G) == 1)) {
    bestG <- length(tabulate(Sumry$cl))
    if (warn) {
      if (bestG == max(G) & warn) 
        warning("optimal number of clusters occurs at max choice")
      else if (bestG == min(G) & warn) 
        warning("optimal number of clusters occurs at min choice")
    }
  }
  oldClass(Sumry) <- NULL
  Sumry$bic <- Sumry$bic[1]
  Sumry$hypvol <- if (is.null(attr(Bic, "Vinv"))) 
    as.double(NA)
  else 1/attr(Bic, "Vinv")
  df <- if (is.null(Sumry$modelName)) 
    NULL
  else with(Sumry, nMclustParams(modelName, d, G, noise = (!is.na(hypvol)), 
                                 equalPro = attr(Sumry, "control")$equalPro))
  ans <- c(list(call = call, data = data, BIC = Bic, df = df), 
           Sumry)
  orderedNames <- c("call", "data", "modelName", "n", "d", 
                    "G", "BIC", "bic", "loglik", "df", "hypvol", "parameters", 
                    "z", "classification", "uncertainty")
  structure(ans[orderedNames], class = "Mclust")
}