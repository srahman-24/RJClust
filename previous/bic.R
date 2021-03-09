bic=function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, 
          ...) 
{
  nparams <- nMclustParams(modelName = modelName, d = d, G = G, noise = noise, equalPro = equalPro)
  2 * loglik - nparams * log(n)
}

nMclustParams= function (modelName, d, G, noise = FALSE, equalPro = FALSE, ...) 
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  checkModelName(modelName)
  if (G == 0) {
    if (!noise) 
      stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- nVarParams(modelName, d = d, G = G) + G * d     #G * (G+1)  for G matrix 
    if (!equalPro) 
      nparams <- nparams + (G - 1)
    if (noise) 
      nparams <- nparams + 2
  }
  return(nparams)
}


nVarParams = function (modelName, d, G, ...) 
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  switch(EXPR = modelName, E = 1, V = G, EII = 1, VII = G, 
         EEI = d, VEI = G + (d - 1), EVI = 1 + G * (d - 1), VVI = G * 
           d, EEE = d * (d + 1)/2, EVE = 1 + G * (d - 1) + d * 
           (d - 1)/2, VEE = G + (d - 1) + d * (d - 1)/2, VVE = G + 
           G * (d - 1) + d * (d - 1)/2, EEV = 1 + (d - 1) + 
           G * d * (d - 1)/2, VEV = G + (d - 1) + G * d * (d - 1)/2, EVV = 1 - G + G * d * (d + 1)/2, 
            VVV = G * d * (d + 1)/2, GMAT = G + G * (G - 1) + G * (G - 1) * (G - 2)/3 + G * G + G + G * (G-1), 
            stop("invalid model name"))
}

