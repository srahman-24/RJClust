mclustBIC= function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), 
          initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), 
          Vinv = NULL, warn = mclust.options("warn"), x = NULL, verbose = interactive(), 
          ...) 
{
  dimData <- dim(data)
  oneD <- (is.null(dimData) || length(dimData[dimData > 1]) == 
             1)
  if (!oneD && length(dimData) != 2) 
    stop("data must be a vector or a matrix")
  if (oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- mclust.options("emModelNames")
        if (n <= d) {
          m <- match(modelNames, c("EII", "VII", "EEI", 
                                   "VEI", "EVI", "VVI", "GMAT"), nomatch = 0)
          modelNames <- modelNames[m]
        }
      }
    }
    if (!is.null(prior)) {
      modelNames <- setdiff(modelNames, c("EVE", "VEE", 
                                          "VVE", "EVV"))
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 
        1:9
      else 0:9
    }
    else {
      G <- sort(as.integer(unique(G)))
    }
    if (is.null(initialization$noise)) {
      if (any(G > n)) 
        G <- G[G <= n]
    }
    else {
      noise <- initialization$noise
      if (is.logical(noise)) 
        noise <- which(noise)
      if (any(match(noise, 1:n, nomatch = 0) == 0)) 
        stop("numeric or logical vector for noise must correspond to row indexes of data")
      initialization$noise <- noise
      nnoise <- length(noise)
      if (any(G > (n - nnoise))) 
        G <- G[G <= n - nnoise]
    }
    if (!is.null(initialization$subset)) {
      subset <- initialization$subset
      if (is.logical(subset)) 
        subset <- which(subset)
      initialization$subset <- subset
      if (any(G > n)) 
        G <- G[G <= n]
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    if (!missing(prior) || !missing(control) || !missing(initialization) || 
        !missing(Vinv)) 
      stop("only G and modelNames may be specified as arguments when x is supplied")
    prior <- attr(x, "prior")
    control <- attr(x, "control")
    initialization <- attr(x, "initialization")
    Vinv <- attr(x, "Vinv")
    warn <- attr(x, "warn")
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) 
      G <- Glabels
    if (is.null(modelNames)) 
      modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      out <- x[as.character(G), modelNames, drop = FALSE]
      mostattributes(out) <- attributes(x)
      attr(out, "dim") <- c(length(G), length(modelNames))
      attr(out, "dimnames") <- list(G, modelNames)
      attr(out, "G") <- as.numeric(G)
      attr(out, "modelNames") <- modelNames
      attr(out, "returnCodes") <- attr(x, "returnCodes")[as.character(G), 
                                                         modelNames, drop = FALSE]
      return(out)
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
  }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x) charmatch("E", 
                                                           x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x) charmatch("V", 
                                                           x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  if (n > .mclust$subset & is.null(initialization$subset)) {
    initialization$subset <- sample(seq.int(n), size = .mclust$subset, 
                                    replace = FALSE)
  }
  l <- length(Gall)
  m <- length(Mall)
  if (verbose) {
    cat("fitting ...\n")
    flush.console()
    pbar <- txtProgressBar(min = 0, max = l * m + 1, style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, dimnames = list(as.character(Gall), 
                                                                  as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]], dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]], dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G), modelNames, drop = FALSE]
    RET <- RET[as.character(G), modelNames, drop = FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1", ] == EMPTY]) {
        out <- mvn(modelName = mdl, data = data, prior = prior)
        BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                             n = n, d = d, G = 1, equalPro = FALSE)
        RET["1", mdl] <- attr(out, "returnCode")
        if (verbose) {
          ipbar <- ipbar + 1
          setTxtProgressBar(pbar, ipbar)
        }
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        if (verbose) {
          ipbar <- l * m + 1
          setTxtProgressBar(pbar, ipbar)
        }
        return(structure(BIC, G = G, modelNames = modelNames, 
                         prior = prior, control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD, returnCodes = RET, 
                         class = "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(data = data, modelName = mclust.options("hcModelNames")[1])
          }
          else {
            hcPairs <- hc(data = data, modelName = "EII")
          }
        }
        else {
          hcPairs <- NULL
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) 
        clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[, g]
        }
        else {
          cl <- qclass(data, as.numeric(g))
        }
        if (verbose) {
          ipbar <- ipbar + 1
          setTxtProgressBar(pbar, ipbar)
        }
        z <- unmap(cl, groups = 1:max(cl))
        if (any(apply(z, 2, max) == 0) & warn) {
          if (warn) 
            warning("there are missing groups")
          small <- sqrt(.Machine$double.neg.eps)
          z[z < small] <- small
          z <- t(apply(z, 1, function(x) x/sum(x)))
        }
        for (modelName in na.omit(modelNames[BIC[g, ] == 
                                             EMPTY])) {
          out <- me(modelName = modelName, data = data, 
                    z = z, prior = prior, control = control, 
                    warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik, n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
          if (verbose) {
            ipbar <- ipbar + 1
            setTxtProgressBar(pbar, ipbar)
          }
        }
      }
    }
    else {
      subset <- initialization$subset
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                          data = data[subset, ])
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data[subset, 
                                                         ])
          }
        }
        else {
          hcPairs <- NULL
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) 
        clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[, g]
        }
        else {
          cl <- qclass(data[subset], as.numeric(g))
        }
        if (verbose) {
          ipbar <- ipbar + 1
          setTxtProgressBar(pbar, ipbar)
        }
        z <- unmap(cl, groups = 1:max(cl))
        if (any(apply(z, 2, max) == 0) & warn) {
          if (warn) 
            warning("there are missing groups")
          small <- sqrt(.Machine$double.neg.eps)
          z[z < small] <- small
          z <- t(apply(z, 1, function(x) x/sum(x)))
        }
        for (modelName in modelNames[!is.na(BIC[g, ])]) {
          ms <- mstep(modelName = modelName, z = z, data = as.matrix(data)[initialization$subset, 
                                                                           ], prior = prior, control = control, warn = warn)
          es <- do.call("estep", c(list(data = data, 
                                        warn = warn), ms))
          out <- me(modelName = modelName, data = data, 
                    z = es$z, prior = prior, control = control, 
                    warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik, n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
          if (verbose) {
            ipbar <- ipbar + 1
            setTxtProgressBar(pbar, ipbar)
          }
        }
      }
    }
  }
  else {
    noise <- initialization$noise
    if (is.null(Vinv) || Vinv <= 0) 
      Vinv <- hypvol(data, reciprocal = TRUE)
    if (is.null(initialization$subset)) {
      if (nnoise == n) 
        stop("All observations cannot be initialised as noise!")
      if (!G[1]) {
        hood <- n * log(Vinv)
        BIC["0", ] <- 2 * hood - log(n)
        if (l == 1) {
          return(structure(BIC, G = G, modelNames = modelNames, 
                           prior = prior, control = control, initialization = list(hcPairs = hcPairs, 
                                                                                   noise = initialization$noise), warn = warn, 
                           n = n, d = d, oneD = oneD, returnCodes = RET, 
                           class = "mclustBIC"))
        }
        G <- G[-1]
        Glabels <- Glabels[-1]
      }
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                          data = data[-noise, ])
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data[-noise, 
                                                         ])
          }
        }
        else {
          hcPairs <- NULL
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) 
        clss <- hclass(hcPairs, G)
      if (verbose) {
        ipbar <- ipbar + 1
        setTxtProgressBar(pbar, ipbar)
      }
      z <- matrix(0, n, max(G) + 1)
      for (g in Glabels) {
        z[] <- 0
        k <- as.numeric(g)
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[, g]
        }
        else {
          cl <- qclass(data[!noise], k = k)
        }
        z[-noise, 1:k] <- unmap(cl, groups = 1:max(cl))
        if (any(apply(z[-noise, 1:k, drop = FALSE], 2, 
                      max) == 0) & warn) {
          if (warn) 
            warning("there are missing groups")
          z[-noise, 1:k] <- max(z[-noise, 1:k], sqrt(.Machine$double.neg.eps))
          z[-noise, 1:k] <- apply(z[-noise, 1:k, drop = FALSE], 
                                  1, function(z) z/sum(z))
        }
        z[noise, k + 1] <- 1
        K <- 1:(k + 1)
        for (modelName in na.omit(modelNames[BIC[g, ] == 
                                             EMPTY])) {
          out <- me(modelName = modelName, data = data, 
                    z = z[, K], prior = prior, Vinv = Vinv, control = control, 
                    warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik, n = n, d = d, G = k, 
                                   noise = TRUE, equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
          if (verbose) {
            ipbar <- ipbar + 1
            setTxtProgressBar(pbar, ipbar)
          }
        }
      }
    }
    else {
      subset <- initialization$subset
      subset <- setdiff(subset, noise)
      initialization$subset <- subset
      if (length(subset) == 0) 
        stop("No observations in the initial subset after removing the noise!")
      if (!G[1]) {
        hood <- n * log(Vinv)
        BIC["0", ] <- 2 * hood - log(n)
        if (l == 1) {
          return(structure(BIC, G = G, modelNames = modelNames, 
                           prior = prior, control = control, initialization = list(hcPairs = hcPairs, 
                                                                                   subset = initialization$subset), warn = warn, 
                           n = n, d = d, oneD = oneD, returnCodes = RET, 
                           class = "mclustBIC"))
        }
        G <- G[-1]
        Glabels <- Glabels[-1]
      }
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                          data = data[subset, ])
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data[subset, 
                                                         ])
          }
        }
        else {
          hcPairs <- NULL
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) 
        clss <- hclass(hcPairs, G)
      if (verbose) {
        ipbar <- ipbar + 1
        setTxtProgressBar(pbar, ipbar)
      }
      for (g in Glabels) {
        k <- as.numeric(g)
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[, g]
        }
        else {
          cl <- qclass(data[subset], k = k)
        }
        z <- unmap(cl, groups = 1:max(cl))
        if (any(apply(z, 2, max) == 0) & warn) {
          if (warn) 
            warning("there are missing groups")
          small <- sqrt(.Machine$double.neg.eps)
          z[z < small] <- small
          z <- t(apply(z, 1, function(x) x/sum(x)))
        }
        for (modelName in na.omit(modelNames[BIC[g, ] == 
                                             EMPTY])) {
          ms <- mstep(modelName = modelName, z = z, data = as.matrix(data)[subset, 
                                                                           ], prior = prior, control = control, warn = warn)
          es <- do.call("estep", c(list(data = data, 
                                        warn = warn), ms))
          if (is.na(es$loglik)) {
            BIC[g, modelName] <- NA
            RET[g, modelName] <- attr(es, "returnCode")
          }
          else {
            es$z <- cbind(es$z, 0)
            es$z[noise, ] <- matrix(c(rep(0, k), 1), 
                                    byrow = TRUE, nrow = length(noise), ncol = k + 
                                      1)
            out <- me(modelName = modelName, data = data, 
                      z = es$z, prior = prior, Vinv = Vinv, control = control, 
                      warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik, n = n, d = d, G = k, 
                                     noise = TRUE, equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
          }
          if (verbose) {
            ipbar <- ipbar + 1
            setTxtProgressBar(pbar, ipbar)
          }
        }
      }
    }
  }
  if (verbose) {
    ipbar <- l * m + 1
    setTxtProgressBar(pbar, ipbar)
  }
  if (!is.null(prior) & any(is.na(BIC))) 
    warning("The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.")
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            Vinv = Vinv, control = control, initialization = list(hcPairs = hcPairs, 
                                                                  subset = initialization$subset, noise = initialization$noise), 
            warn = warn, n = n, d = d, oneD = oneD, criterion = "BIC", 
            returnCodes = RET, class = "mclustBIC")
} 