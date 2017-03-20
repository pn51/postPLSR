
library(pls)
#library(pracma)
#library(R.utils)
#library(parallel)

#' @title .VIPh
#' @description .VIPh
#' @keywords internal
#' @param plsOutput model output of mvr/pls-related functions of the pls package
#' @param h number of components
#' @return VIP scores for all predictors of the PLS model
.VIPh <- function(plsOutput, h) {
  SS <- c(plsOutput$Yloadings)^2 * colSums(plsOutput$scores^2)
  Wnorm2 <- colSums(plsOutput$loading.weights^2)
  SSW <- sweep(plsOutput$loading.weights^2, 2, SS / Wnorm2, "*")
  if (length(SS) > 1) {
    sqrt(nrow(SSW) * apply(SSW, 1, cumsum)[h,] / cumsum(SS)[h])
  } else{
    sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / SS)
  }
}


#' @title .VIPjh
#' @description .VIPjh
#' @keywords internal
#' @param plsOutput model output of mvr/pls-related functions of the pls package
#' @param j predictor index
#' @param h number of components
#' @return VIP score for predictor \code{j} of the PLS model
.VIPjh <- function(plsOutput, j, h) {
  b <- c(plsOutput$Yloadings)[1:h]
  T <- plsOutput$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- plsOutput$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}


#########################################################################
#########################################################################

#' @title .plsrReconstruction
#' @description .plsrReconstruction
#' @keywords internal
#' @param xData predictor expression matrix (samples x genes)
#' @param yData response expression matrix (samples x genes)
#' @param components number of PLSR components to use
#' @param loops if \code{TRUE}, include loops in the regression models
#' @import pls
#' @import R.utils
#' @return A matrix of VIP scores, where $v_{ij}$ is the score corresponding to edge \eqn{i \rightarrow j}
.plsrReconstruction <- function(xData, yData, components=3, loops=FALSE) {

  m <- nrow(xData)
  n <- ncol(yData)

  # compute the VIP scores for each target gene; v_{ji} corresponds to edge (j,i)
  vipMatrix <- sapply(seq(n), function(ii) {
    if (!loops) {
      res <- mvr(yData[, ii, drop = FALSE] ~ xData[, -ii], ncomp = min(n - 1, m, components), method = 'kernelpls', scale = FALSE)
      R.utils::insert(.VIPh(res, res$ncomp), ii, 0)
    } else {
      res <- mvr(yData[, ii, drop = FALSE] ~ xData, ncomp = min(n, m, components), method = 'kernelpls', scale = FALSE)
      .VIPh(res, res$ncomp)
    }
  })

  dimnames(vipMatrix) <- list(colnames(yData), colnames(yData))
  vipMatrix
}

#' @title .plsrReconstructionTarget
#' @description .plsrReconstructionTarget
#' @keywords internal
#' @param xData predictor expression matrix (samples x genes)
#' @param yData response or target expression vector
#' @param source predictor gene of interest
#' @param target response gene of interest
#' @param components number of PLSR components to use
#' @param loops if \code{TRUE}, include loops in the regression models
#' @import pls
#' @return VIP score for the \code{source}-\code{target} gene pair.
.plsrReconstructionTarget <- function(xData, yData, source, target, components=3, loops=FALSE) {


  # compute the VIP score for edge (source, target); assume xData also contains data for the target and yData only contains only data for the target
  if (source == target & !loops) {
    0
  } else{
    maxComp <- min(ncol(xData) - !loops, nrow(xData), components)
    if (!loops) {
      plsrRes <- mvr(yData ~ xData[, -target], ncomp = maxComp, method = 'kernelpls', scale = FALSE)
    } else{
      plsrRes <- mvr(yData ~ xData, ncomp = maxComp, method = 'kernelpls', scale = FALSE)
    }
  }

  vipScore <- .VIPjh(plsrRes,
                     source - ifelse(loops, 0, ifelse(source < target, 0, 1)),
                     plsrRes$ncomp)
  ifelse(is.na(vipScore), 0, vipScore)

}

#########################################################################
#########################################################################

#' @title .generateExprData
#' @description .generateExprData
#' @keywords internal
#' @param coefMatrix edge weight matrix
#' @param x0 initial condition
#' @param type gene expression model to use; 'expr' or 'deltaExpr'
#' @param timepoints number of timepoints
#' @param sigma noise standard deviation
#' @return A list of two generated expression data matrices \code{xData}, \code{yData} containing samples corresponding to the inputs and outputs of the model \eqn{x_i(t+1) = \sum_{j} w_{ji} x}
.generateExprData <- function(coefMatrix, x0, type='expr', timepoints=10, sigma=.1) {

  n <- length(x0)

  # generate data with one of the two models
  if (type == 'expr') {
    xData <- do.call(rbind, Reduce(function(acc, right) {
      c(coefMatrix %*% acc) + rnorm(n, sd = sigma)
    }, seq(timepoints - 1), x0, accumulate = TRUE))
    xData <- scale(xData)
    yData <- xData[-1, ]
    xData <- xData[-nrow(xData), ]
    #colnames(xData) <- rownames(coefMatrix)
    colnames(yData) <- rownames(coefMatrix)
    list(xData = xData, yData = yData)
  } else if (type == 'deltaExpr') {
    xData <- do.call(rbind, Reduce(function(acc, right) {
      c(coefMatrix %*% acc + acc) + rnorm(n, sd = sigma)
    }, seq(timepoints - 1), x0, accumulate = TRUE))
    xData <- scale(xData)
    yData <- diff(xData)
    xData <- xData[-nrow(xData),]
    #colnames(xData) <- rownames(coefMatrix)
    colnames(yData) <- rownames(coefMatrix)
    list(xData = xData, yData = yData)
  } else{
    stop('invalid type')
  }

}

#' @title .simulateVIP
#' @description .simulateVIP
#' @keywords internal
#' @param calculatedVIP matrix of data-derived VIP scores, calculated with \code{.plsrReconstruction}
#' @param adjMatrix prior adjacency matrix
#' @param trueAdjMatrix true adjacency matrix
#' @param components number of PLSR components to use
#' @param loops if TRUE, include loops in the regression models
#' @param type gene expression model to use; 'expr' or 'deltaExpr'
#' @param timepoints number of timepoints to use in the generated expression data
#' @param iterations number of VIP scores to generate
#' @param mu mean edge weight
#' @param sigmaX edge weight standard deviation
#' @param sigmaW initial condition standard deviation
#' @param sigmaE generated data noise
#' @param calcPosterior if \code{TRUE}, calculated the posterior probability
#' @param p prior probability of an edge given an edge in the input network
#' @param q prior probability of an edge given a non-edge in the input network
#' @param cores number of cores for parallelization
#' @import parallel
#' @return a data frame containing the data-derived VIP scores, estimated number of null VIP scores that are less than the data-derived score in the edge and non-edge cases, and (optionally) posterior edge probabilities for each ordered pair of genes.
.simulateVIP <- function(calculatedVIP, adjMatrix, trueAdjMatrix=adjMatrix,
                         components=3, loops=FALSE,
                         type='expr',  timepoints=10, iterations=1000,
                         mu=1, sigmaX=0.1, sigmaW=0.1, sigmaE=0.1,
                         calcPosterior=TRUE, p=0.5, q=0.5,
                         cores=1) {

  # adjMatrix <- abs(adjMatrix) > 0
  adjMatrix <- adjMatrix > 0

  if (loops) {
    diag(adjMatrix) <- FALSE
  }

  n <- sum(adjMatrix) #number of edges.
  m <- nrow(adjMatrix)

  ###############

  # first, generate null VIP scores for the original network and fraction of times is it
  # less than the data-derived score.

  coefMatrix <- t(adjMatrix)

  df <- parallel::mclapply(seq(iterations), mc.cores = cores, FUN = function(jj) {
    coefMatrix[coefMatrix] <- rnorm(n, mu, sigmaW) * sample(c(-1, 1), n, TRUE)

    with(.generateExprData(coefMatrix, rnorm(m, sigmaX), type, timepoints, sigmaE),
         .plsrReconstruction(xData, yData, components, loops) <= calculatedVIP)
  })
  df <- Reduce('+', df) / iterations
  df <- data.frame(edgeNumber = seq(m^2), edge = c(adjMatrix),
                   vip = c(calculatedVIP), originalFraction = c(df))
  if (is.matrix(trueAdjMatrix) && (dim(adjMatrix) == dim(trueAdjMatrix))) {
    df$trueEdge <- c(trueAdjMatrix)
  }

  if (!loops) {
    df <- subset(df, !(edgeNumber %in% (1 + (m + 1)*seq(0, m - 1))))
  }

  ###############

  # now iterate through the edges, swap the edge type, and calcuate for just that edge
  df$swappedFraction <- unlist(parallel::mclapply(df$edgeNumber, mc.cores = cores, FUN = function(ii) {
    curEdge <- adjMatrix[ii]
    adjMatrix[ii] <- !curEdge
    coefMatrix <- t(adjMatrix)

    rowIndex <- (ii - 1) %% m + 1 #source
    colIndex <- (ii - 1) %/% m + 1 #target

    edgeP <- sapply(seq(iterations), function(jj) {
      coefMatrix[coefMatrix] <- rnorm(n + ifelse(curEdge, -1, 1), mu, sigmaW) *
        sample(c(-1, 1), n + ifelse(curEdge, -1, 1), TRUE)

      #find VIP scores
      with(.generateExprData(coefMatrix, rnorm(m, sigmaX), type, timepoints, sigmaE),
           .plsrReconstructionTarget(xData, yData[, colIndex, drop = FALSE],
                                     rowIndex, colIndex, components, loops) <= calculatedVIP[ii])
    })

    sum(edgeP) / iterations

  }))

  # calculate the posterior probability
  if (calcPosterior) {
    calculatePosterior(df, p, q, TRUE)
  } else{
    df
  }

}

#' @title calculatePosterior
#' @description Computes posterior edge probabilities from the simulation results
#' @param simulateVIPResults output of simulateVIP
#' @param p prior probability of an edge given an edge in the input network
#' @param q prior probability of an edge given a non-edge in the input network
#' @param attachResults if \code{TRUE}, add/replace the \code{posterior} column of \code{simulateVIPResults} with the computed posterior probabilities
#' @return a data frame with the computed posterior probablities attached, or a vector of posterior probabilities
#' @export
calculatePosterior <- function(simulateVIPResults, p=0.5, q=0.5, attachResults=TRUE) {

  edgeInd <- simulateVIPResults$edge > 0
  nonEdgeInd <- !edgeInd

  posterior <- rep(0,length(edgeInd))

  posterior[edgeInd] <- simulateVIPResults$originalFraction[edgeInd] * p /
    (simulateVIPResults$originalFraction[edgeInd] * p +
       (1 - simulateVIPResults$swappedFraction[edgeInd]) * (1 - p))
  posterior[nonEdgeInd] <- simulateVIPResults$swappedFraction[nonEdgeInd] * q /
    (simulateVIPResults$swappedFraction[nonEdgeInd] * q +
       (1 - simulateVIPResults$originalFraction[nonEdgeInd]) * (1 - q))

  posterior[edgeInd & is.na(posterior)] <- p
  posterior[nonEdgeInd & is.na(posterior)] <- q

  if (attachResults) {
    simulateVIPResults$posterior <- posterior
    simulateVIPResults
  } else{
    posterior
  }
}


##Returns the AUC for the input network or true network
#posteriorRes: output of simulateVIP with calculatePosterior==TRUE, or output of
#calculatePosterior with attachResults=TRUE
#loops: if FALSE, exclude loops (if any were computed) from the AUC calculation
#useTrueEdge: if
#' @title posteriorAUC
#' @description Computes an AUC for the posterior edge probabilities
#' @param posteriorRes output of \code{simulateVIP} function
#' @param loops if \code{FALSE}, exclude loops (if any were computed) from the AUC calculation
#' @param useTrueEdge if \code{TRUE} and the an additional adjacency matrix was supplied to \code{simulateVIP}, compute the AUC based on this matrix; otherwise, compute based on the prior network
#' @import pracma
#' @return an AUC value
#' @export
posteriorAUC <- function(posteriorRes, loops=FALSE, useTrueEdge=FALSE) {

  # sort posterior probabities in descending order for the ROC curve
  posteriorRes <- posteriorRes[order(posteriorRes$posterior, decreasing = TRUE), ]

  if (loops) {
    m <- ceiling(sqrt(nrow(posteriorRes)))
    posteriorRes <- subset(posteriorRes, !(edgeNumber %in% (1 + (m + 1) * seq(0, m - 1))))
  }

  if (useTrueEdge) {
    tpr <- cumsum(posteriorRes$trueEdge)
    fpr <- cumsum(!posteriorRes$trueEdge)
  } else{
    tpr <- cumsum(posteriorRes$edge)
    fpr <- cumsum(!posteriorRes$edge)
  }
  tpr <- tpr / tail(tpr, 1)
  fpr <- fpr / max(tail(fpr, 1), 1)

  ind <- !duplicated(posteriorRes$posterior, fromLast = TRUE)
  tpr <- c(0, tpr[ind])
  fpr <- c(0, fpr[ind])

  auc <- ifelse(all(fpr == 0), 1, pracma::trapz(fpr, tpr))

  #list(auc=auc,fpr=fpr,tpr=tpr)
  auc
}

###################################
###################################

#' @title .convertDataExpr
#' @description .convertDataExpr
#' @keywords internal
#' @param exprData gene expression matrix (timepoints x genes)
#' @param rescale if \code{TRUE}, rescale the expression data for each gene to mean 0 and standard deviation 1
#' @return A list of two matrices \code{xData}, \code{yData} of samples corresponding to the inputs and outputs of the model \eqn{x_i(t+1) = \sum_{j} w_{ji} x}
.convertDataExpr <- function(exprData, rescale=TRUE) {
  if (rescale) {
    xData <- scale(exprData)
  }
  yData <- xData[-1, ]
  xData <- xData[-nrow(xData), ]
  #colnames(xData) <- colnames(exprData)
  colnames(yData) <- colnames(exprData)
  list(xData = xData, yData = yData)
}

#' @title .convertDataDeltaExpr
#' @description .convertDataDeltaExpr
#' @keywords internal
#' @param exprData expression matrix (timepoints x genes)
#' @param rescale if \code{TRUE}, rescale the expression data for each gene to mean 0 and standard deviation 1
#' @return A list of two matrices \code{xData}, \code{yData} of samples corresponding to the inputs and outputs of the model \eqn{\Delta x_i(t+1) = \sum_{j} w_{ji} x}
.convertDataDeltaExpr <- function(exprData, rescale=TRUE) {
  if (rescale) {
    xData <- scale(exprData)
  }
  yData <- diff(xData)
  xData <- xData[-nrow(xData), ]
  colnames(xData) <- colnames(exprData)
  colnames(yData) <- colnames(exprData)
  list(xData = xData, yData = yData)
}

###################################
###################################

#' @title posteriorPLSRReconstruction
#' @description Computes null distrbution-related results and computes posterior edge probabilities using a prior network and input expression data
#' @param exprData expression matrix with rows corresponding to timepoints in increasing temporal order and columns to genes; or a list containing two matrices \code{xData} (inputs) and \code{yData} (responses) with rows corresponding to samples and columns to genes
#' @param adjMatrix adjacency matrix
#' @param type gene expression model to use; 'expr' or 'deltaExpr'
#' @param components number of PLSR components to use
#' @param timepoints number of timepoints to use in the generated expression data
#' @param iterations number of VIP scores to generate
#' @param mu mean edge weight
#' @param sigmaX edge weight standard deviation
#' @param sigmaW initial condition standard deviation
#' @param sigmaE generated data noise
#' @param p probability of an edge, given that the pair was a prior edge
#' @param q probability of an edge, given that the pair was a prior non-edge
#' @param cores number of cores for parallelization
#' @return a data frame containing the data-derived VIP scores, estimated number of null VIP scores that are less than the data-derived score in the edge and non-edge cases, and (optionally) posterior edge probabilities for each ordered pair of genes.
#' @export
posteriorPLSRReconstruction <- function(exprData, adjMatrix, type='expr', components=3,
                                        timepoints=8, iterations=1000, mu=1, sigmaX=.1, sigmaW=.1, sigmaE=.1,
                                        p=0.5, q=0.5, cores=1) {

  #check data input
  if (!setequal(colnames(adjMatrix), rownames(adjMatrix)) ||
      any(duplicated(rownames(adjMatrix))) || any(duplicated(colnames(adjMatrix)))) {
    stop('invalid adjacency matrix gene names', call. = FALSE)
  }

  if (is.matrix(exprData)) {
    if (!setequal(colnames(exprData), rownames(adjMatrix)) ||
        any(duplicated(colnames(exprData)))) {
      stop('invalid expression matrix gene names', call. = FALSE)
    }

    if (type == 'expr') {
      exprData <- .convertDataExpr(exprData, TRUE)
    } else if (type == 'deltaExpr') {
      exprData <- .convertDataDeltaExpr(exprData, TRUE)
    } else{
      warning('invalid type; using "expr"', call. = FALSE)
      type <- 'expr'
      exprData <- .convertDataExpr(exprData, TRUE)
    }
  } else{
    if (!is.list(exprData) || !('xData' %in% names(exprData)) || !('yData' %in% names(exprData)) ||
        !is.matrix(exprData$xData) || !is.matrix(exprData$yData) ||
        !all(dim(exprData$xData) == dim(exprData$yData))) {
      stop('invalid input data format', call. = FALSE)
    }

    if (!setequal(colnames(exprData$xData), rownames(adjMatrix)) ||
        !setequal(colnames(exprData$xData), colnames(exprData$yData)) ||
        any(duplicated(colnames(exprData$xData))) || any(duplicated(colnames(exprData$yData)))) {
      stop('invalid input data gene names', call. = FALSE)
    }

    if (!(type %in% c('expr', 'deltaExpr'))) {
      warning('invalid type; using "expr"', call. = FALSE)
      type <- 'expr'
    }

    exprData$yData <- exprData$yData[,colnames(exprData$xData)]
  }

  adjMatrix <- adjMatrix[colnames(exprData$yData), colnames(exprData$yData)]

  vipMatrix <- with(exprData, .plsrReconstruction(xData, yData, components))
  simulateVIPResults <- .simulateVIP(vipMatrix, adjMatrix,
                                     components = components, type = type, timepoints = timepoints, iterations = iterations,
                                     mu = mu, sigmaX = sigmaX, sigmaW = sigmaW, sigmaE = sigmaE,
                                     calcPosterior = TRUE, p = p, q = q, cores = cores)

  simulateVIPResults

}
