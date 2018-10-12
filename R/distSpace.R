distSpace <- function(trainingData,
                      testData = NULL,
                      type = "bagdistance",
                      options = NULL){

  ######
  # Check input.
  if (missing(trainingData)) {
    stop("Input argument trainingData is required.")
  }
  if (!is.list(trainingData)) {
    stop(paste("The input argument trainingData should be a list of",
               "matrices for multivariate data, or a list of arrays",
               "for functional data.")
         )
  }

  # Get number of groups in training data.
  nGroups <- length(trainingData)
  if (nGroups < 2) {
    stop("There should be at least two groups in the training data.")
  }

  dataType <- NA
  # browser()
  # Verify the training data.
  dimP <- NA
  for (i in 1:nGroups) {
    currentData <- trainingData[[i]]
    if (is.data.frame(currentData)) {
      currentData <- as.matrix(currentData)
      trainingData[[i]] <- currentData
    }
    if (i == 1 & !is.matrix(currentData)) { 
      dataType <- "functional"
    } 
    if (i == 1 & is.matrix(currentData)) {
      dataType <- "multivariate"
    }
    if (dataType == "multivariate") {
      if (!is.matrix(currentData)) {
        stop(paste("The input argument trainingData should be a list of matrices.",
                   " Argument ", i, " is not a matrix.", sep = "")
        )
      }
      if (!is.numeric(currentData)) {
        stop(paste("The input argument trainingData should be a list of matrices.",
                   " Argument ", i, " is not a numeric matrix.", sep = "")
        )
      }
      dims <- NULL
      dimsTime <- NULL
      if (i == 1) dimP <- ncol(currentData)
      nObs <- nrow(currentData)
      nCol <- ncol(currentData)
    }
    if (dataType == "functional") {
      if (!is.array(currentData))  {
        stop(paste("The input argument trainingData should be a list of arrays.",
                   " Argument ", i, " is not an array.", sep = "")
        )
      }
      if (!is.numeric(currentData))  {
        stop(paste("The input argument trainingData should be a list of arrays.",
                   " Argument ", i, " is not a numeric array.", sep = "")
        )
      }
      dims <- dim(currentData)
      if (i == 1) {
        dimsTime <- dims[1:(length(dims) - 2)] 
        dimP <- dims[length(dims)] 
      }
      nObs <- dims[length(dims) - 1]
      nCol <- dims[length(dims)]
    }
    if (nCol != dimP) {
      stop(paste("The number of variables should be the same across all groups.",
                 " Group ", i, " does not have the same number of variables as",
                 " group 1.", sep = "")
      )
    }
    if (sum(dims[1:(length(dims) - 2)] != dimsTime) > 1) {
      stop(paste("The dimensions of the domain should be the same across all groups.",
                 " Group ", i, " does not have the same domain dimensions as",
                 " group 1.", sep = "")
      )
    }
    if (nObs > sum(complete.cases(currentData))) {
      stop(paste("Missing values in the training data are not allowed.",
                 " Group ", i, " has missing cases.", sep = ""))
    }
    if (nObs <= nCol) {
      stop(paste("Exact fit situations in the training data are not allowed.",
                 " Group ", i, " does not have at least p_i + 1 observations.",
                 sep = ""))
    }
    rm(currentData)
  }


  # Now verify the test data.
  if (!is.null(testData)) {
    if (dataType == "multivariate") {
      if (!is.matrix(testData)) {
        stop("The input argument testData should be a matrix.")
      }
      if (!is.numeric(testData)) {
        stop(paste("The input argument testData should be a numeric matrix.")
        )
      }
      nObsTest <- nrow(testData)
      nColTest <- ncol(testData)
      if (nColTest != dimP) {
        stop(paste("The number of variables should be the same across all groups.",
                   " The testing group does not have the same dimensions as",
                   " group 1 in the training data.", sep = "")
        )
      }
    }
    if (dataType == "functional") {
      if (!is.array(testData)) {
        stop("The input argument testData should be an array.")
      }
      if (!is.numeric(testData)) {
        stop(paste("The input argument testData should be a numeric array.")
        )
      }
      dims <- dim(testData)
      nObsTest <- dims[length(dims) - 1]
      nColTest <- dims[length(dims)]
      if (nColTest != dimP) {
        stop(paste("The number of variables should be the same across all groups.",
                   " The testing group does not have the same dimensions as",
                   " group 1 in the training data.", sep = "")
        )
      }
      if (sum(dims[1:(length(dims) - 2)] != dimsTime) > 1) {
        stop(paste("The dimensions of the domain should be the same across all groups.",
                   " The testing group does not have the same domain dimensions as",
                   " group 1.", sep = "")
        )
      }
    }
    if (nObsTest > sum(complete.cases(testData))) {
      stop(paste("Missing values in the test data are not allowed."))
    }
  } else {
    nObsTest <- 0
  }

  # Check type
  Indtype <- match(type, c("bagdistance", "outlyingness", "adjOutl",
                           "fbd", "fSDO", "fAO"))[1]
  if (is.na(Indtype)) {
    stop(paste("The input argument type should be one of bagdistance",
               "outlyingness or adjOutl for multivariate data and one",
               "of fbd, fSDO or fAO for functional data")
         )
  }
  if (dataType == "multivariate" & Indtype > 3) {
    stop(paste("The input argument type should be one of bagdistance",
               "outlyingness or adjOutl for multivariate data.")
         )
  }
  if (dataType == "functional" & Indtype < 4) {
    stop(paste("The input argument type should be one of fbd,", 
               "fSDO or fAO for functional data")
         )
  }
  
  # Check options
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("Input argument options should be a list.")
  }
  if (!("alpha" %in% names(options))) {
    options$alpha = 0
  }
  if (!("time" %in% names(options))) {
    options$time = NULL
  }

  # Store a vector containing the number of observations in each group.
  nFuncVector <- rep(0, nGroups + 1)
  for (i in 1:nGroups) {
    dims <- dim(trainingData[[i]])
    nFuncVector[i] <- dims[length(dims) - 1]
  }
  if (nObsTest > 0) {
    nFuncVector[nGroups + 1] <- nObsTest
  }

  #######################################################################
  #######################################################################
  
  
  # Merge all data into one structure
  if (dataType == "multivariate") {
    dataToCalcDist <- matrix(0.0,
                             nrow = sum(nFuncVector),
                             ncol = dimP)
    for (i in nGroups:2) {
      dataToCalcDist[(sum(nFuncVector[1:(i - 1)]) + 1):sum(nFuncVector[1:i]), ] <-
        trainingData[[i]]
    }
    dataToCalcDist[1:sum(nFuncVector[1]), ] <-
      trainingData[[1]]
    if (nObsTest > 1) {
      dataToCalcDist[(sum(nFuncVector[1:(nGroups)]) + 1):sum(nFuncVector), ] <-
        testData
    }
  }
  if (dataType == "functional") {
    dataToCalcDist <- trainingData[[1]]
    for (i in 2:length(trainingData)) {
      dataToCalcDist <- abind(dataToCalcDist, trainingData[[i]], along = length(dims) - 1)
    }
    if (nObsTest > 1) {
      dataToCalcDist <- abind(dataToCalcDist, testData, along = length(dims) - 1)
    }
  }

  # Now calculate distances to all training groups
  resultDistance <- matrix(0.0,
                           nrow = sum(nFuncVector),
                           ncol = nGroups + 1)
  resultDistance <- data.frame(resultDistance)
  resultDistance[, nGroups + 1] <- rep(c(paste("TR", 1:nGroups, sep = ""),
                                         paste("TE")),
                                       nFuncVector)
  for (i in 1:nGroups) {
    
    if (dataType == "multivariate") {
      if (type == "bagdistance") {
        tempResult <- bagdistance(x = trainingData[[i]],
                                  z = dataToCalcDist,
                                  options = options)
        if (is.null(tempResult$bagdistance)) {
          warning(paste("Bagdistance to training group", i,
                        "could not be computed.")
          )
        } else {
          resultDistance[, i] <- tempResult$bagdistance
        }
      }
      if (type == "outlyingness") {
        tempResult <- outlyingness(x = trainingData[[i]],
                                   z = dataToCalcDist,
                                   options = options)
        if (is.null(tempResult$outlyingnessZ)) {
          warning(paste("Outlyingness to training group", i,
                        "could not be computed.")
          )
        } else {
          resultDistance[, i] <- tempResult$outlyingnessZ
        }
      }
      if (type == "adjOutl") {
        tempResult <- adjOutl(x = trainingData[[i]],
                              z = dataToCalcDist,
                              options = options)
        if (is.null(tempResult$outlyingnessZ)) {
          warning(paste("Adjusted outlyingness to training group", i,
                        "could not be computed.")
          )
        } else {
          resultDistance[, i] <- tempResult$outlyingnessZ
        }
      }
    }
    else if (dataType == "functional") {
      tempResult <- fOutl(x = trainingData[[i]], 
                          z = dataToCalcDist,
                          type = type,
                          alpha = options$alpha,
                          time = options$time,
                          distOptions = options
                          )
      if (is.null(tempResult$fOutlyingnessZ)) {
        warning(paste("Functional outlyingness to training group", i,
                      "could not be computed.")
        )
      } else {
        resultDistance[, i] <- tempResult$fOutlyingnessZ
      }
    }
    
  }

  return(resultDistance)

}
