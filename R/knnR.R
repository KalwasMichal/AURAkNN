

computeWeights <- function(dfOutput){
  colN <-ncol(dfOutput)
  rowN <- nrow(dfOutput)
  
  
  corrMatrix <- matrix(data=rep(c(0),length.out=colN*colN),nrow=colN,ncol=colN)
  
  for(i in 1:colN)
  {
    for(j in 1:colN)
    {
      if(i == j)
      {
        next
      }
      
      iCol <- dfOutput[,i]
      jCol <- dfOutput[,j]
      validIdx <- complete.cases(iCol,jCol)
      iCol <- iCol[validIdx]
      jCol<- jCol[validIdx]
      
      
      if(length(iCol) <3 || length(unique(iCol)) <2 || length(unique(jCol)) <2 || length(jCol) <3) {
        corrMatrix[i, j] <- 0.0
        next
      }
      
      iColType <- is.numeric(iCol)
      jColType <- is.numeric(jCol)
      
      if(iColType && jColType)
      {
        corrMatrix[i,j]=abs( cor(iCol,jCol,method="spearman"))
      }
      else if(!iColType && !jColType)
      {
        t <- table(iCol,jCol)
        r <- nrow(t)
        k <-ncol(t)
        n <- sum(t)
        if(r<2 || k<2)
        {
          corrMatrix[i,j]=0.0
          next
        }
        chi <- suppressWarnings(chisq.test(t,correct=FALSE)$statistic)
        corrMatrix[i,j] <- as.numeric(sqrt(chi/(n*min(r-1,k-1))))
      }
      else
      {
        if(iColType)
        {
          numVar <-iCol
          catVar <- jCol
        }
        else
        {
          numVar <- jCol
          catVar <- iCol
        }
        corrMatrix[i,j] <- as.numeric(summary(lm(numVar~catVar))$r.squared)
      }
    }
  }
  for(i in 1:colN)
  {
    maxVal <- max(corrMatrix[i,])
    if(maxVal>0)
    {
      corrMatrix[i,]<- corrMatrix[i,]/maxVal
    }
    
  }
  return(corrMatrix)
}




kNN_impute <- function(data,k=5,metric=c("gower","euclidean","manhattan"),
                       mode=c("fast","precise"),num_fun=c("median","mean"),threads=NULL,maxColNa=0.8,maxRowNa=1.0)
{
  chosenMetric <- match.arg(metric)
  metricFlag <- switch (chosenMetric,
                        "euclidean" = 1L,
                        "manhattan" = 2L,
                        "gower"=3L
  )
  chosenMode <- match.arg(mode)
  modeFlag <- switch (chosenMode,
                      "fast" = 0L,
                      "precise" = 1L
  )
  chosenFun <- match.arg(num_fun)
  funFlag <- switch (chosenFun,
                     "median" = 1L,
                     "mean" = 2L
  )
  
  if(!is.data.frame(data) && !is.matrix(data))
  {
    stop("Data must be data frame or matrix")
  }
  if(k<=0){stop("k must be positive")}
  
  if(is.null(threads) || is.na(threads))
  {
    cores <- parallel::detectCores()
    if (is.na(cores)) cores <- 1L
    threads <- max(1L, as.integer(cores) - 1L)
  }
  else
  {
    threads <- as.integer(threads)
  }
  if(maxColNa <0 || maxColNa >1 || maxRowNa <0 || maxRowNa >1 )
  {
    stop("Invalid maxColNa or maxRowNa, both must be <0,1>")
  }
  
  dfOutput <- as.data.frame(data)
  
  # Sprawdzamy liczbę Na w kolumnach i odsiewamy je
  DeletedCols <-names(dfOutput) [colSums(is.na(dfOutput)) / nrow(dfOutput) > maxColNa]
  
  dfOutput <- dfOutput[,!names(dfOutput) %in% DeletedCols,drop=FALSE]
  
  # Sprawdzamy liczbę Na w wierszach i odsiewamy je
  
  keepRows <- rowSums(is.na(dfOutput))/ncol(dfOutput) <= maxRowNa
  
  dfOutput <- dfOutput[keepRows, ,drop=FALSE]
  
  
  shuffle_idx <- sample(nrow(dfOutput))
  dfOutput <- dfOutput[shuffle_idx, , drop = FALSE]
  
  #Sprawdzamy pozostałą populacje
  
  if(nrow(dfOutput)<k){
    stop("Chosen k is too high")
  }
  
  weights <- NULL
  if(modeFlag== 1L)
  {
    weights <- as.numeric(computeWeights(dfOutput))
  }
  
  # One hot encoding i skalowanie dla metryki euklidesowej i manhattan
  
  if(metricFlag %in% c(1L,2L))
  {
    
    NumericCols <- sapply(dfOutput,is.numeric)
    
    # Pomijamy kolumny typu factor etc, które mają tylko jedną możliwą wartość
    numberOfuniqueVal <- sapply(dfOutput[!NumericCols],function(x){length(unique(na.omit(x)))})
    
    factorLikeColToDelete <- names(numberOfuniqueVal[numberOfuniqueVal==1])
    if(length(factorLikeColToDelete) >0)
    {
      dfOutput <- dfOutput[,!names(dfOutput) %in% factorLikeColToDelete]
      NumericCols <- sapply(dfOutput,is.numeric)
    }
    
    
    #Robimy One hot encoding
    options(na.action = "na.pass")
    
    matrixOutput <- model.matrix(~. -1,data=dfOutput)   
    
    options(na.action = "na.omit")
    
    # Skalujemy wartości numeryczne
    numericColsNames <- names(dfOutput[NumericCols])
    
    scale_mean <- numeric(ncol(matrixOutput))
    scale_sd <- rep(1.0,ncol(matrixOutput))
    names(scale_mean) <- colnames(matrixOutput)
    names(scale_sd) <-colnames(matrixOutput)
    
    for(col in numericColsNames)
    {
      colMean <- mean(matrixOutput[,col],na.rm=TRUE)
      colSd <- sd(matrixOutput[,col],na.rm=TRUE)
      if(is.na(colSd) || colSd==0) {colSd <-1.0}
      
      scale_mean[col] <- colMean
      scale_sd[col] <- colSd
      
      
      matrixOutput[,col] <- (matrixOutput[,col]-colMean)/colSd
    }
    
    # Wywołujemy funkcję z C
    imputedMatrixC <- .Call("knn_imputeC",as.matrix(matrixOutput),
                            as.integer(k),as.integer(metricFlag),as.integer(threads))
    
    
    
    # Dekodujemy macierz
    dfFinal <- dfOutput
    
    for(col in numericColsNames)
    {
      naIdx <- is.na(dfFinal[[col]])
      
      if(any(naIdx))
      {
        decodedVal <- imputedMatrixC[,col]*scale_sd[col]+scale_mean[col]
        
        dfFinal[[col]][naIdx] <- decodedVal[naIdx]
      }
    }
    #Odkręcamy One Hot Encoding
    
    factorLikeColsNames <- names(dfOutput[!NumericCols])
    
    for(col in factorLikeColsNames)
    {
      naIdx <- is.na(dfFinal[[col]])
      if(any(naIdx))
      {
        if(!is.factor(dfFinal[[col]]))
        {
          dfFinal[[col]] <- as.factor(dfFinal[[col]])
        }
        
        colLevels <- levels(dfFinal[[col]])
        oheCols <- paste0(col,colLevels)
        
        oheCols <- intersect(oheCols,colnames(imputedMatrixC))
        
        if(length(oheCols) > 0)
        {
          subMatrix <- imputedMatrixC[naIdx,oheCols,drop=FALSE]
          chosenVal <- max.col(subMatrix,ties.method = "random")
          
          dfFinal[[col]][naIdx] <- colLevels[chosenVal]
        }
        
      }
    }
    
  }
  if(metricFlag==3L)
  {
    colRange <- sapply(dfOutput,function(x){
      if(is.numeric(x)){
        range<- max(x,na.rm = TRUE)-min(x,na.rm=TRUE)
        if(range==0 || is.na(range)) {1.0} else {range}
      }
      else{1.0}
    })
    colType <- sapply(dfOutput, function(x){as.integer(!is.numeric(x))})
    catColNames <- names(dfOutput[as.logical(colType)])
    NumericColNames <- names(dfOutput[!as.logical(colType)])
    
    
    levelsList <- list()
    for(col in catColNames)
    {
      dfOutput[[col]] <- as.factor(dfOutput[[col]])
      levelsList[[col]] <- levels(dfOutput[[col]])
      dfOutput[[col]]<- as.numeric(dfOutput[[col]])
    }

    imputedMatrixC <- .Call("knn_imputeC",as.matrix(dfOutput),
                            as.integer(k),as.integer(metricFlag),as.integer(modeFlag),
                            as.integer(funFlag),as.integer(threads),colType,colRange,weights)
    
    
    
    #fallback
    for(col in NumericColNames)
    {
      naIdx <- is.na(imputedMatrixC[,col])
      if(any(naIdx))
      {
        imputedMatrixC[naIdx,col] <- median(imputedMatrixC[,col],na.rm=TRUE)
      }
    }
    for(col in catColNames)
    {
      naIdx <- is.na(imputedMatrixC[,col])
      if(any(naIdx))
      {
        imputedMatrixC[naIdx,col] <- as.numeric(names(which.max(table(imputedMatrixC[,col]))))
      }
    }
    
    
    
    dfFinal <- dfOutput
    
    for(col in catColNames)
    {
      
      imputedIdx <- round(imputedMatrixC[,col])
      imputedIdx <- pmax(1,pmin(imputedIdx,length(levelsList[[col]])))
      dfFinal[[col]] <- factor(levelsList[[col]][imputedIdx], levels = levelsList[[col]])
    }
    for(col in NumericColNames)
    {
      naIdx <- is.na(dfFinal[[col]])
      
      if(any(naIdx))
      {
        dfFinal[[col]][naIdx] <- imputedMatrixC[naIdx,col]
      }
    }
  }
  dfFinal[shuffle_idx, ] <- dfFinal
  return(dfFinal)
}





