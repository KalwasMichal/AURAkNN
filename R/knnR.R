

kNN_impute <- function(data,k=5,metric=c("euclidean","manhattan"),
                       threads=NULL,maxColNa=0.8,maxRowNa=0.5)
{
  chosenMetric <- match.arg(metric)
  metricFlag <- switch (chosenMetric,
                        "euclidean" = 1L,
                        "manhattan" = 2L
  )
  
  if(!is.data.frame(data) && !is.matrix(data))
  {
    stop("Data must be data frame or matrix")
  }
  if(k<=0){stop("k must be positive")}
  
  if(is.null(threads) || is.na(threads))
  {
    threads <- max(1L,as.integer(parallel::detectCores())-1L)
  }
  else
  {
    threads <- as.integer(threads)
  }
  if(maxColNa <=0 || maxColNa >=1 || maxRowNa <=0 || maxRowNa >=1 )
  {
    stop("skibidi")
  }
  
  dfOutput <- as.data.frame(data)
  
  # Sprawdzamy liczbę Na w kolumnach i odsiewamy je
  DeletedCols <-names(dfOutput) [colSums(is.na(dfOutput)) / nrow(dfOutput) > maxColNa]
  
  dfOutput <- dfOutput[,!names(dfOutput) %in% DeletedCols,drop=FALSE]
  
  # Sprawdzamy liczbę Na w wierszach i odsiewamy je
  
  keepRows <- rowSums(is.na(dfOutput))/ncol(dfOutput) <= maxRowNa
  
  dfOutput <- dfOutput[keepRows, ,drop=FALSE]
  
  #Sprawdzamy pozostałą populacje
  
  if(nrow(dfOutput)<k){
    stop("Chosen k is too high")
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
  
  
  return(dfFinal)
  
  
  
  
}