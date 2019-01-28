# This script is copied from the R MetaIntegrator package "immunoStates.R" script. It is thus under a LGPL license. To review the LGPL license, see License_immunoStates_metaintegrator.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#iSdeconvolution
#######################################################################################
#Run Linear Regression model on a Gene Expression MATRIX
#######################################################################################
iSdeconvolution <- function(basisMatrix,#default to immunoStates
                            geneExpressionMatrix){
  
  #format basis matrix and expression data
  basisMatrix          <- data.matrix(basisMatrix)
  geneExpressionMatrix <- data.matrix(geneExpressionMatrix)
  
  #order by row-names
  basisMatrix <- basisMatrix[order(rownames(basisMatrix)),]
  geneExpressionMatrix <- geneExpressionMatrix[order(rownames(geneExpressionMatrix)),]
  
  #convert into real [non-log] space if in log space [be wary of NAs]
  if(max(geneExpressionMatrix,na.rm=T) < 50){
    geneExpressionMatrix <- 2^geneExpressionMatrix
  }
  
  #run quantile normalization on gene expression matrix
  colN <- colnames(geneExpressionMatrix)
  rowN <- rownames(geneExpressionMatrix)
  geneExpressionMatrix <- preprocessCore::normalize.quantiles(geneExpressionMatrix)
  colnames(geneExpressionMatrix) <- colN
  rownames(geneExpressionMatrix) <- rowN
  
  #only pick genes in common between data and basis matrix
  bMgenes  <- row.names(basisMatrix)
  gMgenes  <- row.names(geneExpressionMatrix)
  GMintBM  <- gMgenes %in% bMgenes
  
  if(sum(GMintBM)==0) {
    stop("None of the gene names are present in the the basis matrix. Make sure your gene names are correct and standard.")
  }
  
  geneExpressionMatrix <- geneExpressionMatrix[GMintBM,]
  BMintGM              <- bMgenes %in% row.names(geneExpressionMatrix)
  basisMatrix          <- basisMatrix[BMintGM,]
  
  #standardize basis matrix [rescale globally]
  basisMatrix <- (basisMatrix-mean(basisMatrix,na.rm=T))/stats::sd(as.vector(basisMatrix),na.rm=T)
  
  #declare header for the output matrix
  header <- c('Sample',colnames(basisMatrix),"P-value","Correlation","RMSE")
  
  #declare empty output matrix
  output   <- matrix()
  pvalue   <- 9999
  
  #run for every sample
  sampleIndex <- 1
  while(sampleIndex <= ncol(geneExpressionMatrix)){
    #iterate one variable @ the time
    geneExpressionSample <- geneExpressionMatrix[,sampleIndex]
    
    #remove NAs [this is necessary to avoid issues downstream]
    basisMatrixSample    <- basisMatrix[which(!is.na(geneExpressionSample)),]
    geneExpressionSample <- geneExpressionSample[which(!is.na(geneExpressionSample))]
    
    #scale cell-mixture sample
    geneExpressionSample <- scale(geneExpressionSample)
    
    #run linear regression on a single sample
    decOut  <- DecLinearRegression(basisMatrixSample,geneExpressionSample)
    
    #create output vector
    out <- c(colnames(geneExpressionMatrix)[sampleIndex],
             decOut$props,
             pvalue,
             decOut$r,
             decOut$rmse)
    
    #update output matrix
    if(sampleIndex == 1){
      output <- out
    }else{
      output <- rbind(output, out)
    }
    
    #update sample inder
    sampleIndex <- sampleIndex + 1
  }
  
  #format matrix object containing all results
  outObj <- rbind(header,output)
  outObj <- outObj[,-1]
  outObj <- outObj[-1,]
  outObj <- matrix(as.numeric(unlist(outObj)),nrow=nrow(outObj))
  
  #object is a matrix of samples X cells [+ some output]
  rownames(outObj) <- colnames(geneExpressionMatrix)
  colnames(outObj) <- c(colnames(basisMatrix),"P-value","Correlation","RMSE")
  
  #return object
  return(outObj)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DecLinearRegression
#######################################################################################
#Run Linear Regression model on a Gene Expression MATRIX [genes not probes]
#######################################################################################
DecLinearRegression <- function(xMat,yVector){
  
  #run linear model without intercept
  model    <- stats::lm(yVector ~ xMat -1)
  
  #go from coefficients to proportions
  coeff    <- model$coefficients
  coeff[coeff < 0] <- 0
  props    <- coeff/sum(coeff)
  
  #get RMSE and Correlation for deconvolution
  dec_r    <- stats::cor(model$fitted.values,yVector)
  dec_rmse <- Metrics::rmse(model$fitted.values,yVector)
  
  #return final list in output
  newList <- list("props" = props, "rmse" = dec_rmse, "r" = dec_r)
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("gene","rn","value","variable","keys","immunoStatesMatrix","natural_killer_cell","CD56bright_natural_killer_cell",
                         "CD56dim_natural_killer_cell","monocyte","CD14_positive_monocyte","CD16_positive_monocyte","B_cell","naive_B_cell",
                         "memory_B_cell","plasma_cell","T_cell","CD8_positive_alpha_beta_T_cell","CD4_positive_alpha_beta_T_cell",
                         "gamma_delta_T_cell","granulocyte","neutrophil","eosinophil","basophil","gse","mdh"))