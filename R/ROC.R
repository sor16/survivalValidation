#' ROC curves
#'
#' @param train a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion (k-1)/k of the original data set.
#' @param test a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion 1/k of the original data set.
#' @param FittingFunction a function which returns list of survival probabilities of individuals where each element in the list represent a specific time point.
#' @param covariates character vector specifying the names of covariates
#' @param time a numeric value specifying at what time survival probability is to be calculated.
#' @param by a numeric value specifying the spacing between values in the proportion vector with lower and upper values from xlim
#' @return The function returns a data frame of the proportion that had an event for a given survival probability including confidence intervals
#' @examples
#' ROC(train)
#' @export
#'
ROC  <- function(train,test,FittingFunction,time,covariates){
    require(dplyr)
    if(all(names(train)!=names(test))){
        stop("names of train has to be the same as the names of test.")
    }
    if(all(!c("fu","event") %in% names(train))){
        stop("follow up time in data has to be called fu and the event has to be called event.")
    }
    #Check if patient has had an event at a user specified time
    statusAtTime <- with(test,fu < time & event)
    #Fit a model on train dataset and return prediction on test set
    prediction <- FittingFunction(train=train,test=test,surv=FALSE,time=time,covariates=covariates) %>% unlist()
    threshold <- seq(0,1,length.out=100)
    ROCData <- lapply(threshold,function(i){
        eventPredict = prediction >= i
        TPR=sum(eventPredict & statusAtTime)/sum(statusAtTime)
        FPR=sum(eventPredict & !statusAtTime)/sum(!statusAtTime)
        data.frame(TPR=TPR,FPR=FPR)
    }) %>% bind_rows() %>% arrange(-row_number())
    ROCData$threshold <- threshold
    #trapisunalgun
    AreaOfTrapezoids=with(ROCData,c(0,(FPR[2:(length(TPR))]-FPR[1:(length(TPR)-1)])*0.5*(TPR[1:(length(TPR)-1)]+TPR[2:(length(TPR))])))
    ROCData$AUC=sum(AreaOfTrapezoids)
    ROCData$cumArea=cumsum(AreaOfTrapezoids)
    return(ROCData)
}
