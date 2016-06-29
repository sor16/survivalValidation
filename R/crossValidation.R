#' Cross Validation
#'
#' @param data a data.frame which includes follow up time, event, and covariates.
#' @param k a numeric value specifying k-fold cross validation.
#' @param nrRuns a numeric value indicating the number of runs of cross validation is to be performed.
#' @param ValidationFunciton a function which validates a given model, inbuilt funcitons are IBS,ROC and Calibration.
#' @param ... additional arguments to ValidationFunciton, see \link[survivalValidation]{IBS},
#' \link[survivalValidation]{calibration},\link[survivalValidation]{ROC}. Could also be a userspecified ValidationFunction
#' @return The function returns a large list of length nrRuns. Each element in the list contains all k elements from one run of cross validation.
#' @examples
#' crossValidation(data)
#' @export

crossValidation <- function(data,k,ValidationFunction,nrRuns,...){
    require(dplyr)
    crossValidationList <- list()
    RunsOfCrossValidation <- lapply(1:nrRuns,function(j){
        #k-fold crossvalidation
        indices <- sample(1:k,size=nrow(data),replace=T)
        epsilon <- 10^-1
        lapply(1:k,function(i){
            #split data to test and train
            test <- data %>% filter(abs(indices-i)<epsilon)
            train <- data %>% filter(abs(indices-i)>epsilon)
            ValidationFunction(train=train,test=test,...)
        })
    })
    return(RunsOfCrossValidation)
}
