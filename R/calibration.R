#' Cross Validation
#'
#' @param train a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion (k-1)/k of the original data set.
#' @param test a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion 1/k of the original data set.
#' @param FittingFunction a function which returns list of survival probabilities of individuals where each element in the list represent a specific time point.
#' @param covariates character vector specifying the names of covariates
#' @param time a numeric value specifying at what time survival probability is to be calculated.
#' @param by a numeric value specifying the spacing between values in the proportion vector ranging from 0 to 1.
#' @return The function returns a data frame of the proportion that had an event for a given survival probability including confidence intervals
#' @examples
#' calibration(train,test)
#' @export
#'
calibration <- function(train,test=train,FittingFunction,covariates,time,by=0.1,...){
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
    prediction <- FittingFunction(train=train,test=test,covariates=covariates,time=time,surv=FALSE,...) %>% unlist() %>% as.numeric()
    #Make sequence of proportions
    deciles <- seq(0,1,by=by)
    epsilon <- by/2
    calibrationData <- lapply(deciles,function(i){
        relevantEntries <- statusAtTime[abs(prediction-i)<=epsilon]
        n <- length(relevantEntries)
        if(n<30)
            return(data.frame(deciles=i,proportion=NA,lower=NA,upper=NA,n=n))
        actualDeath <- sum(relevantEntries,na.rm=T)
        estimatedDeath <- round(i*n)
        proportion <- actualDeath/n
        diff <- proportion-i
        p <- prop.test(x=c(actualDeath,estimatedDeath),n=rep(n,2))
        data.frame(deciles=i,proportion=proportion,lower=proportion-abs((diff-p$conf.int[1])),upper=proportion+abs((diff-p$conf.int[2])),n=n)

    }) %>% bind_rows()
    return(calibrationData)
}
