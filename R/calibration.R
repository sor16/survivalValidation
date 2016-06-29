#' Cross Validation
#'
#' @param train a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion (k-1)/k of the original data set.
#' @param test a data.frame which includes follow up time, event, and covariates. In a k-fold cross validation this is a sample proportion 1/k of the original data set.
#' @param FittingFunction a function which returns list of survival probabilities of individuals where each element in the list represent a specific time point.
#' @param formula a formula object corresponding to FittingFunction
#' @param time a numeric value specifying at what time survival probability is to be calculated.
#' @param xlim  a two element vector with values between 0 and 1 specifying lower and upper limit of proportions the function should calculate calibration
#' @param by a numeric value specifying the spacing between values in the proportion vector with lower and upper values from xlim
#' @return The function returns a data frame of the proportion that had an event for a given survival probability including confidence intervals
#' @examples
#' @export
#'
calibration <- function(train,test=train,FittingFunction,formula,time,xlim=c(0,1),by=0.05){
    require(survival)
    require(dplyr)
    #Check if patient has had an event at a user specified time
    statusAtTime <- with(test,fu < time & dead)
    #Fit a model on train dataset and return prediction on test set
    prediction <- FittingFunction(train=train,test=test,formula=formula,time=time,surv=FALSE) %>% unlist() %>% as.numeric()
    #Calculate calibration
    deciles <- seq(xlim[1],xlim[2],by=by)
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

    }) %>% rbind_all()
    return(calibrationData)
}
