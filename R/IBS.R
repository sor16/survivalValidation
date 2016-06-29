#' Integrated Brier Score
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
IBS <- function(train,test,FittingFunction,formula){
    require(survival)
    #Kaplan Meier estimate of censoring survival function on the train set
    censor_object=with(train,Surv(as.numeric(fu),event=!dead))
    censorFit <- survfit(censor_object ~ 1)
    censorSurvPred <- censorFit$surv
    timeVector <- censorFit$time
    #Fit a model on train dataset and return prediction on test set
    survPred <- FittingFunction(train=train,test=test,surv=TRUE,time=timeVector,formula=formula)

    BSVector <- sapply(1:length(timeVector),function(i){
        EventPreTime <- with(test,fu < timeVector[i] & dead)
        AliveAtTime <- with(test,fu > timeVector[i])
        BS=mean((survPred[[i]]^2*EventPreTime/censorSurvPred[i]) + ((1-survPred[[i]]^2)*AliveAtTime/censorSurvPred[i]))
        BS
    })
    #Estimating integral of the Brier score function from 0 to t_max
    h=diff(timeVector)
    areaOfTrapezoids=h*0.5*(BSVector[1:(length(BSVector)-1)]+BSVector[2:(length(BSVector))])
    IntegratedBS <- sum(areaOfTrapezoids,na.rm=T)/max(timeVector)
    return(data.frame(IntegratedBS=IntegratedBS))
}
