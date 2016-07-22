#' Integrated Brier Score
#'
#' @param train a data.frame which includes follow up time (fu), event, and covariates. In a k-fold cross validation this is a sample proportion (k-1)/k of the original data set.
#' @param test a data.frame which includes follow up time (fu), event, and covariates. In a k-fold cross validation this is a sample proportion 1/k of the original data set.
#' @param FittingFunction a function which returns list of survival probabilities of individuals where each element in the list represent a specific time point.
#' @param covariates character vector specifying the names of covariates
#' @examples
#' IBS(train)
#' @export
#'
IBS <- function(train,test,FittingFunction,covariates,...){
    require(survival)
    if(all(names(train)!=names(test))){
        stop("names of train has to be the same as the names of test.")
    }
    if(all(!c("fu","event") %in% names(train))){
        stop("follow up time in data has to be called fu and the event has to be called event.")
    }
    #Kaplan Meier estimate of censoring survival function on the train set
    censor_object=with(train,Surv(as.numeric(fu),event=!event))
    censorFit <- survfit(censor_object ~ 1)
    timeVector <- censorFit$time
    censorSurvPred <- censorFit$surv

    #Fit a model on train dataset and return prediction on test set
    survPred <- FittingFunction(train=train,test=test,surv=TRUE,time=timeVector,covariates=covariates,...)
    #evaluate Brier score at times that survPred allows for
    BSVector <- sapply(1:length(timeVector),function(i){
        EventPreTime <- with(test,fu < timeVector[i] & event)
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
