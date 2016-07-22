#' Predicting Cox Proportional Hazard model
#' @export
#Calculates Cox model and returns list of prediciton for a user specified time
evalCox <- function(train,test=train,covariates,time,surv=TRUE){
    require(survival)
    if(!all(names(train)==names(test))){
        stop("names of train has to be the same as the names of test.")
    }
    if(!all(c("fu","event") %in% names(train))){
        stop("follow up time in data has to be called fu and the event has to be called event.")
    }
    surv_object <- with(train,Surv(as.numeric(fu),event))
    #Making formula from covariates
    formula <- paste("surv_object",paste(covariates,collapse="+"),sep="~")
    coxModel <- coxph(as.formula(formula),data=train)
    #Use model to predict outcome in testset
    baseSurvFit <- survfit(coxModel)
    baseSurv <- approx(baseSurvFit$time,baseSurv$surv,xout=time)$y
    LinearPredictor <- predict(coxModel,newdata=test)
    if(surv){
        predictionFUN <- function(i) i^exp(LinearPredictor)
    }else{
        predictionFUN <- function(i) 1-i^exp(LinearPredictor)
    }
    if(!all(time %in% baseSurv$surv)){
        stop("time argument has to match follow up times from the data.")
    }
    prediction <- lapply(baseSurv,FUN=predictionFUN)
    names(prediction) <-  time
    return(prediction)
}
