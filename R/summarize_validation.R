#' Summarize Validation method
#' Computes aggregated value from one of the validation methods, \link[survivalValidation]{IBS},
#' \link[survivalValidation]{calibration},\link[survivalValidation]{ROC},\link[survivalValidation]{c_statistic}.
#' @export
summarize_validation <- function(cvList,validationMethod,simplify=FALSE){
    require(dplyr)
    if(validationMethod=="IBS"){
        summarizedData <- sapply(cvList,unlist) %>% rowMeans() %>% data.frame(IntegratedBrierScore = (.))
        if(simplify) return(summarizedData %>% summarise(IBS=median(IntegratedBrierScore)))
    }else if(validationMethod=="ROC"){
        summarizedData <- lapply(cvList,bind_rows) %>% bind_rows() %>% group_by(threshold) %>%
            summarise(TPR=mean(TPR,na.rm=T),FPR=mean(FPR,na.rm=T),AUC=mean(AUC,na.rm=T),cumArea=mean(cumArea,na.rm=T)) %>%
            ungroup() %>% distinct(TPR,FPR,.keep_all=TRUE)
        if(simplify) return(summarizedData %>% distinct(AUC))
    }else if(validationMethod=="calibration"){
        summarizedData <- lapply(cvList,bind_rows) %>% bind_rows() %>% group_by(deciles) %>%
        summarise(proportion=mean(proportion,na.rm=T),lower=mean(lower,na.rm=T),upper=mean(upper,na.rm=T))
        linearModel=with(summarizedData,lm(proportion~deciles))
        if(simplify) return(data.frame(A = linearModel$coefficients[[1]], B = linearModel$coefficients[[2]], MSE=mean(linearModel$residuals^2)))
    }else if(validationMethod=="c_statistic"){
        summarizedData <- cvList %>% sapply(unlist) %>% rowMeans() %>% data.frame(c_statistic = (.))
        if(simplify) return(summarizedData %>% summarise(c_statistic=median(c_statistic)))
    }
    return(summarizedData)
}
