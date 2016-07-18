#' Validation plot
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
ValidationPlot <- function(cvList,validationMethod,evaluationFunctions=NULL,legendPosition=NULL,animate=FALSE){
    require(ggplot2)
    require(dplyr)
    if(validationMethod=="IBS"){
        if(is.null(evaluationFunctions)) stop("Please insert character vector for which models were used.")
        IBSData <- sapply(cvList,unlist) %>% rowMeans() %>% data.frame(IntegratedBrierScore = (.))
        ModelType <- factor(evaluationFunctions)
        validationPlot <- ggplot(data = IBSData,aes(ModelType,IntegratedBrierScore)) + geom_boxplot()
    }else if(validationMethod=="ROC"){
        FullROCData <-  lapply(cvList,bind_rows) %>% bind_rows()
        summarizedData <- FullROCData %>% group_by(threshold) %>%
            summarise(TPR=mean(TPR,na.rm=T),FPR=mean(FPR,na.rm=T),AUC=mean(AUC,na.rm=T),cumArea=mean(cumArea,na.rm=T)) %>%
            ungroup() %>% distinct(TPR,FPR,.keep_all=TRUE)
        AUCposition=data.frame(AUC=unique(FullROCData$AUC),x=0.75,y=0.25)
        validationPlot <- ggplot(data=summarizedData,aes(FPR,TPR,label=cumArea)) + geom_line() + geom_abline(slope=1,intercept=0) +
            geom_text(data=AUCposition,aes(x,y,label=paste("AUC=",round(AUC,3)))) + xlim(0,1)+ylim(0,1) + theme_bw() +
            ylab("True positive rate") + xlab("False positive rate")
    }else{
        fullCalibrationData <-  lapply(cvList,bind_rows) %>% bind_rows()
        summarizedData <- fullCalibrationData %>% group_by(deciles) %>% summarise(proportion=mean(proportion,na.rm=T),lower=mean(lower,na.rm=T),upper=mean(upper,na.rm=T))
        linearModel=with(summarizedData,lm(proportion~deciles))
        xlim <- with(summarizedData,range(deciles[!is.na(proportion)])) + c(-0.05,0.05)
        ylim <- with(summarizedData,c(min(lower,na.rm=T),max(upper,na.rm=T)))
        if(is.null(legendPosition)){
            legendPosition <- data.frame(x=rep(xlim[1]+0.1,2),y=c(3*ylim[2]/4.,3*ylim[2]/4 - 0.1))
        }
        textDat= data.frame(legendPosition,text=paste(c("A","B"),format(linearModel$coefficients,digits = 2),sep=" = "))
        validationPlot <- ggplot(data=summarizedData,aes(deciles,proportion)) + geom_point(na.rm=T) +
            geom_smooth(method = "lm",se=F,na.rm=T) + geom_text(data=textDat,aes(x,y,label=text)) +
            geom_errorbar(aes(ymin = lower,ymax = upper,width=0.01,color="red"),na.rm=T) + theme_bw() +
            xlim(xlim) + ylim(ylim) +theme(legend.position="none")
    }
    if(animate){
        require(plotly)
        validationPlot <- ggplotly(validationPlot)
    }
    return(validationPlot)
}





