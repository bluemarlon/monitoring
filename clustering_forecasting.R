PREDICTION.raw <- function(Object,newdata,months, months.pred){
  for(k in 1:6){
    pred <-
      newdata$pred.value[newdata$months==(months+k)] <-
      predict(Object, newdata=newdata[newdata$months==(months+k),])
    newdata$CCF1[newdata$months==(months+k+1)] <- newdata$pred.value[
      newdata$months==(months+k) & newdata$id %in% unique(newdata$id)[which(table(newdata$id)>k)]]
    newdata$CCF6[newdata$months==(months+k+6)] <- pred[
      newdata$Meter[newdata$months==(months+k)] %in% 
        unique(newdata$Meter[newdata$months==(months+k+6)])]# )>=(12-6+k)]
  }
  for(k in 7L:(months.pred-1)){
    pred <- 
      newdata$pred.value[newdata$months==(months+k)] <-
      predict(Object, newdata=newdata[newdata$months==(months+k),])
    newdata$CCF1[newdata$months==(months+k+1)] <- pred[
      newdata$Meter[newdata$months==(months+k)] %in% 
        unique(newdata$Meter[newdata$months==(months+k+1)])]
  }
  newdata$pred.value[newdata$months==(months+months.pred)] <-
    predict(Object, newdata=newdata[newdata$months==(months+months.pred),])
  pred.value <- newdata$pred.value
  pred.value[pred.value<0] <- 0
  pred.value
}

karcher.clustering <- function(df){
  ts.account <- list()
  acc.amp_phase <- acc.fmean <- c()
  n.acc <- 1
  #############################1.cleaning&representation#############################
  #############"fill out the series to be copmlete years"############
  ############Then make the karchar mean####################
  ###with added months
  for(k in unique(df$CNGAccount)){
    # k <- unique(df$Account)[n.acc]
    ts.account[[n.acc]] <- df$CCF[df$CNGAccount==k]#normalized
    ###fill up to full year
    acc.period <- sum(df$CNGAccount==k) %/% 12
    #####series should be more than one year
    if(acc.period>2){
      acc.mov <- (sum(df$CNGAccount==k)+1) %% 12
      if(acc.mov!=0){
        ts.account[[n.acc]] <- 
          c(mean(ts.account[[n.acc]][1:(acc.period)*12]),
            ts.account[[n.acc]],
            colMeans(t(matrix(ts.account[[n.acc]][rep((acc.mov):11, acc.period) + 
                                                    rep(0:(acc.period-1)*12, each=11-acc.mov+1)], 
                              nrow = 11-acc.mov+1))) )
      }else{
        ts.account[[n.acc]] <- 
          c(mean(ts.account[[n.acc]][1:(acc.period)*12]),
            ts.account[[n.acc]] )
      }
      ###karcher mean to every normalized into one year
      karcher.acc <- matrix(ts.account[[n.acc]], nrow = 12)
      acc.out <- time_warping(karcher.acc, time = 1:12, showplot= FALSE) 
      acc.fmean <- cbind(acc.fmean, acc.out$fmean)
      acc.amp_phase <- rbind(acc.amp_phase, c(acc.out$amp.var, acc.out$phase.var))
      n.acc <- n.acc + 1}else{
        # acc.less.two.years <- c(acc.less.two.years, k)
        df <- df[-which(k==df$CNGAccount), ]
      }
  }
  fmean.CID <- diss(t(acc.fmean), METHOD = "CID")
  fmean.hc <- hclust(fmean.CID, method = "ward.D2")
  fmean.memb <- cutree(fmean.hc, k = 6)
  return(fmean.memb)
}

acc.forecasting.change <- function(fmean.memb, df){
  ####################################3.ARIMA################################
  #########################a.monthly pattern#########################
  #########################b.lags1,6,12#########################
  #########################c.temp, snow, CDD#########################
  #########################d.iid lags random effect#########################
  #########################e.predictions&PI on each clust#########################
  cluster.use <- fmean.memb
  mae = c()
  df <- cbind(df, dummy(df$Month))
  colnames(df)[(dim(df)[2]-(length(unique(df$Month))-2) ):dim(df)[2]] <- paste0("L", 2:12)
  
  year <- min(df$Year)
  year.pred <- max(df$Year) #+1 #if keep all
  
  df$CCF1 <- NA
  df$CCF6 <- NA
  df$CCF12 <- NA
  for(meter in unique(df$Meter)){
    df$CCF1[df$Meter==meter] = lag(df$CCF[df$Meter==meter], 1)
    df$CCF6[df$Meter==meter] = lag(df$CCF[df$Meter==meter], 6)
    df$CCF12[df$Meter==meter] = lag(df$CCF[df$Meter==meter], 12)
  }
  df <- df[complete.cases(df[, c("CCF1", "CCF6", "CCF12")]),]
  # df1$CCF1 <- lag(df$CCF, k=1)
  # df1$CCF6 <- lag(df$CCF, k=6)
  # df1$CCF12 <- lag(df$CCF, k=12)
  
  df.pred <- df[df$Year==year.pred, ]
  df.pred$months <- (df.pred$Year-year)*12+df.pred$Month-1
  df <- df[df$Year>=year & df$Year < year.pred, ]
  df$months <- (df$Year-year)*12+df$Month-1
  df <- df[df$months>0,]
  n.acc <- length(unique(df$Meter))
  residual.lr <- res <- c()
  for(clu in unique(cluster.use)){
    idd <- 0L
    df1 <- df[df$Meter %in% unique(df$Meter)[which(cluster.use==clu)], ]
    
    
    df1.pred <- df.pred[df.pred$Meter %in% unique(df$Meter)[which(cluster.use==clu)], ]
    
    for(meter in unique(df1$Meter)) {
      df1.pred$id[df1.pred$Meter==meter] <- df1$id[df1$Meter==meter] <- idd <- idd+1
    }
    remove(idd)
    lr1<-lmer(CCF~(L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+
                     +0)*(CCF1+CCF6+CCF12)+avgt+snow+cdd +(1|id)+(0+CCF1|id)+(0+CCF6|id)+(0+CCF12|id), data = df1)
    residual.lr <- rbind(residual.lr, cbind(residuals(lr1), df1$Meter, df1$months))
    df$residual[df$Meter %in% df1$Meter & df$months %in% df1$months] <- residuals(lr1)
    df$pred[df$Meter %in% df1$Meter & df$months %in% df1$months] <- predict(lr1)
    ####preditcion, PI on each clusters
    df1.pred$CCF1[df1.pred$Month>1] <- NA
    df1.pred$CCF6[df1.pred$Month> 6] <- NA
    pred.value <- PREDICTION.raw(lr1, newdata = df1.pred, months = max(df1$months), months.pred = max(table(df1.pred$Meter)) )
    res <- rbind(res, 
                 cbind(pred.value, residua <- pred.value-df1.pred$CCF, df1.pred$Meter, df1.pred$months))
    df.pred$residual[df.pred$Meter %in% df1.pred$Meter & df.pred$months %in% df1.pred$months] <- residua
    df.pred$pred[df.pred$Meter %in% df1.pred$Meter & df.pred$months %in% df1.pred$months] <- pred.value
  }
  ########################################4.change point#########################
  #########################residuals, original#########################
  ###########try to do the change points on R-ecp
  # K > floor(nrow(Z)/(delta + 1)
  change.series <- change.series.org <-  list()
  for(k in unique(df$Meter)){
    location.id <- df$Meter==k
    delta <- floor(sum(location.id)/3 )
    change.id <- e.cp3o_delta(t(t(df$residual[location.id])), 
                              K=floor(sum(location.id)/(delta+1) ), delta=delta, alpha=1, verbose=FALSE)
    change.id.org <- e.cp3o_delta(t(t(df$CCF[location.id])), K=floor(sum(location.id)/(delta+1) ), 
                                  delta=delta, alpha=1, verbose=FALSE)
    change.series <- c(change.series, list(change.id$estimates))
    change.series.org <- c(change.series.org, list(change.id.org$estimates))
  }
  ####give some plots show
  data.whole <- rbind(df, df.pred)
  data.whole<- data.whole[order(data.whole$UCONNSpace,
                                data.whole$Meter,data.whole$Year, data.whole$Month),]
  data.use <- data.whole
  remove(data.whole)
  
  data.use <- data.use %>% tbl_df() %>% mutate(clu = 1)
  Meters <- unique(df$Meter)
  for(k in cluster.use){
    ss = which(cluster.use==k)
    data.use$clu[ data.use$Meter %in% Meters[ss] ] = k
  }
  change.use <- c() 
  change.use <- change.use %>% tbl_df()
  for(m in 1:length(Meters)){
    change.use <- rbind(change.use, c(Meters[m], 
            data.use$READ.TO[data.use$Meter==Meters[m]][change.series.org[[m]]], 
        data.use$READ.TO[data.use$Meter==Meters[m]][change.series[[m]]]))
  }
  colnames(change.use) <- c("Meter", "change", "res.change")
  change.use <- change.use %>% mutate(Meter.clu = paste0(Meter, ":", fmean.memb),
                                      change = as_date(change),
                                      res.change = as_date(res.change))
  data.use <- data.use %>%select(Meter, CCF, pred, Year, Month, READ.TO, CNGAccount,
                                 clu, AccountType, Facility_Name) %>%
    mutate(Meter.clu = paste0(Meter, ":", clu)) %>%
    gather(key = "Type", value = "val", c("CCF", "pred"))
  
  return(c(list(data.use), list(change.use)))
}

