############ model for diff coefs, sharing sigma over buildings
fun.esti <- function(df.test){
  ###### linear model estimations
  ####df.test should have: Year, Month, CNGAccount, CCF, it it tbl_df()
  N <- df.test %>% select(Year) %>% unique() %>% unlist() %>% length()
  K <- df.test %>% select(CNGAccount) %>% unique() %>% unlist() %>% length()
  beta.hat <- df.test %>% group_by(Month, CNGAccount) %>% summarise(mean.ccf=mean(CCF))
  sigma.hat <- df.test %>% 
    left_join(beta.hat, by = c("Month","CNGAccount")) %>%
    mutate(s.dif = (CCF-mean.ccf)^2) %>%
    group_by(Month) %>%
    summarise(v = mean(s.dif)/(N*K)) 
  coef.sigma <- left_join(beta.hat, sigma.hat, by = "Month") %>% 
    rename(coefs = mean.ccf, sigma = v) %>%
    arrange(CNGAccount, Month)
  return(coef.sigma)
}

############# the testing
test.procedure.p <- function(df.future, df1.history,
                             N, v, alpha, B1){
  #####df.future and df1.history are simulated tbl_df() data;
  ###df.future each column is an account CCF; 
  ###df1.history first two column is the Year, Month then account CCF;
  ### the order of accounts of df.future and df1.history are the same as well.
  ###simulate coefs 1K times
  K <- ncol(df.future)
  p.p <- c()
  for(k in 1:B1){
    #### the procedure (e)
    sigma.b <- v/rchisq(n = nrow(df.future), df = v) * df1.history$sigma
    
    coefs.b <- apply(df1.history[, -c(1,2)], 2, FUN = function(x){
      rnorm(length(x), 
            mean = x,
            sd = sqrt(sigma.b/N))
    })
    p.p <- c(p.p, max(apply(rbind(df.future, coefs.b), 2, FUN = function(x){
      ###x is vector with the future values and the estimates
      ## the first half are future values and the second half are estimates
      x.half = length(x)/2
      sum((x[1:x.half] - x[-c(1:x.half)])*x[-c(1:x.half)]/sigma.b)/
        sqrt(sum(x[-c(1:x.half)]^2/sigma.b))
    })) < qnorm((1-alpha)^(1/K)) )
  }
  return(mean(p.p))
}

statistic <- function(df.future, df1.history,
                      N, v, alpha, B1){
  p.p <- c()
  for(k in 1:B1){
    #### the procedure (e)
    sigma.b <- v/rchisq(n = nrow(df.future), df = v) * df1.history$sigma
    
    coefs.b <- apply(df1.history[, -c(1,2)], 2, FUN = function(x){
      rnorm(length(x), 
            mean = x,
            sd = sqrt(sigma.b/N))
    })
    p.p <- cbind(p.p, apply(rbind(df.future, coefs.b), 2, FUN = function(x){
      ###x is vector with the future values and the estimates
      ## the first half are future values and the second half are estimates
      x.half = length(x)/2
      sum((x[1:x.half] - x[-c(1:x.half)])*x[-c(1:x.half)]/sigma.b)/
        sqrt(sum(x[-c(1:x.half)]^2/sigma.b))
    }) )
  }
  return(p.p)
}
########### monitoring new data 
real.monitoring <- function(df.test, df.future, iter, monitor.year){
  alpha = 0.05
  alpha.gamma <- 0.05
  SM = 12
  B1 <- B2 <- 1e2
  
  N <- df.test %>% select(Year) %>% unique() %>% unlist() %>% length()
  K <- df.test %>% select(CNGAccount) %>% unique() %>% unlist() %>% length()
  coef.sigma <- fun.esti(df.test)
  p.alpha <- c()
  ####history y (b)
  df1 = spread(coef.sigma, CNGAccount, coefs)
  # ####generate historical y (b)
  df.history <- apply(df1[,-c(1,2)], 2, FUN = function(x){
    c(replicate(N, rnorm(12, 
                         mean = x,
                         sd = sqrt(df1$sigma))))
  })
  df.history <- df.history %>% tbl_df() %>%
    mutate(Year = rep(1:N, each = 12), 
           Month = rep(1:12, times = N)) 
  df.history <- df.history %>%
    gather("CNGAccount", "CCF", -c(Year,Month))
  
  v <- 12*N - 12
  # ####all the beta, sigma estimates (c)
  coef.sigma1 <- fun.esti(df.history);
  df1.history = spread(coef.sigma1, CNGAccount, coefs)
  
  p.1 <- replicate(B2, test.procedure.p(
    apply(df1.history[,-c(1,2)], 2, FUN = function(x){
      rnorm(SM, 
            mean = x,
            sd = sqrt(df1.history$sigma))
    }), df1.history[1:SM, ],
    N, v, alpha, B1) )
  
  gamma.p <- quantile(p.1, probs = alpha.gamma)
  # pp.1 <- rbind(pp.1, p.1)
  
  ### test for true type I error
  # p.p <- test.procedure.p(df.future, df1.history[1:SM, ],
  #                         N, v, alpha, B1)
  stat.future <- statistic(df.future, df1[1:SM, ],
                           N, v, alpha, B1)
  
  p.p <- mean(apply(stat.future, 2, max) < qnorm((1-alpha)^(1/K)) )
  
  p.alpha <- p.p < gamma.p
  
  p.each <- rowMeans(stat.future< qnorm((1-alpha)^(1/K)))<gamma.p
  
  
  return(c(list(coef.sigma), list(p.alpha), list(p.each)))
}
