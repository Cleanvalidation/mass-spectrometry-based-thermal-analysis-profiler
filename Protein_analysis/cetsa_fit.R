#' curve fitting function
#'
#' CETSA curve fitting function
#' Fit a curve to a set of data
#'
#' @param d.  data frame containing temperature and value or norm_value columns
#' @param norm.  Boolean.  If true then use norm_value column otherwise use value column
#'
#' @return nls2 curve model
#'
#' @import nls2
#'
## fit curve using mstherm equation
cetsa_fit <- function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    #c(Pl=0, a = 550, b = 10
    fine_start <- expand.grid(p=c(0,0.1),k=seq(0,100),m=seq(5,45,by=10))
    new_start <- nls2::nls2(y ~ fit.cetsa(p, k, m, t),
                            data = myData,
                            start = fine_start,
                            algorithm = "grid-search",#note: check other ones
                            control = nls.control(warnOnly=T,maxiter=5000))
    nls2::nls2(y ~ fit.cetsa(p, k, m, t),
               data = myData,
               start = new_start,
               control = nls.control(warnOnly=F),
               algorithm = "grid-search",
               lower = c(0, 1, 10),
               upper = c(0.4, 100000, 100))
  }, error = function(err) {
    return(NA)
  })
  return(result)
}


cetsa_fit2 <- function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    #c(Pl=0, a = 550, b = 10
    fine_start <- expand.grid(p=c(0,0.2),k=seq(0,1000),m=seq(5,100,by=5))
    new_start <- nls2::nls2(y ~ fit.cetsa(p, k, m, t),
                            data = myData,
                            start = fine_start,
                            algorithm = "grid-search",#note: check other ones
                            control = nls.control(warnOnly=T,maxiter=5000))
    nls2::nls2(y ~ fit.cetsa(p, k, m, t),
               data = myData,
               start = new_start,
               control = nls.control(warnOnly=F),
               algorithm = "grid-search",
               lower = c(0, 1, 10),
               upper = c(0.4, 100000, 100))
  }, error = function(err) {
    return(NA)
  })
  fgh2 <- deriv(value ~ (1 - p)/(1 + exp(-k*(1/t - 1/m))) + p, c("p","k","m"), function(p,k,m){} ) 
  #create a grid
  t<-seq.int(37,67,0.05)
  beta2.est <- coef(new_start)
  f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3])
  
  g.new <- base::attr(f.new,"gradient")
  g.new1 <- t(as.vector(attr(f.new,"gradient")))
  V.beta2 <- vcov(new_start)
  GS<-rowSums((g.new%*%V.beta2)*g.new)
  #95% CI
  alpha <- 0.05
  deltaf <- sqrt(GS)*qt(1-alpha/2,summary(new_start)$df[2])
  df.delta <- data.frame(temperature=t, value=f.new, lwr.conf=f.new-deltaf, upr.conf=f.new+deltaf)
  sigma2.est <- summary(new_start)$sigma
  deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,summary(new_start)$df[2])
  df.delta[c("lwr.pred","upr.pred")] <- cbind(f.new - deltay,f.new + deltay)
  #pl<-ggplot(d)+geom_point(mapping=aes(x=temperature,y=value))
  # 
  # pl + geom_ribbon(data=df.delta, aes(x=t, ymin=lwr.pred, ymax=upr.pred), alpha=0.1, fill="blue") +
  #   geom_ribbon(data=df.delta, aes(x=t, ymin=lwr.conf, ymax=upr.conf), alpha=0.2, fill="#339900") +
  #   geom_line(data=df.delta, aes(x=t, y=f.new), colour="#339900", size=1)
  return(list(result,df.delta))
}
