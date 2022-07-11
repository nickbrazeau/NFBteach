niameyA2 <- pomp(
  data=subset(niamey,community=="A",select=-community),
  times="biweek",t0=0,
  skeleton=vectorfield(
    Csnippet("
      double incidence;
      incidence = b*S*I;
      DS = -incidence;
      DI = incidence-gamma*I;")),
  initializer=Csnippet("
      S = S_0;
      I = I_0;"),
  paramnames=c("b","gamma","S_0","I_0"),
  statenames=c("S","I"))



poisson.loglik <- function (params) {
  x <- trajectory(niameyA2,params=params)
  sum(dpois(x=obs(niameyA2),lambda=params["p"]*x["I",,],log=TRUE))
}


logit <- function (p) log(p/(1-p))    # the logit transform
expit <- function (x) 1/(1+exp(-x))   # inverse logit

f5 <- function (par) {
  params <- c(S_0=20000,I_0=1,gamma=1,b=exp(par[1]),p=expit(par[2]))
  -poisson.loglik(params)
}

fit5 <- optim(f5,par=c(log(0.0001),logit(0.2)))
fit5
