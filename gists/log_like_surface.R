b <- seq(from=0.0000975,to=0.0000993,length=50)
p <- seq(0.36,0.41,length=50)
grid <- expand.grid(b=b,p=p)
grid <- ddply(grid,~b+p,mutate,loglik=f6(b,p))
grid <- subset(grid,is.finite(loglik))
ggplot(grid,aes(x=b,y=p,z=loglik,fill=loglik))+
  geom_tile()+geom_contour(binwidth=1)
