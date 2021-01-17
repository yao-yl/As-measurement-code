# model for the 271 resampled wells
######################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("/Users/YulingYao/Documents/Research/arsenic")
data_repeated308=read.csv("./data/4. repeated308.csv")
data_repeated14_15=read.csv("./data/2014vs15.csv",sep=',') 

head(data_repeated308)
head(data_repeated14_15)

y1= data_repeated308$Arsenic_lab2000
y2= data_repeated308$Arsenic_lab2014
y3=data_repeated14_15$As2015
depth=data_repeated14_15$Depth.2000
loc_x=data_repeated14_15$X2000
loc_y=data_repeated14_15$Y2000
scale_map=max(loc_x)-min(loc_x)
loc_x=(loc_x-min(loc_x))/scale_map
loc_y=(loc_y-min(loc_y))/scale_map

loc_x=loc_x[complete.cases(y3)]
loc_y=loc_y[complete.cases(y3)]
y1=log(y1[complete.cases(y3)])
y2=log(y2[complete.cases(y3)])
y3=log(y3[complete.cases(y3)])


# hyper-parameter is from optimization *(Type-2 MAP)
stan_fit_mu=stan(file="gp_map_mu.stan", iter=2000,  
								 data = list(N=length(y1),y1=y1, y2=y2, y3=y3,
								 						loc_X=cbind(loc_x,loc_y),  
								 						rho= 0.04169228,
								 						alpha=1.06221590, sigma=1.16859570))

# hyper-parameter is from optimization *(Type-2 MAP)
stan_fit=stan(file="gp_resampled.stan", iter=2000,  init=ttt,
									data = list(N=length(y1),y1=y1, y2=y2, y3=y3, 
															loc_X=cbind(loc_x,loc_y), spline_degree=3,
															rho= 0.04169228,
															alpha=1.06221590, sigma=1.16, 
															num_knots=10, knots= unname( quantile(c(y1,y2,y3) ,probs=seq(from=0, to=1, length.out = 10))) ,  
															N_grid=24, 	f_grid= sort( c( seq(-0.8, 6.6, length.out = 20) ,  log (c(40,50,60,80)))) ))
 

 
sss=extract(stan_fit, pars=c("mu", "mu0"))


sss2=extract(stan_fit_mu, pars=c("mu", "mu0"))
f_grid= (sort( c( seq(-0.8, 6.6, length.out = 20) ,  log (c(40,50,60,80))))) 
f12= sss$change_grid
ppb50=c()
ppb100=c()
ppb10=c()
for( i in 1:length(f_grid)){
	ppb10[i]= mean(f_grid[i]+ f12[,i] +   rnorm(dim(f12)[1],  0, 1.16859570) >  log(10))
ppb50[i]= mean(f_grid[i]+ f12[,i] +  rnorm(dim(f12)[1],  0, 1.16859570) >  log(50))
ppb100[i]= mean(f_grid[i]+ f12[,i] +  rnorm(dim(f12)[1],  0, 1.16859570) >  log(100))
}
for( i in 1:24){
	f12[,i]=f12[,i]+sss$sigma_change*rnorm(4000,0,1)
}
library(RColorBrewer)
red_c=brewer.pal(n = 4, name = "YlOrRd")[2:4]
# graphing 
pdf("~/Desktop/spline.pdf", height=2.6, width=7.5)
layout(matrix(c(1:3),nrow=1), width = c(1.2,1,0.6),height = c(1))
f_grid= (sort( c( seq(-0.8, 6.6, length.out = 20) ,  log (c(40,50,60,80))))) 

par(oma=c(1.5,2.7,2,0), pty='m',mar=c(1,1.5,1,1) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=1, cex.lab=0.6, cex.main=0.7) 
plot(f_grid,  colMeans(sss$change_grid),  axes=F,ylab="",xlab="",type='l',col="#8F2727",lwd=1.5, ylim=c(-1.39,1.39 ) , xlim=  log(c(0.9,1000)) , xaxs='i' )
axis(2, las=2,lwd=0.5,at=log(c(1/4, 1/2,  1,2, 4)),  labels  = c("1/4" , "1/2"  ,"1\n no change"  ,2,  4  ) )
abline(h=0,col= 'darkgreen',  lwd=0.8 )
zz=apply(sss$change_grid, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
polygon(x=c(f_grid,rev(f_grid)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=F)
polygon(x=c(f_grid,rev(f_grid)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=F)
axis(1, padj=-0.1,lwd=0.5,at=log(c(0.1,1,10,50,100,1000)),  labels  = c(expression(10^-1), 1  , 10, 50,expression(10^2),expression(10^3) ) )
lines(f_grid,  colMeans(sss$change_grid),col="#8F2727",lwd=1.5, xpd=F)
box(lwd=0.5,bty='l')
abline(h=log(c(1/4,1/2,2 ,4)), lty=2, col='grey40')
abline(v= log( c(10,50, 100)), lty=2, col='grey40')

mtext(1, line=1.3,text="baseline As in 2000 (ppb)",cex=0.8)
mtext(3, text="expected multiplicative change in 2014/15",cex=0.8,line=1)
 abline(v=  ( c(100, 200,300, 400)), lty=2, col='grey40')
plot(f_grid, ppb10, col=red_c[1], lwd=1.5, type='l' , axes=F,ylab="",xlab="",yaxs='i', xpd=T, ylim=c(0,1), xlim=log( c(1,800)))
lines(f_grid, ppb50, col=red_c[2], xpd=T,lwd=1.5)
lines(f_grid, ppb100, col=red_c[3], xpd=T,lwd=1.5)
axis(1, padj=-0.1,lwd=0.5,at=log(c(0.1,1,10,50,100,1000)),  labels  = c(expression(10^-1), 1  , 10, 50,expression(10^2),expression(10^3) ) )
axis(2, las=2,lwd=0.5,at=c(0,0.5,1))
abline(h=c(0.25,0.5,0.75), lty=2, col='grey40')
#text(x= c(1.8, 3.5, 5.5), y=c(0.45,0.54,0.55), col= red_c, labels = c(">10", ">50", ">100") , cex=1)
text(x=   (c(1.8, 3.6, 5.8)), y=c(0.45,0.54,0.55), col= red_c, labels = c(">10", ">50", ">100 ppb") , cex=1)
mtext(1, line=1.3,text="baseline As in 2000 (ppb)",cex=0.8)
mtext(3, text="probablity of excessing safe levels \n in 15 years",cex=0.8, line=0.5)
abline(v= log( c(10,50, 100)), lty=2, col='grey40')
abline(v=  ( c(100,200,300,400)), lty=2, col='grey40')

plot(c(0,0), axes=F,ylab="",xlab="",xaxs='i', xpd=T, 
		 xlim=c(0.5, 2.2), ylim=c(30,75))
points(c(1,2), y=  c(mean(  exp(sss2$mu0)) , mean( exp(sss2$mu) )  ), pch=18, col="#8F2727", cex=2)
abline(h=c(30,40,50,60,70), col='gray40', lty=2)
lines(x=c(1,1), y= as.vector(   quantile(  exp(sss2$mu0), c(0.975, 0.025) )), lwd=1, col="#B97C7C" )
lines(x=c(1,1), y= as.vector( quantile(exp(sss2$mu0), c(0.75, 0.25) )), lwd=3, col="#8F2727" )
lines(x=c(2,2), y= as.vector( quantile(exp(sss2$mu), c(0.975, 0.025) )), lwd=1, col="#B97C7C" )
lines(x=c(2,2), y= as.vector( quantile(  exp(sss2$mu), c(0.75, 0.25) )), lwd=3, col="#8F2727" )
axis(2, las=2, at=c(30,50,70))
mtext(3, text="spatial mean",cex=0.8, line=0.5)
text(x= c(1.3, 1.7), y=c(53,41), col= 1, labels = c(" 2000", "14-15 ") , cex=1)
dev.off()


 



 