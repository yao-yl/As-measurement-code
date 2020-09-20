# model for the full data
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data_repeated308=read.csv("./data/4. repeated308.csv")
data_repeated14_15=read.csv("./data/2014vs15.csv",sep=',') 
data=read.csv("./data/1. survey2000-01.csv")
data2=read.csv("./data/2. survey2012-13.csv")
 
kit2=as.integer(factor(  data2$Arsenic_kit, levels=sort(unique(data2$Arsenic_kit ))))
kit2_numberic=  data2$Arsenic_kit
kit2_numberic[kit2_numberic==0]=0.001 # to avoid numerical issue in graphing
individual_y1=data$Arsenic_lab


library(splines2)
x1=data$X
x2=data2$X
y1=data$Y
y2=data2$Y
scale_x=max(x1, x2)-min(x1, x2)
mmx=min(x1, x2)
mmy=min(y1, y2)
x1= (x1-	mmx)/scale_x
x2= (x2-	mmx)/scale_x 
y1= (y1-	mmy)/scale_x
y2= (y2-	mmy)/scale_x 
 knots_x=seq(max(min(x1),min(x2)) ,min(max(x1),max(x2)), length.out =30)
knots_y=seq(max(min(y1),min(y2)),min(max(y1),max(y2)), by=knots_x[2]-knots_x[1])

 

Bx1 <- t(bSpline(x1, knots=knots_x, degree=3,    Boundary.knots = c(-0.05,1.05),  intercept = FALSE)) 
Bx2 <- t(bSpline(x2, knots=knots_x, degree=3, Boundary.knots = c(-0.05,1.05),  intercept = FALSE)) 
By1 <- t(bSpline(y1, knots=knots_y, degree=3,  Boundary.knots = c(-0.04,0.73),intercept = FALSE)) 
By2 <- t(bSpline(y2, knots=knots_y, degree=3,  Boundary.knots = c(-0.04,0.73),intercept = FALSE)) 
Bdx1 <- t(dbs(x=x1,derivs=2, knots=knots_x,   Boundary.knots = c(-0.05,1.05), degree=3, intercept = FALSE)) 
Bdx2 <- t(dbs(x=x2,derivs=2,  knots=knots_x,  Boundary.knots =  c(-0.05,1.05), degree=3, intercept = FALSE))
Bdy1 <- t(dbs(x=y1,derivs=2,  knots=knots_y,  Boundary.knots = c(-0.04,0.73),  degree=3, intercept = FALSE)) 
max(Bdy1)
Bdy2 <- t(dbs(x=y2,derivs=2,  knots=knots_y,  Boundary.knots = c(-0.04,0.73), degree=3, intercept = FALSE)) 
#dim(Bx1)=  33 4574
#dim(By1)=  33 4574

# trim knots 
assign_mat1=matrix(1:(nrow(Bx1)*nrow(By1)), nrow(Bx1), nrow(By1), byrow = TRUE)
B1=Laplacian1=matrix(NA, nrow(Bx1)*nrow(By1), ncol(Bx1) )
B2=Laplacian2=matrix(NA, nrow(Bx2)*nrow(By2), ncol(Bx2) )
id_keep=rep(0,(nrow(Bx1)*nrow(By1)) )
for(i in 1:nrow(Bx1)){
	flag=rep(0,nrow(By1))
	for(j in 1:nrow(By1))
	{
		id=assign_mat1[i,j]
		B1[id, ]= Bx1[i,]*By1[j,] 
		B2[id, ]= Bx2[i,]*By2[j,] 
		Laplacian1[id,]= Bdx1[i,]*By1[j,] + Bx1[i,]*Bdy1[j,] 
		Laplacian2[id,]= Bdx2[i,]*By2[j,] + Bx2[i,]*Bdy2[j,] 
		flag[j]=sum(B1[id, ])
	}
	if(sum(flag)>0)
	id_keep [assign_mat1[i, min (which(flag!=0))]:assign_mat1[i, max(which(flag!=0) )]]=1
}
trim_index=which( id_keep==1)
num_basis=length(trim_index)
B1_trim=B1[trim_index, ]
B2_trim=B2[trim_index, ]
Laplacian2_trim=Laplacian2[trim_index, ]/1000

col_mat=matrix( adjustcolor(  map2color(log(individual_y1),col_palette,limits=log(c(0.1,1000))), alpha.f = 0.3))
col_mat2=matrix( adjustcolor(  map2color(log(kit2_numberic),col_palette,limits=log(c(0.1,1000))), alpha.f = 0.3))
pdf("~/Desktop/spline_basis.pdf", height=2.5, width=7)
layout(matrix(c(1:3, 4,4,4,5,5,5),nrow=3), width = c(1,1,1),height = c(1))
par(oma=c(1.5,0.3,1,0), pty='m',mar=c(1,2,1,0.5) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=1, cex.lab=0.6, cex.main=0.7) 
x_grid=seq(0,1,0.005)
knots_x2=seq(0.1,0.9,0.1)
cc=brewer.pal("Paired", n=11)
Bx_grid<- t(bSpline(x_grid, knots=knots_x2, degree=3, intercept = FALSE)) 
plot(x_grid,Bx_grid[1,], type='n', ylim=c(0,1),axes=F,ylab="",xlab="")
for( i in 1:(dim(Bx_grid)[1]-1))
	lines(x_grid,Bx_grid[i,], col=cc[i], lwd=1)
axis(2, at=0,  las=2,lwd=0.5)
box(lwd=0.5, bty='l')
mtext(3, text="basis function",cex=0.8,  line=-1)
Bx_grid<- t(dbs(x_grid, knots=knots_x2, degree=3, derivs = 1,  intercept = FALSE)) 
plot(x_grid,Bx_grid[1,], type='n', ylim=c(-20,20),axes=F,ylab="",xlab="")
for( i in 1:(dim(Bx_grid)[1]-1))
	lines(x_grid,Bx_grid[i,], col=cc[i], lwd=1)
axis(2, at=0,  las=2,lwd=0.5)
		 box(lwd=0.5, bty='l') 
Bx_grid<- t(dbs(x_grid, knots=knots_x2, degree=3, derivs = 2,  intercept = FALSE)) 
mtext(3, text="1st derivative",cex=0.8,  line=-1)
plot(x_grid,Bx_grid[1,], type='n', ylim=c(-500,500),axes=F,ylab="",xlab="")
for( i in 1:(dim(Bx_grid)[1]-1))
	lines(x_grid,Bx_grid[i,], col=cc[i], lwd=1)
axis(2, at=0,  las=2,lwd=0.5)
box(lwd=0.5, bty='l') 
mtext(3, text="2rd derivative",cex=0.8, line=-1)



plot(x1,y1,xlim=c(-0.05,1.05), ylim=c(-0.04,0.73), xaxs='i',yaxs='i', pch=19, cex=0.5, col=col_mat,  axes=F,ylab="",xlab="")
for(i in 1:nrow(Bx1))
	for(j in 1:nrow(By1)){
		if(id_keep[assign_mat1[i,j]]==1)
			points(knots_x[i],knots_y[j], pch=4, cex=0.5)
	}
axis(1, at=c(0,4000,8000)/scale_x, labels =c(0,4,8), padj=-0.2, lwd=0, lwd.tick=0.5, lty=5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =c(0,3,6), las=2 , lwd=0, lwd.tick=0.5, lty=5 )
mtext(3, text="2000",cex=0.8, line=0.5)
box(lwd=0.8, lty=5)
mtext(1, line=1, cex=0.7, text="East (km)")
mtext(2, line=1, cex=0.7, text="North (km)")
plot(x2,y2,xlim=c(-0.05,1.05), ylim=c(-0.04,0.73), xaxs='i',yaxs='i', pch=19, cex=0.5, col=col_mat2,  axes=F,ylab="",xlab="")
for(i in 1:nrow(Bx1))
	for(j in 1:nrow(By1)){
		if(id_keep[assign_mat1[i,j]]==1)
			points(knots_x[i],knots_y[j], pch=4, cex=0.5)
	}
axis(1, at=c(0,4000,8000)/scale_x, labels =c(0,4,8), padj=-0.2, lwd=0, lwd.tick=0.5, lty=5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =c(0,3,6), las=2 , lwd=0, lwd.tick=0.5, lty=5 )
mtext(3, text="2014/15",cex=0.8, line=0.5)
legend(x="bottomright", pch=c(4,20, NA), lty=c(NA,NA,5),  legend = c("spline knots", "samples", "boundary"), cex=0.7,  box.lwd =0, xpd=T)
box(lwd=0.8, lty=5)
mtext(1, line=1, cex=0.7, text="East (km)")
mtext(2, line=1, cex=0.7, text="North (km)")
dev.off()



#### calibration 
cal=read.csv("intercal.csv")
head(cal)
cal_kit=cal$kit
cal_lab=cal$lab
cal_kit[cal_kit==999]=1000
data_cal=data.frame(kit=as.integer(factor(cal_kit)), cal_lab, log_lab=log(cal_lab) )
data_cal=data_cal[data_cal$cal_lab!=0,]
n_cal=nrow(data_cal)

stan_cal_data=list(
	w=9,
	N_cal=n_cal, 
	kit=data_cal$kit,
	log_lab=data_cal$log_lab
)
fit_cal=stan("./code/sep_cal.stan", data=stan_cal_data, iter = 2000, chains=3)
cal_samples=extract(fit_cal)


logistic_beta_mean=mean(cal_samples$beta)
logistic_beta_sd=sd(cal_samples$beta)
c_mean=colMeans(cal_samples$c)
c_sd=apply(cal_samples$c, 2, sd)

stan_data=list( n1=length(individual_y1),n2=length(kit2),
								w=9, 
								individual_y1=log(individual_y1),
								z2= kit2,  
								c=c_mean,
								logistic_beta=logistic_beta_mean,
								num_basis=dim(B1_trim)[1],
								B1=t(B1_trim),
								B2=t(B2_trim),
								Bd2=t(Laplacian2_trim),
								B1_m=sum(B1_trim!=0),
								B2_m=sum(B2_trim!=0),
								Bd2_m=sum(Laplacian2_trim!=0),
					      num_knots_2=3,
								knots=c(2,3.5,4.5))

fit_spline=stan("spline.stan", data=stan_data, init=ttt,
								iter = 2000, chains=4)

sss=extract(fit_spline)

coef_mat=matrix(NA, 4000,100) # coef before delta
y_grid= log( seq(1,1000, length.out=100))
for(i in 1:length(y_grid)){
		yy= y_grid[i]
		coef_mat[,i]=   (sss$beta_y  * exp(yy/2)   +sss$beta_y_log * yy   +  sss$beta_l) 
}
laplician_grid=seq(-3,3, length.out=100)
y_grid2=log(c(100,250))
change_lap=array(NA, c(3,2,4000,100))

for(i in 1:2)
	for(j in 1:100){
	yy= y_grid2[i]
change_lap[i,1,,j]= sss$beta_itercept+  (sss$beta_y  * exp(yy/2)   +sss$beta_y_log * yy   +  sss$beta_l) *laplician_grid[j] + sss$sigma2*rnorm(4000) 
change_lap[i,2,,j]= sss$beta_itercept+  (sss$beta_y  * exp(yy/2)   +sss$beta_y_log * yy   +  sss$beta_l) *laplician_grid[j] + (sss$sigma2+  sss$sigma) *rnorm(4000) 
	}

 pdf("~/Desktop/mixture_coef_log.pdf", width=7, height=2.2)
 layout(matrix(c(1,1,2,4,3,5),nrow=2), width = c(2,1,1),height = c(1,1))
 par( oma=c(1.5,1.5,1.5,0), pty='m',mar=c(1,1,1,4) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.9, cex.lab=0.6, cex.main=0.7, yaxs='i' )
 
 y_grid1=y_grid
plot( (y_grid1),  colMeans(coef_mat),  axes=F,ylab="",xlab="",type='l',col="#8F2727",lwd=2, ylim=c(-1,2.5) , xlim= c(-0.1,log(1000)), xaxs='i', yaxs='i' )
 abline(h=c(1,2,-1, -0.5), lty=2, col='grey40')
 axis(2, las=2,lwd=0.5,at=c(1,2,0,-1))
 abline(h=0,col= 'grey30',  lwd=1 )
 zz=apply(coef_mat, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 1000)), labels =c(1, 10,100, 1000)  )
 mtext(3, text="diffusion effect: \n coefficien on Laplacian", cex=0.8)
 mtext(1, text="well As (ppb) in 2000", cex=0.8, line=1)
 mtext(2, text="regression coefficient", cex=0.8, line=1)
 box(lwd=0.5,bty='l')

 
 
 
 par(mar=c(1,1,1,1))
 xxx=change_lap[2,1,,]
 plot(laplician_grid, colMeans(xxx), type='l', lwd=1, col="darkred",axes=F,ylab="",xlab="", ylim=c(log(1/4), log(4)))
 mtext(2, text="If y1\n = 250", las=2,cex=0.7,line=1.2)
 mtext(3, text="without observational noise", cex=0.7)
 box(lwd=0.5,bty='l')
 abline(h=log(c(1/4,1/2,2 ,4)), lty=2, col='grey60')
 axis(2, las=2,lwd=0.5,at=log(c(1/4, 1/2,  1,2, 4)),  labels  = c("1/4" , "1/2"  ,"1"  ,2,  4  ) )
 abline(h=0,col= 'darkgreen',  lwd=0.8 )
 zz=apply(xxx, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 lines(laplician_grid, colMeans(xxx), lwd=1, col="darkred")
 box(lwd=0.5,bty='l')
 axis(1, at=c(-3,0,3))
 
 
 
 xxx=change_lap[2,2,,]
 plot(laplician_grid, colMeans(xxx), type='l', lwd=1, col="darkred",axes=F,ylab="",xlab="", ylim=c(log(1/30), log(30)))
 box(lwd=0.5,bty='l')
 abline(h=log(c(1/25, 1/5,  1,5, 25)), lty=2, col='grey60')
 axis(2, las=2,lwd=0.5,at=log(c(1/25, 1/5,  1,5, 25)),  labels  = c("1/25" , "1/5"  ,"1"  ,5,  25  ) )
 zz=apply(xxx, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 abline(h=0,col= 'darkgreen',  lwd=0.8 )
 lines(laplician_grid, colMeans(xxx), lwd=1, col="darkred")
 box(lwd=0.5,bty='l')
 axis(1, at=c(-3,0,3))
 mtext(3, text="with observational noise", cex=0.7)
 mtext(3, text="predicted multiplicative changes from y1 to y2             ", cex=0.77 ,outer=T, adj=1,  padj = 0)
 
 xxx=change_lap[1,1,,]
 plot(laplician_grid, colMeans(xxx), type='l', lwd=1, col="darkred",axes=F,ylab="",xlab="", ylim=c(log(1/4), log(4)))
 box(lwd=0.5,bty='l')
 abline(h=log(c(1/4,1/2,2 ,4)), lty=2, col='grey60')
 axis(2, las=2,lwd=0.5,at=log(c(1/4, 1/2,  1,2, 4)),  labels  = c("1/4" , "1/2"  ,"1"  ,2,  4  ) )
 abline(h=0,col= 'darkgreen',  lwd=0.8 )
 zz=apply(xxx, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 lines(laplician_grid, colMeans(xxx), lwd=1, col="darkred")
 box(lwd=0.5,bty='l')
 axis(1, at=c(-3,0,3))
 mtext(2, text="If y1\n = 100", las=2,cex=0.7,line=1.2)
 mtext(1, text="Laplacian", cex=0.7,line=1)
 
 
 
 xxx=change_lap[1,2,,]
 plot(laplician_grid, colMeans(xxx), type='l', lwd=1, col="darkred",axes=F,ylab="",xlab="", ylim=c(log(1/30), log(30)))
 box(lwd=0.5,bty='l')
 abline(h=log(c(1/25, 1/5,  1,5, 25)), lty=2, col='grey60')
 axis(2, las=2,lwd=0.5,at=log(c(1/25, 1/5,  1,5, 25)),  labels  = c("1/25" , "1/5"  ,"1"  ,5,  25  ) )
 zz=apply(xxx, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(laplician_grid,rev(laplician_grid)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 abline(h=0,col= 'darkgreen',  lwd=0.8 )
 lines(laplician_grid, colMeans(xxx), lwd=1, col="darkred")
 box(lwd=0.5,bty='l')
 axis(1, at=c(-3,0,3))
 mtext(1, text="Laplacian", cex=0.7,line=1)
 
 
 
 
 
 dev.off()




 
 pdf("~/Desktop/mixture_coef_log_2.pdf", width=6.8, height=2)
 layout(matrix(c(1,1,2,2,3,4,5,6),nrow=2), width = c(1,1,0.9,0.7),height = c(1.1,1))
 
 par( oma=c(1.5,1.5,1.5,0), pty='m',mar=c(1,1,1,1) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.9, cex.lab=0.6, cex.main=0.7, yaxs='i' )
 
 y_grid1=y_grid
 plot( (y_grid1),  colMeans(coef_mat),  axes=F,ylab="",xlab="",type='l',col="#8F2727",lwd=2, ylim=c(-1,1.5) , xlim= c(-0.1,log(1000)), xaxs='i', yaxs='i' )
 abline(h=c(1,2,-1, -0.5), lty=2, col='grey40')
 axis(2, las=2,lwd=0.5,at=c(1,2,0,-1))
 abline(h=0,col= 'grey30',  lwd=1 )
 zz=apply(coef_mat, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 1000)), labels =c(1, 10,100, 1000)  )
 mtext(3, text="regress on y", cex=0.8)
 mtext(3, text="coefficien on Laplacian", cex=0.9,line=0.5, outer = T, adj=0.2,padj=0.7)
 mtext(1, text="well As (ppb) in 2000", cex=0.7, line=1)
 mtext(2, text="regression coefficient", cex=0.8, line=1)
 box(lwd=0.5,bty='l')
 
 
 
 y_grid1=y_grid
 plot( (y_grid1),  colMeans(coef_mat2),  axes=F,ylab="",xlab="",type='l',col="#8F2727",lwd=2, ylim=c(-1,1.5) , xlim= c(-0.1,log(1000)), xaxs='i', yaxs='i' )
 abline(h=c(1,2,-1, -0.5), lty=2, col='grey40')
 axis(2, las=2,lwd=0.5,at=c(1,2,0,-1))
 abline(h=0,col= 'grey30',  lwd=1 )
 zz=apply(coef_mat2, 2, quantile, c(0.975, 0.75, 0.25, 0.025))
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[1,],rev(zz[4,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.4), xpd=T)
 polygon(x=c(y_grid1,rev(y_grid1)), y=c(zz[2,],rev(zz[3,] )  ), border=NA,col=  adjustcolor( "#B97C7C", alpha.f = 0.45), xpd=T)
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 1000)), labels =c(1, 10,100, 1000)  )
 mtext(3, text="constant modeling", cex=0.8)
 mtext(1, text="well As (ppb) in 2000", cex=0.7, line=1)
 box(lwd=0.5,bty='l')
 par(mar=c(1,5,0.5,0))
 hist(vv, breaks=seq(0.5,6,by=0.2), ylim=c(0,6e6), axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050", xlim=c(0,log(350)))
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 300)), labels =c(1, 10,100, 300) )
 axis(1, padj=0,lwd=0.5,at=log( c( 300)), labels =c( 300) , cex=0.8)
 mtext(3, text="post.dis. of \n baseline of y", cex=0.8,line=-0.5)
 mtext(2, text="all draws", las=2, cex=0.7)
 hist(colMeans(vv), breaks=seq(0.5,6,by=0.2),  axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050", xlim=c(0,log(350)))
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 300)), labels =c(1, 10,100, 300) )
 mtext(2, text="well-wise\n mean", las=2, cex=0.7)
 axis(1, padj=0,lwd=0.5,at=log( c( 300)), labels =c( 300) , cex=0.8)
 par(mar=c(0.5,1,0.5,0))
 hist(vvv, breaks=seq(-8,8.8,by=0.4), ylim=c(0,6e6), axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050", xlim=c(-3, 8.5))
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 1000)), labels =c(1, 10,100, 1000) )
 mtext(3, text="baseline +\n obs.noise", cex=0.8,line=-0.5)
 hist(colMeans(vvv), breaks=seq(-8,8.8,by=0.4), axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050", xlim=c(-3, 8.5))
 axis(1, padj=-0.5,lwd=0.5,at=log( c(1, 10,100, 1000)), labels =c(1, 10,100, 1000) )
 dev.off()

