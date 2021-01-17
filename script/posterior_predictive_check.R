# Posterior predictive check
y21_median= apply(sss$y21, 2, median)
y22_median= apply(sss$y2_individual, 2, median)
n2=length(y21_median)
change=matrix(NA, 4000,n2)
y_grid= ( seq(0.001,1000, length.out=100))
y2_individual = beta_itercept + y21 + (beta_y*y21_exp+ beta_y_log*y21)  .* laplician+ beta_l* laplician + (sigma+sigma2)* obs_error_2; 
uu1=( y2_207 -y1_207) 
uu2=( y3_207 -y1_207) 
uu=( (y2_207 + y3_207)/2 -y1_207) 
change= sss$y2_individual -sss$y21
for(i in 1:n2){
	change[,i]=  	change[,i]+ sss$sigma*rnorm(4000) 
}
change_trim=matrix(NA, 4000, 271)
for( i in 1:4000){
	temp=sss$y21[i,] 
	id=sample(which(temp <range(y1_207)[2] & temp >range(y1_207)[1] ) , size=271)
  change_trim[i,]=change[i,id] 
}
mean(rowMeans(change_trim)< mean(uu2))
change2=change
for(i in 1:n2){
	yy=sss$y21[,i]
	change2[,i]= sss$beta_itercept+  (sss$beta_y  * exp(yy/2)   +sss$beta_y_log * yy   +  sss$beta_l) * sss$laplician[,i]  
}
change_trim2=change_trim
for( i in 1:4000){
	temp=sss$y21[i,] 
	id=sample(which(temp <range(y1_207)[2] & temp >range(y1_207)[1] ) , size=271)
	change_trim2[i,]=change2[i,id] 
}
pdf("~/Desktop/ppc.pdf", width=6.5, height=2)
layout(matrix(c(1,1,2,3,4,4,5,5),nrow=2), width = c(1,1,1),height = c(1,1))
par(oma=c(1.5,0,1.5,0), pty='m',mar=c(0.5,1,0.5,0.5) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.8, cex.lab=0.9, cex.main=0.9)   
 
hist(colMeans(change), breaks=seq(-3.4,3,by=0.2), xlim=c(-3.4,3),axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050",ylim=c(0,1000))
axis(1, at=c(-3,0,3), padj=-0.5,lwd=0.5)
mtext(3,  text="posterior mean of \n well-wise log change", cex=0.7, line=-0.4)

mtext(1,  text="log As change", cex=0.7, line=1)
abline(v=0,lty=2)
vv=y2_207-y1_207
hist(vv[vv>-3.2 & vv<3],  breaks=seq(-3.2,3,by=0.2), xlim=c(-3.2,3),axes = FALSE, main="", xlab="", ylab="", col=alpha(alpha=0.5,"forestgreen"), border = "forestgreen",ylim=c(0,100))
axis(1, at=c(-3,0,3), padj=-0.5,lwd=0.5)
mtext(3,  text="observed  well-wise \n log change 14", cex=0.7, line=0)

abline(v=0,lty=2)
vv=y3_207-y1_207
hist(vv[vv>-3.2 & vv<3],  breaks=seq(-3.2,3,by=0.2), xlim=c(-3.2,3),axes = FALSE, main="", xlab="", ylab="", col=alpha(alpha=0.5,"forestgreen"), border = "forestgreen",ylim=c(0,100))
axis(1, at=c(-3,0,3), padj=-0.5,lwd=0.5)
mtext(1,  text="log As change", cex=0.7, line=1)
abline(v=0,lty=2)
mtext(3,  text="log change 15", cex=0.7, line=-1)

hist(  rowMeans( change_trim), breaks=30, axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050",ylim=c(0,400)) 
abline(v=  (mean(uu1)), col="forestgreen", lwd=1)
abline(v=  (mean(uu2)), col="forestgreen", lwd=1)
abline(v=0, lty=2)
axis(1, at=c(-0.6,-0.3, 0), padj=-0.5,lwd=0.5)
mtext(1,  text="mean of log As change", cex=0.7, line=1)
mtext(3,  text="PPC. mean of 271 random wells \n p=0.07, 0.09 ", cex=0.7, line=0)


hist(  apply( change_trim2,1, sd ), breaks=30, axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050",ylim=c(0,500)) 
abline(v=  (sd(uu1)), col="forestgreen", lwd=1)
abline(v=  (sd(uu2)), col="forestgreen", lwd=1)
abline(v=0, lty=2)
axis(1, at=c(0.5,0.6,0.7,0.8), padj=-0.5,lwd=0.5)
mtext(3,  text="PPC. sd of 271 random wells \n p=0.28, 0.31 ", cex=0.7, line=0)
mtext(1,  text="sd. of log As change", cex=0.7, line=1)
dev.off()



change3=change
for(i in 1:n2){
	yy=sss$y21[,i] + (sss$sigma)* rnorm(4000)
	change3[,i]= exp(sss$y2_individual[,i]) - exp(yy)
}
change_trim3=change_trim
for( i in 1:4000){
	temp=sss$y21[i,] 
	id=sample(which(temp <range(y1_207)[2] & temp >range(y1_207)[1] ) , size=271)
	change_trim3[i,]=change3[i,id] 
}
linear_full=rowMeans(change3)
linear_inter=rowMeans(change4)
pdf("~/Desktop/change_all.pdf",width=6.5,height=2.4)
layout(matrix(c(1,1,2,2,3,3,4,5),nrow=2), width = c(1,1,1,1),height = c(1,1))
par(oma=c(1,2,2,0), pty='m',mar=c(1,1,1,0.5) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.9, cex.lab=0.9, cex.main=0.9)   
par(mar=c(0,1,0.5,0.5))
plot(c(0,0), axes=F,ylab="",xlab="",xaxs='i', xpd=T,  type='n',
		 xlim=c(0.5, 2.2), ylim=c(-23,3))
points(c(1,2), y=  c(mean( linear_full), mean( linear_inter )  ), pch=18, col="#8F2727", cex=2)
abline(h=c(0,-10,-20,-30), col='gray40', lty=2)
lines(x=c(1,1), y= as.vector(   quantile( linear_full, c(0.975, 0.025) )), lwd=1, col="#B97C7C" )
axis(2,   at=c(mean( linear_full),as.vector( quantile(linear_full, c(0.975, 0.025) ))), col="darkred", labels = round(c(mean( linear_full),as.vector( quantile(linear_full, c(0.975, 0.025) )))), las=2 )
lines(x=c(-4,1), y=rep(-7.221511, 2), col="#B97C7C", lty=2)
lines(x=c(-4,1), y=rep(-21.126504 , 2), col="#B97C7C", lty=2)
lines(x=c(-4,1), y=rep(-14.01132, 2), col="#B97C7C", lty=2)

lines(x=c(1,1), y= as.vector( quantile(linear_full, c(0.75, 0.25) )), lwd=3, col="#8F2727" )
lines(x=c(2,2), y= as.vector( quantile(linear_inter, c(0.975, 0.025) )), lwd=1, col="#B97C7C" )
lines(x=c(2,2), y= as.vector( quantile(linear_inter, c(0.75, 0.25) )), lwd=3, col="#8F2727" )
axis(2, las=2, at=c(0,-10,-20,-30))
mtext(3, text="average As change",cex=0.7, line=0.5)
text(x= c(1.25, 1.7), y=c(-22,-13), col= 1, labels = c("overall", "intercept") , cex=1)
mtext(2, line=1.8, text="As change (ppb)",cex=0.7)
par(mar=c(1,1,1,0.5))
vv=apply( change_trim3,1, mean )
hist(vv[vv<45], xlim=c(-85, 45), xaxs="i", breaks=seq(-85,45,by=5),  axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050",ylim=c(0,800))
lines( x=rep( mean( exp(y2_207)- exp(y1_207) ), 2), y=c(0,600)  ,  col='forestgreen')
#axis(1, at=c(-8, -12), labels = c("obs. 14", "obs. 15"), lwd=0.5, col="forestgreen")
text(x=c(-3, -20),  y=c(600, 630),  labels = c("14", "obs.\n 15"), col="forestgreen", cex=0.8)

lines( x=rep( mean( exp(y3_207)- exp(y1_207) ), 2), y=c(0,600),  col='forestgreen')
axis(1, at=c(0,-80,-40,40), padj=-0.5, lwd=0.5)
abline(v=c(0), col=1, lty=2)
mtext(3, line=-1, text="posterior predictive dis. \n of mean As change  when \n sample size = 271" ,cex=0.7)
mtext(1, line=1, text="As change (ppb)",cex=0.7)
axis(1, at=c(-400, -200,0, 200, 400), padj=-0.5, lwd=0.5)
vv=apply( change3,2, mean)
hist(vv[vv<400], xlim=c(-400, 400), xaxs="i", breaks=seq(-400,400,by=25),  axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050", ylim=c(0,2500))
abline(v=c(0), col=1, lty=2)
mtext(1, line=1, text="As change (ppb)",cex=0.7)
axis(1, at=c(-400, -200,0, 200, 400), padj=-0.5, lwd=0.5)
mtext(3, text ="posterior mean \n of well-wise changes",line=0.5,cex=0.7)
vv=exp(y3_207)- exp(y1_207)
hist( vv[vv<400], breaks=seq(-400,400,by=25), col=alpha("forestgreen", alpha=0.5), border = "forestgreen",xlab="", ylab="",axes = FALSE, main="", ylim=c(0,150))
abline(v=c(0), col=1, lty=2)
mtext(3, text="observed changes in \n validation set 2015\n mean = -12",line=-0.5,ce=0.7)
vv=exp(y2_207)- exp(y1_207)
hist( vv[vv<400], breaks=seq(-400,400,by=25), col=alpha("forestgreen", alpha=0.5), border = "forestgreen",xlab="", ylab="", main="",  axes = FALSE, ylim=c(0,150))
abline(v=c(0), col=1,lty=2)
mtext(3, text ="validation set 2014 \n mean = -8",line=-1,cex=0.7)
axis(1, at=c(-400, -200,0, 200, 400), padj=-0.5, lwd=0.5)
mtext(1, line=1, text="As change (ppb)",cex=0.7)
dev.off()
pal_lap2=rev( brewer.pal(n = 11, name = "PRGn"))
col_change= adjustcolor(  map2color(change_median_1,  pal_lap2 ,limits= c(-3,3)) , alpha.f = 0.4)
col_change2= adjustcolor(  map2color(y2_207-y1_207,  pal_lap2 ,limits= c(-3,3)) , alpha.f = 0.7)
col_change3= adjustcolor(  map2color(y3_207-y1_207,  pal_lap2 ,limits= c(-3,3)) , alpha.f = 0.7)
loc_x=data_repeated14_15$X2000
loc_y=data_repeated14_15$Y2000
loc_x=(loc_x-min(loc_x))/scale_x
loc_y=(loc_y-min(loc_y))/scale_x
loc_x=loc_x[complete.cases(y3)]
loc_y=loc_y[complete.cases(y3)]
pdf("~/Desktop/spatio_change.pdf",width=6.5,height=2.5)
layout(matrix(c(1,1,2,3,10,10, 4,7,5,8,6,9 ),nrow=2), width = c(2,1,0.2,1,1,1),height = c(1,1))
par( oma=c(0,0,0.5,0), pty='m',mar=c(0.5,0.3,0.5,0.3) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.6, cex.lab=0.9, cex.main=0.9, xpd=F) 
par( mar=c(2,2,1,0.2))
plot(x2,y2,
		 xlim=c(0,1), ylim=c(0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.5, col=col_change,  axes=F,ylab="",xlab="", xpd=T)
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
box(lwd=0.5)
axis(1, at=c(0,4000,8000)/scale_x, labels =c(0,4,8), padj=-0.2, lwd=0, lwd.tick=0.5, lty=5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =c(0,3,6), las=2 , lwd=0, lwd.tick=0.5, lty=5 )
mtext(3, text="posterior mean of well-wise As change",cex=0.6, line=0.5)
mtext(1, line=1, cex=0.5, text="East (km)")
mtext(2, line=1, cex=0.5, text="North (km)")

par( mar=c(0.5,0.5,1,0.2))

plot(loc_x,loc_y,
		 xlim=c(-0,1), ylim=c(-0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change2,  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
mtext(3, text="observed changed in \n validation set (2014)",cex=0.5, line=-0.1)

par( mar=c(0.5,0.5,0.5,0.2))
plot(loc_x,loc_y,
		 xlim=c(-0,1), ylim=c(-0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change3,  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
mtext(3, text=" validation set (2015)",cex=0.5, line=-0.1)

par( mar=c(0.5,0.5,1,0.2))

### change>2
 aa=2
id_t=which(change_median_1>log(aa))
id=id_t[(order(change_median_1[id_t]))]
plot(x2[id],y2[id],
		 xlim=c(0,1), ylim=c(0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.6, col=col_change[id],  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
mtext(3, text="only show points \n with change > 2",cex=0.5, line=-0.1)


id=which(y2_207-y1_207>log(aa))
plot(loc_x[id],loc_y[id],
		 xlim=c(0,1), ylim=c(0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change2[id],  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
mtext(3, text=" validation set (2014)",cex=0.5, line=-0.1)

id=which(y3_207-y1_207>log(aa))
plot(loc_x[id],loc_y[id],
		 xlim=c(0,1), ylim=c(0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change3[id],  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
mtext(3, text=" validation set (2015)",cex=0.5, line=-0.1)

### change<1/2
par( mar=c(0.5,0.5,0.5,0.2))

id_t=    which(change_median_1<log(1/aa))
id=id_t[rev(order(change_median_1[id_t]))]
plot(x2[id],y2[id],
		 xlim=c(0,1), ylim=c(0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.6, col=col_change[id],  axes=F,ylab="",xlab="")
mtext(3, text="change < 1/2",cex=0.5, line=-0.1)


abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
id=which(y2_207-y1_207<log(1/aa))
plot(loc_x[id],loc_y[id],
		 xlim=c(-0,1), ylim=c(-0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change2[id],  axes=F,ylab="",xlab="")
id=which(y3_207-y1_207<log(1/aa))
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)
plot(loc_x[id],loc_y[id],
		 xlim=c(-0,1.05), ylim=c(-0,0.73), xaxs='i',yaxs='i', pch=19, cex=0.8, col=col_change3[id],  axes=F,ylab="",xlab="")
abline(v=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
abline(h=seq(0,8000, by=1000)/scale_x,lty=2,col='grey')
axis(1, at=c(0,4000,8000)/scale_x, labels =rep(NA, 3), padj=-0.2, lwd=0, lwd.tick=0.5 )
axis(2, at=c(0,3000,6000)/scale_x, labels =rep(NA, 3), las=2 , lwd=0, lwd.tick=0.5)
box(lwd=0.5)



legend_image <- as.raster(matrix( rev(pal_lap2[1:11]), nc=1))
par(mar=c(0.3,0,0,0),pty='m')
plot(c(-1,4),c(-1,4),type = 'n', axes = F,xlab = '', ylab = '', main = "")
rasterImage(legend_image,xleft=0 , -1, 3,xright=1)
text(x=1.3, y = c(-1.1,1,3.1), labels = c(1/20,1,20),cex=0.8 ,xpd=T)
dev.off()


