### This RScript reproduces analyses from Mason et al. manuscript on stable isotopes in horned larks of the imperial valley ###
### Contact mason@lsu.edu for further information or to report bugs / errors ###

setwd("~/.IVLarkIsotopesJoO")

### Package Requirements ###
#install.packages(c("gplots","SIBER","rjags","ellipse","car","plotrix"))
require(gplots)
require(SIBER)
require(rjags)
require(ellipse)
require(car)
require(plotrix)
require(nlme)
citation("nlme")
packageVersion("nlme")
require(viridis)
require(MuMIn)
citation("MuMIn")
packageVersion("MuMIn")

### Read in feather isotope data ###
feathers<-read.csv("./Appendices/IsotopeVoucherTable_v4.csv",row.names=2,stringsAsFactors=F)

feathers$breeding<-c("breeding","winter")[as.numeric(feathers$Month %in% c("10","11","12","1","2"))+1]
feathers$year<-as.numeric(feathers $Year)
feathers$histfactor<-feathers $Year > 1930 #Set up historical and contemporary factor
feathers$histfactor<-(as.numeric(feathers$histfactor)+1)
feathers$community<-1 #For SIBER ellipse analysis later on in script

feathers$Locality<-factor(feathers$Locality)

### How many samples for historical and contemporary feather series ###
table(feathers$histfactor)

### Compare 13C and 15N values between historical and contemporary ###
### Apply correction to C13 values (suess effects) ###
correction1<-(1960-as.numeric(feathers$Year))
correction1[correction1<0]<-0

correction2<-(2017-as.numeric(feathers$year))-correction1
correction2

feathers$δ13C.corrected <-as.numeric(feathers$δ13C) + (correction1 * -0.005) + (correction2 * -0.022)

### Get summary stats ###
feathers$δ15N<-as.numeric(feathers$δ15N)
feathers$δ13C <-as.numeric(feathers$δ13C)
feathers$histfactor <-as.factor(feathers$histfactor)

### NITROGEN ANALYSES ###
### Get means and standard deviations for N ###
### Contemporary 
mean(as.numeric(feathers$δ15N[feathers$histfactor==2]))
sd(feathers$δ15N[feathers$histfactor==2])

### Historical
mean(as.numeric(feathers$δ15N[feathers$histfactor==1]))
sd(feathers$δ15N[feathers$histfactor==1])

### LMM for bird differences in 15N ###
bird_N_lmm<-lme(fixed=δ15N~histfactor,random=~1|factor(Locality),data=feathers)
summary(bird_N_lmm)
r.squaredGLMM(bird_N_lmm)
hist(residuals(bird_N_lmm))

### CARBON ANALYSES ###
### Get means and standard deviations for C ###
### Contemporary 
mean(feathers$δ13C[feathers$histfactor==2])
sd(feathers$δ13C[feathers$histfactor==2])

### Historical
mean(feathers$δ13C[feathers$histfactor==1])
sd(feathers$δ13C[feathers$histfactor==1])

### LMM for bird differences in 13C ###
bird_C_lmm<-lme(δ13C~histfactor,random=~1|Locality,data=feathers)
summary(bird_C_lmm)
r.squaredGLMM(bird_C_lmm)
hist(residuals(bird_C_lmm))

### HYDROGEN ANALYSES ###
### Get means and standard deviations for N ###
### Contemporary 
mean(feathers$δ2H.corrected[feathers$histfactor==2])
sd(feathers$δ2H.corrected[feathers$histfactor==2])

### Historical
mean(feathers$δ2H.corrected[feathers$histfactor==1])
sd(feathers$δ2H.corrected[feathers$histfactor==1])

### LMM for bird differences in 2H ###
bird_H_lmm<-lme(δ2H.corrected~histfactor,random=~1|Locality,data=feathers)
summary(bird_H_lmm)
r.squaredGLMM(bird_H_lmm)
hist(residuals(bird_H_lmm))

d2h_hist_test<-t.test(feathers$δ2H.corrected[feathers$histfactor==1],mu=-65)
d2h_hist_test

d2h_cont_test<-t.test(feathers$δ2H.corrected[feathers$histfactor==2],mu=-65)
d2h_cont_test

### Generate plot of just bird data with ellipses and boxplots ###
hist_list<-split(feathers, feathers$histfactor)

## Add correction fo d2H values to account for fractionation between precipitation and feathers ###
feathers$δ2H.corrected<-feathers$δ2H + 25

## Create figure 3, which displays histograms of d2h values for contemporary and historical feathers ##
#png(file="./Fig3_DeuteriumHistogram_v10.png",width=3.25,height=3.5,units="in",res=300)
quartz(width=3.25,height=3.5)

par(mar=c(2.5,2,2,0.5))
par(lwd=0.5)
hist(feathers$δ2H.corrected[feathers$histfactor==2],col=viridis(2,alpha=0.6)[1],xlim=c(-90,-10),ylim=c(0,8),main="",xlab="",ylab="",axes=F,breaks=10)
hist(feathers$δ2H.corrected[feathers$histfactor==1],col=viridis(2,alpha=0.6)[2],xlim=c(-90,-10),ylim=c(0,8),main="",xlab="",ylab="",axes=F,breaks=10,add=T)
box(lwd=0.5)
axis(1,cex.axis=0.65,mgp=c(0,0.25,0),tck=-0.025,lwd=0.5)
mtext(expression(paste(delta^"2","H (‰)",sep="")),side=1,line=1.5,cex=0.8)
axis(2,cex.axis=0.65,mgp=c(0,0.25,0),tck=-0.025,lwd=0.5)
mtext("Frequency",side=2,line=1.15,cex=0.8)

par(xpd=F)
abline(v=-65,lwd=3,lty=2) #48 based on atmospheric
abline(v= mean(feathers$δ2H.corrected[feathers$histfactor==2]),lwd=3,lty=2,col=viridis(2)[1]) #48 based on atmospheric
abline(v= mean(feathers$δ2H.corrected[feathers$histfactor==1]),lwd=3,lty=2,col=viridis(2)[2]) #48 based on atmospheric

par(xpd=T)
lines(y=rep(8.6,2),x=c(mean(feathers$δ2H.corrected[feathers$histfactor==2]),-65),lwd=2,col=viridis(2)[1])
text(label="***",y=8.8,x=mean(c(mean(feathers$δ2H.corrected[feathers$histfactor==2]),-65)),col=viridis(2)[1])

lines(y=rep(8.6,2),x=c(mean(feathers$δ2H.corrected[feathers$histfactor==1]),-65),lwd=2,col=viridis(2)[2])
text(label="NS",y=8.8,x=mean(c(mean(feathers$δ2H.corrected[feathers$histfactor==1]),-65)),col=viridis(2)[2],cex=0.5)

lines(y=rep(9.2,2),x=c(mean(feathers$δ2H.corrected[feathers$histfactor==1]),mean(feathers$δ2H.corrected[feathers$histfactor==2])),lwd=2,col="black")
lines(y=c(9.1,9.3),x=rep(mean(feathers$δ2H.corrected[feathers$histfactor==1]),2),lwd=2,col=viridis(2)[2])
lines(y=c(9.1,9.3),x=rep(mean(feathers$δ2H.corrected[feathers$histfactor==2]),2),lwd=2,col=viridis(2)[1])
text(label="NS",y=9.4,x=mean(c(mean(feathers$δ2H.corrected[feathers$histfactor==1]),mean(feathers$δ2H.corrected[feathers$histfactor==2]))),col="black",cex=0.5)

legend("topright",legend=c("Contemporary","Historical"),pch=22,pt.bg=viridis(2,alpha=0.6),cex=0.65,pt.cex=1)
dev.off()

### Means, SDs, and two-sample t-test for plant differences in 15N and 13C ###
### PLANT ###
### Read in plant data ###
plants<-read.csv(file="./Appendices/PlantIsotopeData_v1.csv",row.names=NULL,stringsAsFactors=F)
plants$Locality<-as.numeric(factor(plants$Longitude))

## Nitrogen ##
##Contemporary / Agricultural
mean(plants$δ15N[plants$Species=="Taraxacum officinale"])
sd(plants$δ15N[plants$Species=="Taraxacum officinale"])

##Historical / Mesquite Scrub
mean(plants$δ15N[plants$Species=="Helianthus niveus tephrodes"])
sd(plants$δ15N[plants$Species=="Helianthus niveus tephrodes"])

##Plant N LMM
plant_N_lmm<-lme(δ15N~Species,random=~1|Locality,data=plants)
summary(plant_N_lmm)
r.squaredGLMM(plant_N_lmm) #Report in MS
hist(residuals(plant_N_lmm))

## Carbon ##
##Contemporary / Agricultural
mean(plants$δ13C[plants$Species=="Taraxacum officinale"])
sd(plants$δ13C[plants$Species=="Taraxacum officinale"])

##Historical / Mesquite Scrub
mean(plants$δ13C[plants$Species=="Helianthus niveus tephrodes"])
sd(plants$δ13C[plants$Species=="Helianthus niveus tephrodes"])

##Plant C LMM
plant_C_lmm<-lme(δ13C~Species,random=~1|Locality,data=plants)
summary(plant_C_lmm)
r.squaredGLMM(plant_C_lmm) #Report in MS
hist(residuals(plant_C_lmm))

### SOIL ###
soil<-read.csv(file="./Appendices/SoilIsotopeData_v1.csv",row.names=NULL,stringsAsFactors=F)
soil$Locality<-factor(as.numeric(factor(soil$Latitude)))

## Nitrogen ##
##Contemporary
mean(soil$δ15N[soil$Habitat.Type=="Agricultural Field"])
sd(soil$δ15N[soil $Habitat.Type=="Agricultural Field"])

##Historical
mean(soil$δ15N[soil $Habitat.Type=="Creosote Bush Scrub"])
sd(soil$δ15N[soil $Habitat.Type=="Creosote Bush Scrub"])

##Soil N LMM
soil_N_lmm<-lme(δ15N~Habitat.Type,random=~1|Locality,data=soil)
summary(soil_N_lmm)
r.squaredGLMM(soil_N_lmm) #Report in MS
hist(residuals(soil_N_lmm))

## Carbon ##
##Contemporary
mean(soil$δ13C[soil$Habitat.Type=="Agricultural Field"])
sd(soil$δ13C[soil $Habitat.Type=="Agricultural Field"])

##Historical
mean(soil$δ13C[soil $Habitat.Type=="Creosote Bush Scrub"])
sd(soil$δ13C[soil $Habitat.Type=="Creosote Bush Scrub"])

## Soil C LMM
soil_C_lmm<-lme(δ13C~Habitat.Type,random=~1|Locality,data=soil)
summary(soil_C_lmm)
r.squaredGLMM(soil_C_lmm) #Report in MS
hist(residuals(soil_C_lmm))

### SIBER ANALYSIS FOR FEATHERS ###
### Organize data frame to match createSiberObject demands ###
siber.data<-data.frame(N15=feathers$δ15N,C13= feathers$δ13C.corrected,community=feathers$community,group= feathers$histfactor,stringsAsFactors=F)

siber.data<-siber.data[,c(2,1,4,3)]
colnames(siber.data)<-c("iso1", "iso2", "group", "community")
siber.data<-siber.data[order(siber.data$community,siber.data$group),]

siber.data$iso1<-as.numeric(siber.data$iso1)
siber.data$iso2<-as.numeric(siber.data$iso2)
siber.data<-createSiberObject(siber.data)

### Try standard siber plot ###
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.90, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

quartz(width=8,height=4)
par(mfrow=c(1,2))
siberObject<-plotSiberObject(siber.data,
                  ax.pad = 2, 
                  hulls = F, community.hulls.args, 
                  ellipses = T, group.ellipses.args,
                  group.hulls = T, group.hull.args,
                  bty = "L",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~'\u2030'),
                  ylab = expression({delta}^15*N~'\u2030')
                  )
              
group.ML <- groupMetricsML(siber.data)
print(group.ML)

### Fitting Bayesian Models to the data ###
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.data, parms, priors)
head(siber.data)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)
head(SEA.B)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                xlab = c("Community | Group"),
                ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                bty = "L",
                las = 1,
                main = "SIBER ellipses on each group"
                )

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)
SEA.B.credibles

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)
  
names(SEA.B.modes)<-c("feather historical","feather contemp")
SEA.B.credibles

### Include Plant and Soil data in one big scatterplot ###
plants_list<-split(plants, plants$Species)
soil_list<-split(soil, soil$Habitat.Type)

#png(file="./Fig2_CompleteIsotopeBivariateplot_v10.png",width=3.25,height=3.5,units="in",res=300)

quartz(width=3.25,height=3.5)

par(xpd=F)
par(mgp=c(1.8,0.75,0))
par(mar=c(5,3,0.5,1))

plot(feathers$δ13C,feathers$δ15N,pch=21,col=rev(viridis(2))[as.numeric(feathers$histfactor)],bg="transparent",xlab=expression(paste(delta^"13","C (‰)",sep="")),ylab= expression(paste(delta^"15","N (‰)",sep="")),cex.lab=0.7,cex=0.7,cex.axis=0.7,xlim=range(c(soil$δ13C, plants$δ13C, feathers$δ13C.corrected)),ylim=range(c(soil$δ15N, plants$δ15N, feathers$δ15N)))

with(hist_list[[1]], dataEllipse(δ13C,δ15N, level = 0.90, add = TRUE,plot.points=F,col=viridis(2)[2],lwd=1,center.pch=F,center.cex=1))
with(hist_list[[2]], dataEllipse(δ13C,δ15N, level = 0.90, add = TRUE,plot.points=F,col=viridis(2)[1],lwd=1,center.pch=F,center.cex=1))

points(plants$δ13C[plants$Species=="Helianthus niveus tephrodes"], plants$δ15N[plants$Species=="Helianthus niveus tephrodes"],pch=23,col=viridis(2)[2],cex=0.7)
points(plants$δ13C[plants$Species=="Taraxacum officinale"], plants$δ15N[plants$Species=="Taraxacum officinale"],pch=23,col=viridis(2)[1],cex=0.7)

# with(plants_list[[1]], dataEllipse(δ13C,δ15N, level = 0.40, add = TRUE,plot.points=F,col="#bababa",lwd=1,center.pch=NA,center.cex=1))
# with(plants_list[[2]], dataEllipse(δ13C,δ15N, level = 0.40, add = TRUE,plot.points=F,col="#404040",lwd=1,center.pch=NA,center.cex=1))

points(soil$δ13C[soil$Habitat.Type=="Creosote Bush Scrub"], soil$δ15N[soil$Habitat.Type=="Creosote Bush Scrub"],pch=24,col=viridis(2)[2],cex=0.7)
points(soil$δ13C[soil$Habitat.Type=="Agricultural Field"], soil$δ15N[soil$Habitat.Type=="Agricultural Field"],pch=24,col=viridis(2)[1],cex=0.7)

# # with(soil_list[[1]], dataEllipse(δ13C,δ15N, level = 0.40, add = TRUE,plot.points=F,col="#404040",lwd=1,center.pch=NA,center.cex=1))
# with(soil_list[[2]], dataEllipse(δ13C,δ15N, level = 0.40, add = TRUE,plot.points=F,col="#bababa",lwd=1,center.pch=NA,center.cex=1))

par(xpd=NA)

points(x=-37,y=-10,pch=21,cex=0.6,col=viridis(2)[1],bg="transparent")
text(x=-36.5,y=-10,cex=0.4,label=paste("Contemporary larks (28)"),adj=c(0,0.5))

points(x=-37,y=-11.5,pch=21,cex=0.6,col=viridis(2)[2],bg="transparent")
text(x=-36.5,y=-11.5,cex=0.4,label=paste("Historical larks (16)"),adj=c(0,0.5))

points(x=-22.5,y=-10,pch=23,cex=0.6,col=viridis(2)[1],bg="transparent")
text(x=-22,y=-10,cex=0.4,label="Agricultural Asteraceae (15)",adj=c(0,0.5))

points(x=-22.5,y=-11.5,pch=23,cex=0.6,col=viridis(2)[2],bg="transparent")
text(x=-22,y=-11.5,cex=0.4,label="Desert Asteraceae (15)",adj=c(0,0.5))

points(x=-10.5,y=-10,pch=24,cex=0.6,col=viridis(2)[1],bg="transparent")
text(x=-10,y=-10,cex=0.4,label="Agricultural soil (15)",adj=c(0,0.5))

points(x=-10.5,y=-11.5,pch=24,cex=0.6,col=viridis(2)[2],bg="transparent")
text(x=-10,y=-11.5,cex=0.4,label="Desert soil (15)",adj=c(0,0.5))

dev.off()
