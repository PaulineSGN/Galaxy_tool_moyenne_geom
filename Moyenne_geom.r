#### loading required R libraries
#### chargement des packages R utilisés
library(gdata)
library(XLConnect)
library(rms)

###### overall parameters and settings
###### paramètres globaux utilisés

args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0)
{
    stop("This tool needs at least one argument")
}else{
    data <- args[1]
    nrep <- as.numeric(args[2])
    sep <- args[3]
    HR <- args[4]
     
}

if (HR =="false"){HR<-FALSE} else {HR<-TRUE}

###nrep: number of samples used to calculate geometric means
###nrep: nombre d'échantillons utilisés pour calculer les moyennes géométriques
nrep<-nrep

###### common functions
###### fonction utiles pour la suite

		convert.to.numeric<-function(x){
		t(apply(x,1,function(x){as.double(sub(" ","",as.character(x)))}))}

		
		### calculus of the logarithm of nrep geometric means, sampling based on a lognormal distribution with the same moments as the empirical ones (means & Ics)
			#to prevent negative values
		### calcul du logarithme de nrep moyennes géométriques, l'échantillonnage étant fait avec la distribution lognormale de mêmes moments que les momenst empriques (means et ICs)
			#pour éviter d'avoir des valeurs négatives

		lgeomean<-function(means,ICs,nrep)
		{#means: vector: mean estimates for the different categories 
		#ICs: vector: in proportion to the mean, difference between the extremum of the 95% confidence interval and the mean
		require(mvtnorm)
		#calculation of the parameters of the log normal distribution (on the log scale)
		#cf. http://127.0.0.1:26338/library/stats/html/Lognormal.html
		logsigma<-sqrt(log((ICs/qnorm(0.975)/means)^2+1))
		logmean<-log(means)-1/2*logsigma^2

		#gaussian sampling on the log scale then taking exponential
		temp<-exp(rmvnorm(nrep,mean=logmean,sigma=diag(logsigma*logsigma)))

		#taking geometric mean over categories, but kept on the log scale
		geomm.rep<-apply(temp,1,function(x){(mean(log(x),na.rm=TRUE))})
		#c(mean(geomm.rep),sd(geomm.rep))
		geomm.rep}



###### importation des données
###### importation of data
temp<-read.csv(file=data,sep=sep,header=HR,encoding="UTF-8")
data2008_2012<-temp[4:14,]
data2013_2017<-temp[21:31,]

meandata2008_2012<-convert.to.numeric(data2008_2012[,c(3,6,9)])
ICdata2008_2012<-convert.to.numeric(data2008_2012[,c(5,8,11)])
meandata2013_2017<-convert.to.numeric(data2013_2017[,c(3,6,9)])
ICdata2013_2017<-convert.to.numeric(data2013_2017[,c(5,8,11)])

####### code to calculate (nrep) logarithms of geometric means by region (Greco)
####### code pour calculer les nrep logarithmes de moyennes géométriques par région (GRECO)

set.seed(1)
#first period
#première période
rest2008_2012<-sapply(1:dim(data2008_2012)[1],function(region){lgeomean(meandata2008_2012[region,],ICdata2008_2012[region,],nrep)})

set.seed(3)
#first period but with different seed
#première période mais avec une graine différente
rest2008_2012_s3<-sapply(1:dim(data2008_2012)[1],function(region){lgeomean(meandata2008_2012[region,],ICdata2008_2012[region,],nrep)})

set.seed(2)
#second period
#seconde période
rest2013_2017<-sapply(1:dim(data2013_2017)[1],function(region){lgeomean(meandata2013_2017[region,],ICdata2013_2017[region,],nrep)})


####### code to summarize the above nrep logarithms of geometric means by region into the statistics of an overall geometric mean across regions, taking the first period as reference
###### code pour passer des nrep logarithmes de moyenne géométrique par région aux statistiques de la moyenne géométrique globale, en prennat la première période comme référence

#for the first period
#pour la première période
res2008_2012_scaled<-{temp<-apply(rest2008_2012_s3,1,function(x){mean(x)})-apply(rest2008_2012,1,function(x){mean(x)});c(mean(exp(temp)),sd(exp(temp)),quantile(exp(temp),prob=c(0.025,0.975)))}

#for the second period
#pour la seconde période
res2013_2017_scaled<-{temp<-apply(rest2013_2017,1,function(x){mean(x)})-apply(rest2008_2012,1,function(x){mean(x)});c(mean(exp(temp)),sd(exp(temp)),quantile(exp(temp),prob=c(0.025,0.975)))}

############### NATIONAL OUPUTS:
############### SORTIES NATIONALES:

res2008_2012_scaled
res2013_2017_scaled

write(res2008_2012_scaled, file = "res2008_2012_scaled.txt")
write(res2013_2017_scaled,file= "res2013_2017_scaled.txt")

############### REGIONAL OUPUTS:
############### SORTIES REGIONALES (GRECO):

regres2008_2012_scaled<-apply(rest2008_2012_s3-rest2008_2012,2,function(x){temp<-x;c(mean=mean(exp(temp)),sd=sd(exp(temp)),quantile(exp(temp),prob=c(0.025,0.975)))})
regres2013_2017_scaled<-apply(rest2013_2017-rest2008_2012,2,function(x){temp<-x;c(mean=mean(exp(temp)),sd=sd(exp(temp)),quantile(exp(temp),prob=c(0.025,0.975)))})
dimnames(regres2008_2012_scaled)[[2]]<-as.character(data2008_2012[,2])
dimnames(regres2013_2017_scaled)[[2]]<-as.character(data2013_2017[,2])

write(regres2008_2012_scaled, file = "regres2008_2012_scaled.txt")
write(regres2013_2017_scaled, file = "regres2013_2017_scaled.txt")
