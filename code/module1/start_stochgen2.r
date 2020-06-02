##############################################################
# This is starting script #2 to the module 
#'Stochastic effects on the genetic structure of populations'
##############################################################

# implements the basic model of stochastic population genetics with recombination

# frequency[1], frequency[2], frequency[3], frequency[4] denote the frequencies of genotypes ab, Ab,aB and AB respectively. The same applies to the fitnesses

#--set parameters
populationSize<-1000
mutationRate<-0.0001
recombinationRate<-0.1
selectionStrength<-0.05
mutationMatrix<-rbind(c((1-mutationRate)^2,mutationRate*(1-mutationRate),mutationRate*(1-mutationRate),mutationRate^2),c(mutationRate*(1-mutationRate),(1-mutationRate)^2,mutationRate^2,mutationRate*(1-mutationRate)),c(mutationRate*(1-mutationRate),mutationRate^2,(1-mutationRate)^2,mutationRate*(1-mutationRate)), c(mutationRate^2,mutationRate*(1-mutationRate),mutationRate*(1-mutationRate),(1-mutationRate)^2))
fitness<- c(1,1+selectionStrength ,1+selectionStrength,(1+selectionStrength)^2)
 
# initial condition 
frequency<-c(1,0,0,0)

numberOfGenerations<-1000

frequencies<-matrix(0,nrow=numberOfGenerations,ncol=4)
time<-1:numberOfGenerations


for(i in time){	
	##----mutation
	frequency<-(mutationMatrix %*% frequency)
	##----selection
	frequency<- frequency*fitness
	frequency<- frequency/sum(frequency)
	##---- recombination 
	LinkageDisequilibrium<- frequency[1]* frequency[4]-frequency[2]* frequency[3]
	frequency[1]<-frequency[1]-recombinationRate*LinkageDisequilibrium
	frequency[2]<-frequency[2]+recombinationRate*LinkageDisequilibrium
	frequency[3]<-frequency[3]+recombinationRate*LinkageDisequilibrium
	frequency[4]<-frequency[4]-recombinationRate*LinkageDisequilibrium
	##---sampling 
	frequency<- ((rmultinom(1,populationSize,frequency))/populationSize)
	
	##---saving values
	frequencies[i,]<-frequency
	}
	
##---print the allele frequencies

plot(frequencies[,4]~time,type="l",ylab="frequencies")
lines(frequencies[,3]~time ,col=2)
lines(frequencies[,2]~time ,col=3)
lines(frequencies[,1]~time ,col=4)

legend(numberOfGenerations*0.7,0.7, legend=c("AB","Ab","aB","ab"), col=1:4, lty=1)
