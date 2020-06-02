##############################################################
# This is starting script #1 to the module 
#'Stochastic effects on the genetic structure of populations'
##############################################################

# implements the basic model of stochastic population genetics

#---set parameters

populationSize<-10000 # corresponds to N in the reader
mutationRate<-0.001   # corresponds to m
selectionStrength<-0  # corresponds to s

# Note: this module uses long, self-explanatory names for the parameters in the scripts.
# This is a good strategy in programming, but also a matter of personal taste.

#--- set initial condition
frequency1<-0.1  # frequency of the allele to be followed (A)
	
numberOfGenerations<-40000

frequencies<-numeric(numberOfGenerations)
time<-1:numberOfGenerations

for(i in time){	
	##----mutation
	frequency1<- (frequency1*(1-mutationRate) + mutationRate *(1-frequency1))
	##----selection
	frequency1<- ( frequency1 *(1-selectionStrength) / ( frequency1 *(1-selectionStrength) +(1-frequency1)) )
	##---sampling 
	frequency1<- ((rbinom(1,populationSize,frequency1))/populationSize)
	
	##---saving values
	frequencies[i]<-frequency1
	}

plot(frequencies~time,type="l")

# Note: the last command is equivalent to 'plot(time,frequencies,type="l")'