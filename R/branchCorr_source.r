#' Expit function
#'
#' Evaluates the expit function
#' @param x A number or vector
#' @return The value of expit(x)
#' @export
#' @examples
#' x = .25; expit(x)
#' x = 1:5; expit(x)
expit = function(x){
  return( exp(x)/(1+exp(x)) )
}

#' Forward simulate from a stochastic compartmental model
#'
#' Uses Gillespie forward simulation from a stochastic compartmental model, i.e. representing a hematopoietic tree.
#'
#' @param t.end The total time until end of simulation
#' @param initPopulation A vector of length equal to total number of compartments/types; entries contain the initial population of each type
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param maxEvents The maximum number of events to simulate
#' @return A vector containing the mature type populations at end of simulation period (time t=t.end).
#' @export
#' @examples
#' progStructure <- c(1,1,1,2,2)
#' initPopulation <- c(10,0,0,0,0,0,0,0)
#' rates <- c(.3,.2,.5,.06, .03, 2, 4, 5, 3, 1, .15,.5,.8,.1,.05)
#' simCompartments(5, initPopulation, rates,progStructure)
simCompartments <- function(t.end, initPopulation, rates, progStructure, maxEvents = 999999999999){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  pop <- initPopulation    # tracks current populations in each compartment
  t.cur <- 0               # current time
  i = 0
  while(i < maxEvents){
    i = i +1
    #print(pop)
    #for(i in 1:maxEvents){
    if( sum( pop < 0 ) > 0 ){
      print("Negative population?")
      return(-99)
    }
    if( sum(pop) == 0 ){
      #print("Extinction")
      return( rep(0,length(pop)) )   #extinction
    }
    #line up population multipliers with per-particle rate vector
    popVec <- c( pop[1], rep(pop[1], numProgs), pop[unique(progStructure)+1], pop[progStructure+1], tail(pop, matureTypes) )
    eventRates <- rates*popVec

    t.next <- rexp(1, sum(eventRates))
    t.cur <- t.cur + t.next

    if(t.cur > t.end){            #end of simulation period
      return(pop)
    }

    #easy to debug this way: can also use 'switch'
    decision <- sample( length(eventRates), 1, prob = eventRates/sum(eventRates) )
    #early fate decisions in hidden compartments
    if(decision == 1){
      pop[1] <- pop[1] + 1
    }
    for(j in 1:numProgs){
      if(decision == (j+1)){
        pop[1] <- pop[1]-1; pop[j+1] <- pop[j+1]+1 #hsc differentiation
      } else if(decision ==(j+1+numProgs) ){
        pop[j+1] <- pop[j+1] - 1 #prog death
      }
    }
    #mature cell production/death decisions:
    for(k in 1:matureTypes){
      if( decision == (k+1+2*numProgs) ){ #differentiation, go up 1
        pop[1+numProgs+k] <- pop[1+numProgs+k] + 1
      }else if ( decision == (k+1+2*numProgs+matureTypes) ){ #death
        pop[1+numProgs+k] <- pop[1+numProgs+k] - 1
      }
    }
  }
  return(-99)
}

#' Simulate discretely observed data from a stochastic compartmental model
#'
#' This code simulates data analogously to \code{\link{simCompartments}}, but records populations at a specified list of observation times.
#'
#' @param t.end The total time until end of simulation
#' @param initPopulation A vector of length equal to total number of compartments/types; entries contain the initial population of each type
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param obsTimes A vector containing the observation times
#' @param maxEvents The maximum number of events to simulate
#' @param vecOutput Logical, whether to return a vectorized representation of the matrix of observation times by population counts
#' @return A matrix of population sizes of each compartment at each observaiton time. Rows index cell type; columns index observation times. Returns
#' the vectorized form of thie matrix if vecOutput=TRUE.
#' @export
#' @examples
#' progStructure <- c(1,1,1,2,2)
#' initPopulation <- c(10,0,0,0,0,0,0,0)
#' obsTimes <- 1:5
#' rates <- c(.3,.2,.5,.06, .03, 2, 4, 5, 3, 1, .15,.5,.8,.1,.05)
#' simCompObserved(5, initPopulation, rates,progStructure,obsTimes)
simCompObserved <- function(t.end, initPopulation, rates, progStructure, obsTimes, vecOutput = FALSE, maxEvents = 999999999999){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  pop <- initPopulation    # tracks current populations in each compartment
  t.cur <- 0               # current time
  i = 0

  obs <- matrix(0, length(pop), length(obsTimes))   # Matrix, will hold all observations of process
  m <- 1; obsTime <- obsTimes[m]          # j will keep track of the index of which sample time to return

  while(i < maxEvents){
    i = i +1
    if( sum( pop < 0 ) > 0 ){
      print("Negative population?")
      return(-99)
    }

    #line up population multipliers with per-particle rate vector
    popVec <- c( pop[1], rep(pop[1], numProgs), pop[unique(progStructure)+1], pop[progStructure+1], tail(pop, matureTypes) )
    eventRates <- rates*popVec

    if( sum(pop) == 0 ){
      #print("Extinction")
      t.next <- 99999999
      #return( rep(0,4) )   #extinction
    } else {  t.next <- rexp(1, sum(eventRates)) }

    t.cur <- t.cur + t.next

    while( t.cur > obsTime){
      obs[,m] <- pop
      m <- m+1
      if(m < length(obsTimes) + 1 ){
        obsTime <- obsTimes[m]
      } else { obsTime <- 999999999}  #after we exhaust the list, just set samp to be very large
    }

    if(t.cur > t.end){            #end of simulation period
      if(t.cur > obsTime){
        obs[,m] <- pop
      }
      if( vecOutput){
        #return( c(obs[3,], obs[4,]) )
        return( as.vector( t( obs[ (2+numProgs):length(pop), ]) ) ) #add hsc + index 1, then skip over progenitors
      } else {
        return( obs ) }
    }

    decision <- sample( length(eventRates), 1, prob = eventRates/sum(eventRates) )

    #early fate decisions in hidden compartments
    if(decision == 1){
      pop[1] <- pop[1] + 1
    }
    for(j in 1:numProgs){
      if(decision == (j+1)){
        pop[1] <- pop[1]-1; pop[j+1] <- pop[j+1]+1 #hsc differentiation
      } else if(decision ==(j+1+numProgs) ){
        pop[j+1] <- pop[j+1] - 1 #prog death
      }
    }
    #mature cell production/death decisions:
    for(k in 1:matureTypes){
      if( decision == (k+1+2*numProgs) ){ #differentiation, go up 1
        pop[1+numProgs+k] <- pop[1+numProgs+k] + 1
      }else if ( decision == (k+1+2*numProgs+matureTypes) ){ #death
        pop[1+numProgs+k] <- pop[1+numProgs+k] - 1
      }
    }
  }
  return(-99)
}

#' Simulate discretely observed data from a stochastic compartmental model with error check
#'
#' Wrapper for \code{\link{simCompObserved}} that checks for no error code
#' @inheritParams simCompObserved
#' @export
sim.once <- function(t.end, initPopulation, rates, progStructure, obsTimes, vecOutput=FALSE, maxEvents = 999999999999){
  res = -99 # error catch
  while(res[1] == -99){
    res <- simCompObserved(t.end, initPopulation, rates, progStructure, obsTimes, vecOutput, maxEvents)  }
  return(res)
}

#' Multivariate hypergeometric sampling of data
#'
#' Samples hypergeometrically from a matrix containing discretely observed data from a stochastic compartmental model.
#'
#' @param data A matrix containing discretely observed data, in the format returned by \code{\link{simCompObserved}}
#' @param sampSize The sample size or number of draws n in a hypergeometric distribution
#' @export
hyperGeoSample <- function(data, sampSize){
  sampleMatrix <- matrix(NA, dim(data)[1], dim(data)[2])
  for(j in 1:dim(data)[2]){
    sampleMatrix[,j] <- samp_mvhGeo(data[,j], sampSize)
  }
  return(sampleMatrix)
}

#' Sample multivariate hypergeometric distribution
#'
#' This function generates a single random vector from the multivariate hypergeometric distribution.
#'
#' @param colorPops A vector containing the number of each ``color" or type of success in the population
#' @param sampSize A number, the size of the sample to be drawn
#' @return A vector drawn without replacement from colorPops
#' @export
#' @examples
#' samp_mvhGeo(c(100,500,200,300,300,500,50),1000)
samp_mvhGeo <- function(colorPops, sampSize){
    K <- length(colorPops)
    N <- sum(colorPops)
    if(sampSize > N){
      res <- colorPops
    }else{
      res <- rep(0, K)
      subSample <- table(sample(rep(1:K, colorPops), sampSize))
      res[as.numeric(names(subSample))] <- subSample
    }
    return(res)
}

####################################
######### Moment equations #########
####################################

#' This function computes the model-based mean of a mature cell type compartment, given rates and that the process begins with one initial
#' progenitor.
#' @param t The length of time
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param progType An index indicating the initial progenitor type
#' @param type An index indicating the mature cell type
#' @return The mean population of compartment `type' after time t
M_2x <- function(t, rates, progType, type, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)

  if(progStructure[type] != progType){ return(0) }

  nu <- rates[1 + 2*numProgs + type]; mu <- rates[1 + 2*numProgs + matureTypes + type];
  mu0 <- rates[1+numProgs+progType]
  return( nu*(exp(-mu0*t) - exp(-mu*t) )/(mu - mu0) )
}

#' Mean population starting from HSC
#'
#' This function computes the model-based mean of a mature cell type compartment, given rates and that the process begins with one initial
#' hematopoietic stem cell (compartment 1).
#' @param t The length of time
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param progType An index indicating the progenitor type that the mature cell type can descend from
#' @param type An index indicating the mature cell type
#' @return The mean population of compartment `type' after time t
M_1x <- function(t, rates, progType, type, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  lam = rates[1]; nu0 <- rates[1+progType]; mu0 <- rates[1+numProgs+progType]; nu <- rates[1 + 2*numProgs + type]
  mu <- rates[1 + 2*numProgs + type + matureTypes]; sumProg <- sum(rates[2:(1+numProgs)])
  return( exp(t*(lam - sumProg)) * (nu0*nu)/(mu-mu0) * ( exp(t*(sumProg - lam - mu0))/(sumProg - lam - mu0) - exp(t*(sumProg - lam - mu))/(sumProg - lam - mu)
                                                         + 1/(sumProg - lam - mu) - 1/(sumProg - lam - mu0) ) )
}

#' Second moments starting from HSC
#'
#' This function computes the model-based second moments, denoted U_mm|0 in the manuscript, of a mature cell type, given rates and that the process begins with one initial
#' hematopoietic stem cell (compartment 1).
#' @inheritParams M_1x
#' @return Value of second moment at time t
U_xx <- function(t, rates, progType, type, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  lam <- rates[1]; nu0 <- rates[1+progType]; mu0 <- rates[1+numProgs+progType]
  nu <- rates[1 + 2*numProgs + type]; mu <- rates[1 + 2*numProgs + matureTypes + type]
  sumProg <- sum(rates[2:(1+numProgs)])

  I1 <- 2*nu0*nu^2/(mu-mu0)*( (mu0-mu)*exp((sumProg-lam-mu0)*t)/(mu*(mu0-2*mu)*(sumProg-lam-mu0))
                              - exp((sumProg-lam-mu0-mu)*t)/(mu*(sumProg-lam-mu0-mu)) - exp((sumProg-lam -2*mu)*t)/((mu0-2*mu)*(sumProg-lam-2*mu))
                              + (mu-mu0)/(mu*(mu0-2*mu)*(sumProg-lam-mu0))  + 1/(mu*(sumProg-lam-mu0-mu)) + 1/((mu0-2*mu)*(sumProg-lam-2*mu)) )

  I2 <- 2*lam*nu0^2*nu^2/(mu-mu0)^2*( exp((sumProg-lam-2*mu0)*t)/((sumProg-lam-mu0)^2*(sumProg-lam-2*mu0))
                                      - 2*exp((sumProg-lam-mu0-mu)*t)/((sumProg-lam-mu0)*(sumProg-lam-mu)*(sumProg-lam-mu0-mu)) + 2*(mu0-mu)*exp(-mu0*t)/(mu0*(sumProg-lam-mu)*(sumProg-lam-mu0)^2)
                                      + exp((sumProg-lam-2*mu)*t)/((sumProg-lam-mu)^2*(sumProg-lam-2*mu)) + 2*(mu-mu0)*exp(-mu*t)/(mu*(sumProg-lam-mu)^2*(sumProg-lam-mu0))
                                      + (mu-mu0)^2*exp((lam-sumProg)*t)/((lam-sumProg)*(sumProg-lam-mu)^2*(sumProg-lam-mu0)^2) - 1/((sumProg-lam-mu0)^2*(sumProg-lam-2*mu0))
                                      + 2/((sumProg-lam-mu0)*(sumProg-lam-mu)*(sumProg-lam-mu0-mu)) - 2*(mu0-mu)/(mu0*(sumProg-lam-mu)*(sumProg-lam-mu0)^2)
                                      - 1/((sumProg-lam-mu)^2*(sumProg-lam-2*mu)) - 2*(mu-mu0)/(mu*(sumProg-lam-mu)^2*(sumProg-lam-mu0))
                                      - (mu-mu0)^2/((lam-sumProg)*(sumProg-lam-mu)^2*(sumProg-lam-mu0)^2) )

  return( exp((lam-sumProg)*t)* (I1+I2) )
}

#' Second moments starting from progenitor
#'
#' This function computes the model-based second moments, denoted U_mm|a in the manuscript, of a mature cell type, given rates and that the process begins with one initial
#' progenitor cell.
#' @inheritParams M_2x
#' @return Value of second moment at time t
V_xx <- function(t, rates, progType, type, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  if(progStructure[type] != progType){ return(0) }

  nu <- rates[1+2*numProgs+type]; mu <- rates[1+2*numProgs+matureTypes+type]; mu0 <- rates[1+numProgs+progType]
  return( 2*nu^2/(mu-mu0) * exp(-mu0*t) * ( (mu0-mu)/(mu*mu0 - 2*mu^2) - exp(-mu*t)/mu - exp( (mu0 - 2*mu)*t )/(mu0 - 2*mu) ) )
}

#' Second cross-moments starting from HSC
#'
#' This function computes the model-based second cross-moments, denoted U_mn|0 in the manuscript, between two mature cell types,
#' given rates and that the process begins with one initial hematopoietic stem cell (compartment 1).
#' @param t The length of time
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param progType1 An index indicating the progenitor type that mature cell 1 can descend from
#' @param progType2 An index indicating the progenitor type that mature cell 2 can descend from
#' @param type1 An index indicating the type of mature cell 1
#' @param type2 An index indicating the type of mature cell 2
#' @return Value of second cross-moment at time t
U_xy <- function(t,rates, progType1, progType2, type1, type2, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  lam <- rates[1]; nuProg1 <- rates[1+progType1]; nuProg2 <- rates[1+progType2]
  muProg1 <- rates[1+numProgs+progType1]; muProg2 <- rates[1+numProgs+progType2] #progenitor rates
  nu1 <- rates[1+2*numProgs+type1]; nu2 <- rates[1+2*numProgs+type2]  #final type rates
  mu1 <- rates[1+2*numProgs+type1+matureTypes]; mu2 <- rates[1+2*numProgs+type2+matureTypes]
  sumProg <- sum(rates[2:(1+numProgs)])

  if(progType1 != progType2){
    I1 <- 0
    I2 <- (2*lam*nuProg1*nuProg2*nu1*nu2/((mu1-muProg1)*(mu2-muProg2)) *
             ( exp((sumProg-lam-muProg1-muProg2)*t)/( (sumProg-lam-muProg1-muProg2)*(sumProg-lam-muProg1)*(sumProg-lam-muProg2) )
               - exp((sumProg-lam-muProg1-mu2)*t)/( (sumProg-lam-muProg1)*(sumProg-lam-mu2)*(sumProg-lam-muProg1-mu2) )
               + (muProg2-mu2)*exp(-muProg1*t)/(muProg1*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu2))
               - exp((sumProg-lam-muProg2-mu1)*t)/((sumProg-lam-muProg2)*(sumProg-lam-mu1)*(sumProg-lam-muProg2-mu1))
               + exp((sumProg-lam-mu1-mu2)*t)/((sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu1-mu2))
               + (mu2-muProg2)*exp(-mu1*t)/(mu1*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-muProg2))
               + (muProg1-mu1)*exp(-muProg2*t)/(muProg2*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu1))
               + (mu1-muProg1)*exp(-mu2*t)/(mu2*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-muProg1))
               + (mu1-muProg1)*(mu2-muProg2)*exp((lam-sumProg)*t)/((lam-sumProg)*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu1)*(sumProg-lam-mu2))

               #fix the following constant
               - 1/( (sumProg-lam-muProg1-muProg2)*(sumProg-lam-muProg1)*(sumProg-lam-muProg2) )+ 1/( (sumProg-lam-muProg1)*(sumProg-lam-mu2)*(sumProg-lam-muProg1-mu2) )
               - (muProg2-mu2)/(muProg1*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu2))
               + 1/((sumProg-lam-muProg2)*(sumProg-lam-mu1)*(sumProg-lam-muProg2-mu1))
               - 1/((sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu1-mu2)) + (muProg2-mu2)/(mu1*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-muProg2))
               + (mu1-muProg1)/(muProg2*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu1)) + (muProg1-mu1)/(mu2*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-muProg1))
               - (mu1-muProg1)*(mu2-muProg2)/((lam-sumProg)*(sumProg-lam-muProg1)*(sumProg-lam-muProg2)*(sumProg-lam-mu1)*(sumProg-lam-mu2)) ) )
  } else {
    mu0 <- muProg1
    nu0 <- nuProg1
    I1 <- (nu0*nu1*nu2/(mu2-mu0)*(  (mu0-mu2)*exp( (sumProg-lam-mu0)*t )/(mu1*(mu0-mu1-mu2)*(sumProg-lam-mu0) )
                                    - exp((sumProg-lam-mu1-mu0)*t)/(mu1*(sumProg-lam-mu1-mu0) )
                                    - exp((sumProg-lam-mu1-mu2)*t)/( (mu0-mu1-mu2)*(sumProg-lam-mu1-mu2) )
                                    + (mu2-mu0)/(mu1*(mu0-mu1-mu2)*(sumProg-lam-mu0)) + 1/(mu1*(sumProg-lam-mu1-mu0)) + 1/((mu0-mu1-mu2)*(sumProg-lam-mu1-mu2))  )
           + nu0*nu1*nu2/(mu1-mu0)*(  (mu0-mu1)*exp( (sumProg-lam-mu0)*t )/(mu2 * (mu0-mu1-mu2)*(sumProg-lam-mu0) )
                                      - exp( (sumProg-lam - mu2 - mu0)*t)/(mu2*(sumProg-lam - mu2 - mu0) )
                                      - exp((sumProg-lam - mu1 - mu2)*t)/( (mu0-mu1-mu2)*(sumProg-lam - mu1 - mu2) )
                                      + (mu1-mu0)/(mu2*(mu0-mu1-mu2)*(sumProg-lam-mu0)) + 1/(mu2*(sumProg-lam-mu2-mu0)) + 1/((mu0-mu1-mu2)*(sumProg-lam-mu1-mu2)) ))

    I2 <- (2*lam*nu0^2*nu1*nu2/((mu1-mu0)*(mu2-mu0)) * ( exp((sumProg-lam-2*mu0)*t)/( (sumProg-lam-2*mu0)*(sumProg-lam-mu0)^2 )
                                                         - exp((sumProg-lam-mu0-mu2)*t)/( (sumProg-lam-mu0)*(sumProg-lam-mu2)*(sumProg-lam-mu0-mu2) )
                                                         + (mu0-mu2)*exp(-mu0*t)/(mu0*(sumProg-lam-mu0)^2*(sumProg-lam-mu2)) - exp((sumProg-lam-mu0-mu1)*t)/((sumProg-lam-mu0)*(sumProg-lam-mu1)*(sumProg-lam-mu0-mu1))
                                                         + exp((sumProg-lam-mu1-mu2)*t)/((sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu1-mu2)) + (mu2-mu0)*exp(-mu1*t)/(mu1*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu0))
                                                         + (mu0-mu1)*exp(-mu0*t)/(mu0*(sumProg-lam-mu0)^2*(sumProg-lam-mu1)) + (mu1-mu0)*exp(-mu2*t)/(mu2*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu0))
                                                         + (mu1-mu0)*(mu2-mu0)*exp((lam-sumProg)*t)/((lam-sumProg)*(sumProg-lam-mu0)^2*(sumProg-lam-mu1)*(sumProg-lam-mu2))
                                                         - 1/((sumProg-lam-mu0)^2*(sumProg-lam-2*mu0)) + 1/( (sumProg-lam-mu0)*(sumProg-lam-mu2)*(sumProg-lam-mu0-mu2))
                                                         - (mu0-mu2)/(mu0*(sumProg-lam-mu0)^2*(sumProg-lam-mu2)) + 1/((sumProg-lam-mu0)*(sumProg-lam-mu1)*(sumProg-lam-mu0-mu1))
                                                         - 1/((sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu1-mu2)) + (mu0-mu2)/(mu1*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu0))
                                                         + (mu1-mu0)/(mu0*(sumProg-lam-mu0)^2*(sumProg-lam-mu1)) + (mu0-mu1)/(mu2*(sumProg-lam-mu1)*(sumProg-lam-mu2)*(sumProg-lam-mu0))
                                                         - (mu1-mu0)*(mu2-mu0)/((lam-sumProg)*(sumProg-lam-mu0)^2*(sumProg-lam-mu1)*(sumProg-lam-mu2))  ) )
  }
  return( exp((lam-sumProg)*t)*(I1 + I2) )
}


#' Second cross-moments starting from progenitor
#'
#' This function computes the model-based second cross-moments, denoted U_mn|a in the manuscript, between two mature cell types,
#' given rates and that the process begins with one progenitor cell.
#' @param t The length of time
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param progType1 An index indicating the progenitor type that mature cell 1 can descend from
#' @param progType2 An index indicating the progenitor type that mature cell 2 can descend from
#' @param type1 An index indicating the type of mature cell 1
#' @param type2 An index indicating the type of mature cell 2
#' @return Value of second cross-moment at time t
V_xy <- function(t,rates, progType1, progType2, type1,type2, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)
  muProg1 <- rates[1+numProgs+progType1]; muProg2 <- rates[1+numProgs+progType2]; #progenitor rates
  nu1 <- rates[1+2*numProgs+type1]; nu2 <- rates[1+2*numProgs+type2];  #final type rates
  mu1 <- rates[1+2*numProgs+type1+matureTypes]; mu2 <- rates[1+2*numProgs+type2+matureTypes]

  if(progType1 != progType2 || progType1 != progStructure[type1] || progType2 != progStructure[type2]){
    return(0)
  }
  else{
    mu0 <- muProg1
    return( nu1*nu2/(mu2-mu0) * exp(-mu0*t) * ( (mu0-mu2)/( mu1*(mu0-mu1-mu2) ) - exp(-mu1*t)/mu1 - exp( (mu0-mu1-mu2)*t)/(mu0-mu1-mu2) )
            + nu1*nu2/(mu1-mu0) * exp(-mu0*t) * ( (mu0-mu1)/( mu2*(mu0-mu1-mu2) ) - exp(-mu2*t)/mu2 - exp( (mu0-mu1-mu2)*t)/(mu0-mu1-mu2) ) )
  }
}

#' Model-based variance from HSC
#'
#' This function computes the model-based variance of a mature cell type given that the process begins with one HSC.
#' The function simply converts the raw moments U, M into the corresponding variance expression.
#' @inheritParams U_xx
#' @return Variance of specified mature cell compartment at time t
#' @export
Var_xFrom1 <- function(t,rates,progType,type,progStructure){
  return( U_xx(t,rates,progType,type,progStructure) + M_1x(t,rates,progType,type,progStructure) - M_1x(t,rates,progType,type,progStructure)^2 )
}

#' Model-based variance from progenitor
#'
#' This function computes the model-based variance of a mature cell type given that the process begins with one progenitor.
#' The function simply converts the raw moments U, M into the corresponding variance expression.
#' @inheritParams V_xx
#' @return Variance of specified mature cell compartment at time t
#' @export
Var_xFrom2 <- function(t,rates,progType,type,progStructure){
  return( V_xx(t,rates,progType,type,progStructure) + M_2x(t,rates,progType,type,progStructure) - M_2x(t,rates,progType,type,progStructure)^2 )
}

#' Model-based covariance from HSC
#'
#' This function computes the model-based covariance between two mature cell types given that the process begins with one HSC.
#' The function simply converts the raw moments U, M into the corresponding covariance expression.
#' @inheritParams U_xy
#' @return Covariance between specified mature type compartments at time t
#' @export
Cov_xyFrom1 <- function(t,rates, progType1, progType2, type1, type2, progStructure){
  return( U_xy(t,rates, progType1, progType2, type1, type2, progStructure) - M_1x(t,rates,progType1,type1, progStructure)*M_1x(t,rates,progType2,type2, progStructure))
}

#' Model-based covariance from progenitor
#'
#' This function computes the model-based covariance between two mature cell types given that the process begins with one progenitor.
#' The function simply converts the raw moments U, M into the corresponding covariance expression.
#' @inheritParams V_xy
#' @return Covariance between specified mature type compartments at time t
#' @export
Cov_xyFrom2 <- function(t,rates, progType1, progType2, type1, type2, progStructure){
  if(progType1!=progType2){return(0)}
  return( V_xy(t,rates, progType1, progType2, type1, type2, progStructure) - M_2x(t,rates,progType1,type1, progStructure)*M_2x(t,rates,progType2,type2, progStructure))
}

#' Model-based correlation from HSC
#'
#' This function computes the model-based correlation between two mature cell types given that the process begins with one HSC.
#' The function is basically a wrapper around those used to compute variance and covariances.
#' @inheritParams U_xy
#' @return Correlation between specified mature type compartments at time t
#' @export
Cor_xyFrom1 <- function(t,rates, progType1, progType2, type1, type2, progStructure){
  denom <- sqrt(Var_xFrom1(t,rates,progType1,type1,progStructure)*Var_xFrom1(t,rates,progType2,type2,progStructure))
  if(denom==0){return(0)}
  Cov_xyFrom1(t,rates, progType1, progType2, type1, type2, progStructure)/denom
}

#' Model-based correlation from progenitor
#'
#' This function computes the model-based correlation between two mature cell types given that the process begins with one progenitor.
#' The function is basically a wrapper around those used to compute variance and covariances.
#' @inheritParams V_xy
#' @return Correlation between specified mature type compartments at time t
#' @export
Cor_xyFrom2 <- function(t,rates, progType1, progType2, type1, type2, progStructure){
  if(progType1!=progType2){return(0)}
  denom <- sqrt(Var_xFrom2(t,rates,progType1,type1,progStructure)*Var_xFrom2(t,rates,progType2,type2,progStructure))
  if(denom==0){return(0)}
  Cov_xyFrom2(t,rates, progType1, progType2, type1, type2, progStructure)/denom
}

#' Marginalized model-based correlations
#'
#' Computes the marginalized correlation between two mature cell types, given the initial probability distribution.
#' @param t The time length
#' @param initProb vector containing initial probabilities of beginning in HSC comparmtent and each progenitor compartment
#' @param type1 Index of mature type 1
#' @param type2 Index of mature type 2
#' @param rates A rate vector containing self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @return The model-based marginalized correlation between types at time t
#' @export
marginalizedCor_xy <- function(t,rates, initProb, type1, type2, progStructure){
  progType1 <- progStructure[type1] #true corresponding progenitors according to model
  progType2 <- progStructure[type2]

  Cov <- initProb[1]^2*Cov_xyFrom1(t,rates, progType1, progType2, type1, type2, progStructure) +
    initProb[1]*(1- initProb[1])*U_xy(t,rates, progType1, progType2, type1, type2, progStructure) +
    initProb[1+progType1]^2*Cov_xyFrom2(t,rates, progType1, progType2, type1, type2, progStructure) +  #will be zero if don't share progenitor
    initProb[1+progType1]*(1-initProb[1+progType1])*V_xy(t,rates, progType1, progType2, type1, type2, progStructure) -
    initProb[1]*initProb[1+progType1]*M_2x(t,rates,progType1,type1, progStructure)*M_1x(t,rates,progType2,type2, progStructure) -
    initProb[1]*initProb[1+progType2]*M_2x(t,rates,progType2,type2, progStructure)*M_1x(t,rates,progType1,type1, progStructure) -
    (progType1!=progType2)*initProb[progType1+1]*initProb[progType2+1]*M_2x(t,rates,progType1,type1, progStructure)*M_2x(t,rates,progType2,type2, progStructure) #pairs only have one nonzero term here

  Var_x <- initProb[1]*( U_xx(t,rates,progType1,type1,progStructure) + M_1x(t,rates,progType1,type1, progStructure) ) -
    initProb[1]^2*( M_1x(t,rates,progType1,type1, progStructure)^2 ) +
    initProb[1+progType1]*(V_xx(t,rates,progType1,type1,progStructure) + M_2x(t,rates,progType1,type1, progStructure)) -
    initProb[1+progType1]^2*(M_2x(t,rates,progType1,type1, progStructure)^2) -
    2*initProb[1]*initProb[1+progType1]*M_1x(t,rates,progType1,type1, progStructure)*M_2x(t,rates,progType1,type1, progStructure)

  Var_y <- initProb[1]*( U_xx(t,rates,progType2,type2,progStructure) + M_1x(t,rates,progType2,type2, progStructure) ) -
    initProb[1]^2*( M_1x(t,rates,progType2,type2, progStructure)^2 ) +
    initProb[1+progType2]*(V_xx(t,rates,progType2,type2,progStructure) + M_2x(t,rates,progType2,type2, progStructure)) -
    initProb[1+progType2]^2*(M_2x(t,rates,progType2,type2, progStructure)^2) -
    2*initProb[1]*initProb[1+progType2]*M_1x(t,rates,progType2,type2, progStructure)*M_2x(t,rates,progType2,type2, progStructure)

  return( Cov/(sqrt(Var_x*Var_y)) )
}

#' Marginalized model-based correlations with sampling
#'
#' Computes the marginalized correlation between two mature cell types, given the initial probability distribution, incorporating
#' the effect of hypergeometric sampling.
#' @inheritParams marginalizedCor_xy
#' @param pop1 Sample size of mature type 1
#' @param pop2 Sample size of mature type 2
#' @param total1 Total population of mature type 1 cells
#' @param total2 Total population of mature type 2 cells
#' @return Correlation between types at time t after marginalizing over sampling distribution and initial distribution
#' @export
samplingCor_xy <- function(t,rates, initProb, type1, type2, pop1, pop2, total1, total2, progStructure ){
  progType1 <- progStructure[type1] #true corresponding progenitors according to model
  progType2 <- progStructure[type2]

  Cov <- initProb[1]^2*Cov_xyFrom1(t,rates, progType1, progType2, type1, type2, progStructure) +
    initProb[1]*(1- initProb[1])*U_xy(t,rates, progType1, progType2, type1, type2, progStructure) +
    initProb[1+progType1]^2*Cov_xyFrom2(t,rates, progType1, progType2, type1, type2, progStructure) +  #will be zero if don't share progenitor
    initProb[1+progType1]*(1-initProb[1+progType1])*V_xy(t,rates, progType1, progType2, type1, type2, progStructure) -
    initProb[1]*initProb[1+progType1]*M_2x(t,rates,progType1,type1, progStructure)*M_1x(t,rates,progType2,type2, progStructure) -
    initProb[1]*initProb[1+progType2]*M_2x(t,rates,progType2,type2, progStructure)*M_1x(t,rates,progType1,type1, progStructure) -
    (progType1!=progType2)*initProb[progType1+1]*initProb[progType2+1]*M_2x(t,rates,progType1,type1, progStructure)*M_2x(t,rates,progType2,type2, progStructure) #pairs only have one nonzero term here

  var_x1 <- initProb[1]*( U_xx(t,rates,progType1,type1,progStructure) + M_1x(t,rates,progType1,type1, progStructure) ) -
    initProb[1]^2*( M_1x(t,rates,progType1,type1, progStructure)^2 ) +
    initProb[1+progType1]*(V_xx(t,rates,progType1,type1,progStructure) + M_2x(t,rates,progType1,type1, progStructure)) -
    initProb[1+progType1]^2*(M_2x(t,rates,progType1,type1, progStructure)^2) -
    2*initProb[1]*initProb[1+progType1]*M_1x(t,rates,progType1,type1, progStructure)*M_2x(t,rates,progType1,type1, progStructure)

  var_x2 <- initProb[1]*( U_xx(t,rates,progType2,type2,progStructure) + M_1x(t,rates,progType2,type2, progStructure) ) -
    initProb[1]^2*( M_1x(t,rates,progType2,type2, progStructure)^2 ) +
    initProb[1+progType2]*(V_xx(t,rates,progType2,type2,progStructure) + M_2x(t,rates,progType2,type2, progStructure)) -
    initProb[1+progType2]^2*(M_2x(t,rates,progType2,type2, progStructure)^2) -
    2*initProb[1]*initProb[1+progType2]*M_1x(t,rates,progType2,type2, progStructure)*M_2x(t,rates,progType2,type2, progStructure)


  Mean_x1 <- initProb[1]*M_1x(t,rates,progType1,type1,progStructure) + initProb[1+progType1]*M_2x(t,rates,progType1,type1, progStructure)
  Mean_x2 <- initProb[1]*M_1x(t,rates,progType2,type2,progStructure) + initProb[1+progType2]*M_2x(t,rates,progType2,type2, progStructure)

  Cov <- pop1*pop2*Cov/(total1*total2) #- n*(N-n)*( CovX34 + MeanX3*MeanX4 )/( N^2*(N-1) )
  var_y1 <- pop1*(total1-pop1)*Mean_x1/( total1*(total1-1) ) - pop1*(total1-pop1)*(var_x1 + Mean_x1^2)/(total1^2*(total1-1)) + pop1^2*var_x1/total1^2
  var_y2 <- pop2*(total2-pop2)*Mean_x2/( total2*(total2-1) ) - pop2*(total2-pop2)*(var_x2 + Mean_x2^2)/(total2^2*(total2-1)) + pop2^2*var_x2/total2^2
  return( Cov/(sqrt(var_y1*var_y2)) )
}

#####################################################################
##### Loss functions, objectives for optimization and inference #####
#####################################################################

#' Compute observed pairwise correlations
#'
#' Returns the empirical pairwise correlations between each pair of types, given a dataset of observed counts
#' and a vector containing the corresponding observation times. The data matrix should be in the same format
#' as returned by \code{\link{simCompObserved}} with vecOutput=TRUE
#'
#' @param data A matrix of mature cell type counts, in format produced by \code{\link{simCompObserved}} with vecOutput=TRUE
#' @param obsTimes A vector of corresponding observation times
#' @return Matrix containing pairwise correlations at each observation time
#' @export
getObsCorr <- function(data,obsTimes){
  len <- length(obsTimes)
  numTypes <- dim(data)[2]/len
  pairs <- combn(seq(1:numTypes),2) - 1
  obsCorr <- matrix(NA, nrow = dim(pairs)[2], ncol = len ) #will store observed pairwise corr
  for(j in 1:dim(pairs)[2]){
    for(i in 1:len)
      obsCorr[j,i] <- cor(data[, (pairs[1,j]*len + i) ], data[, (pairs[2,j]*len+i) ]   )
  }
  return(obsCorr)
}

#' Correlation loss function objective
#'
#' Computes the objective function of the correlation-matching loss function estimator. Takes in model parameters and a matrix
#' of empirical correlations in the format produced by \code{\link{getObsCorr}}. This can then be plugged into a generic optimization
#' routine such as \code{nlminb}. Natural constraints such as positivity of parameters are enforced via a log-barrier penalty.
#'
#' @param par A vector containing the process rates followed by all but the first component of the initial distribution vector. The rates
#' should be ordered as before, beginning with self-renewal rates of first compartment, followed by differentiation rates of intermediate types,
#' death rates of intermediate types, differentiation rates of mature types, and death rates of mature types, i.e.
#' (renewHSC, diffProgA, diffProgB, ..., deathProgA, deathProgB, ... diffMature1, diffMature2, ..., deathMature1, deathMature2, ...)
#' @param obsTimes A vector of observation times
#' @param obsCorr Matrix of empirical correlations returned by \code{\link{getObsCorr}}
#' @param nSample A vector containing the sample sizes of each mature type as arguments to hypergeometric sampling; this is the number of draws n
#' @param nTotal A vector containing the total sizes of each mature population; this is the total population N in hypergeometric sampling
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @return Objective value
#' @export
samplingCorrObjective <- function( par, obsTimes, obsCorr, nSample, nTotal, progStructure){
  numProgs <- length(unique(progStructure)) #number of unique progenitors
  matureTypes <- length(progStructure)

  rates <- par[1:(length(par)-numProgs)]
  initProbTemp <- par[(length(par)-numProgs+1):length(par)]
  p1 <- 1/(1+sum(exp(initProbTemp)))
  initProb <- c(p1, p1*exp(initProbTemp)) #using multinomial logistic parametrization
  numTypes <- matureTypes
  pairs <- combn(seq(1:numTypes),2) #indices of final types, laid out in an order of all combinations

  Corr <- obsCorr #matrix of same size to hold analytical correlation expressions

  for( j in 1:dim(pairs)[2]){ #loop over number of pairwise combinations:
    Corr[j,] <- samplingCor_xy(obsTimes,rates,initProb, pairs[1,j], pairs[2,j],
                               nSample[pairs[1,j]], nSample[pairs[2,j]], nTotal[pairs[1,j],], nTotal[pairs[2,j],], progStructure)
  }
  obj <- Corr - obsCorr
  # we constrain rates to be positive, reserve growth rate to be positive (lam - sum of nu's) > 0
  penalty <- .01*log( rates[1] - sum(rates[2:(numProgs+1)])) #log barrier penalty
  return( sum( obj^2 ) - penalty ) #frobenius norm: kind of like GMM matrix norm with weight matrix W = Identity
}

#' Optimize correlation loss function
#'
#' This function uses the generic optimization package \code{nlminb} to optimize the objective function given by
#' \code{\link{samplingCorrObjective}}. Takes additional initial guess parameter and allows specification of true death rates,
#' which are often known and can be fixed throughout optimization.
#' @inheritParams samplingCorrObjective
#' @param initGuess Vector containing initial guess for par
#' @param trueDeaths Vector containing the true death rates to be fixed
#' @param max Max iterations, default to 5000
#' @return An \code{nlminb} object containing optimal parameters
#' @export
inferNLMINBSamplingCorr <- function( initGuess, obsTimes, obsCorr, nSample, nTotal, trueDeaths, progStructure, max=5000){   #trueDeaths=c(.26,.13,.11,.16,.09)

  numProgs <- length(unique(progStructure)) #number of unique progenitors
  #rates <- par[1:(length(par)-numProgs)]

  #constraints for nlminb
  low <- rep(0, length(initGuess)-numProgs)
  up <- rep(Inf, length(initGuess)-numProgs)

  #fix death rates
  low[(length(low)-length(trueDeaths)+1):length(low)] <- up[(length(up)-length(trueDeaths)+1):length(up)] <- trueDeaths
  #add initProb constraints (in current setup, unconstrained)
  low <- c(low, rep(-Inf,numProgs)); up <- c(up, rep(Inf,numProgs))

  return( nlminb(initGuess, samplingCorrObjective, obsTimes = obsTimes, obsCorr = obsCorr,
                 nSample = nSample, nTotal = nTotal, progStructure = progStructure,
                 lower = low, upper = up, control = list(rel.tol = 1e-15, iter.max = max)) )
}

#' Infer optimal rates via correlation loss function
#'
#' This function is a wrapper for \code{\link{inferNLMINBSamplingCorr}} that takes a dataset in the format produced by
#' \code{\link{simCompObserved}} and infers its most likely parameters. Computes the empirical correlations and optimizes
#' the correlation loss function over many random initializations, and returns the best estimate in terms of lowest objective
#' function value.
#'
#' @param obsTimes A vector of observation times
#' @param trueDeaths Vector containing the true death rates to be fixed
#' @param nSample A vector containing the sample sizes of each mature type as arguments to hypergeometric sampling; this is the number of draws n
#' @param nTotal A vector containing the total sizes of each mature population; this is the total population N in hypergeometric sampling
#' @param progStructure A vector of length equal to number of mature types whose i'th entry contains the corresponding hidden progenitor type from which mature type i descends
#' @param data A matrix of observed counts in the same format as produced by \code{\link{simCompObserved}}
#' @param numInits Number of random restarts
#' @param initMean The mean for the random initial parameters
#' @return Matrix whose rows contain the best three solutions
#' @export
#'
optimizeSamplingCorrNLMINB <- function(data, numInits, initMean, obsTimes, nSample, nTotal, trueDeaths, progStructure){
  len <- length(obsTimes)
  obsCorr <- getObsCorr(data,obsTimes)
  numProgs <- length(unique(progStructure)) #number of unique progenitors

  estimates <- matrix(NA,numInits,(length(initMean) + numProgs + 1))
  for(i in 1:numInits){
    tempRates <- runif(length(initMean), min = initMean - .4*abs(initMean), max = initMean + 1.3*abs(initMean))
    initGuess <- c( tempRates, runif(numProgs,1,10) ) #fix range for init guess for initprob
    #error handling:
    tryCatch({
      sol <- inferNLMINBSamplingCorr( initGuess, obsTimes, obsCorr, nSample, nTotal, trueDeaths, progStructure)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")} )
    param <- sol$par
    #label switching: let's put bigger of the two further diff rates  first wlog
    #     if(param[4] < param[5]){
    #       #if(param[6] < param[7]){
    #       temp <- param[4]; param[4] <- param[5]; param[5] <- temp
    #       temp <- param[6]; param[6] <- param[7]; param[7] <- temp
    #     }
    estimates[i,] <- c(param, sol$obj)
    est <- estimates[i,]
    initProbTemp <- est[(length(initMean)+1):(length(initMean)+numProgs)]
    p1 <- 1/(1+sum(exp(initProbTemp)))
    estimates[i,(length(initMean)+1):(length(initMean)+numProgs)] <- p1*exp(initProbTemp)
    estimatesBest <- estimates[which(estimates[,dim(estimates)[2]] %in% sort(estimates[,dim(estimates)[2]])[1:1]),] #best three for robustness
  }
  return(estimatesBest)
}
#' Generate random initial population vector
#'
#' This function converts an initial distribution to an initial population indicator vector. It samples from
#' the initial distribution, and returns an initial population vector with an indicator
#' of 1 in the compartment from which the process will begin, with zeros elsewhere.
#'
#' @param initProbs The initial distribution vector
#' @param totalTypes The total number of compartments or types
#' @return A vector of zeros in all but one entry, indicating the initial type
#' @export
randInit <- function(initProbs,totalTypes){
  initInd <- sample(1:length(initProbs),size=1,prob=initProbs)
  initPop <- rep(0,totalTypes)
  initPop[initInd] <- 1
  return(initPop)
}
