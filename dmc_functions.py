# -*- coding: utf-8 -*-
import math
import random

def localEnergy(r1,r2,r12):
    return(-4 + potential(r1,r2,r12))

def potential(r1,r2,r12):
    return(1/(r12))

def displacement(tau):
    return( -2*tau + math.sqrt(tau)*random.random() )
    
def importanceSampling(r1,r2,r12,tau,referenceEnergy):
   energycomparison =  tau*(referenceEnergy - localEnergy(r1,r2,r12))
   return math.exp(energycomparison) 
