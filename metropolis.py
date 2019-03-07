# -*- coding: utf-8 -*-

import math
import random
import dmc_functions as dmc

def metropolisStep(x):
    rmax = .7
    return x + rmax*(random.random() - .5)
def metropolis_step_list(arrayofparticles):
    return (list(map(metropolisStep,arrayofparticles)))
    
def metropolis(totalWalker,referenceEnergy): 
    counter = 0
    r = [.1,.1,.1,-.1,-.1,-.1]
    
    for walk in range(100): 
        rnew=metropolis_step_list(r)
        r1 = math.sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
        r2 = math.sqrt((r[3]**2)+(r[4]**2)+(r[5]**2))
        
        rnew1 = math.sqrt((rnew[0]**2)+(rnew[1]**2)+(rnew[2]**2))
        rnew2 = math.sqrt((rnew[3]**2)+(rnew[4]**2)+(rnew[5]**2))
        
        s1 = 4*(rnew1-r1)
        s2 = 4*(rnew2-r2)
        dS = s1+s2
        metropolisRndCheck = random.random()
        
        if dS <= 0:
            r=rnew
            counter+=1
        elif metropolisRndCheck < math.exp(-dS):
            r=rnew
            counter+=1
            
    print ("Acceptance rate =  " + str(counter) + "%" ) 
    variationalTotalEnergy = 0
    variationalSquaredTotalEnergy=0
   
    for i in range(400):
        variationalEnergy=0
        for j in range(2500):
            rnew=metropolis_step_list(r)
            r1 = math.sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
            r2 = math.sqrt((r[3]**2)+(r[4]**2)+(r[5]**2))
        
            rnew1 = math.sqrt((rnew[0]**2)+(rnew[1]**2)+(rnew[2]**2))
            rnew2 = math.sqrt((rnew[3]**2)+(rnew[4]**2)+(rnew[5]**2))
            
            r12 = math.sqrt( ((r[0]-r[3])**2)+((r[1]-r[4])**2)+((r[2]-r[5])**2) )
            s1 = 4*(rnew1-r1)
            s2 = 4*(rnew2-r2)
            dS = s1+s2
            metropolisRndCheck = random.random()
            
            if dS <= 0:
                r=rnew
            elif metropolisRndCheck < math.exp(-dS):
                r=rnew
                
            variationalEnergy += dmc.localEnergy(r1,r2,r12)
        
        variationalTotalEnergy += variationalEnergy/2500
        variationalSquaredTotalEnergy+=(variationalEnergy/2500)**2
        
        totalWalker.append(r)
        
    referenceEnergy.append(variationalTotalEnergy/400)     
    standard_error = math.sqrt( ( (variationalSquaredTotalEnergy/400) - ( (variationalTotalEnergy/400)**2) ) /400)  
    print("Current variational Energy : "+str( variationalTotalEnergy / 400 ))
    print("Standard Error is " + str(standard_error))
