# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Monte Carlo Methods in Ab initio Quantum chemistry hammond lester reynolds
import math
import random
import matplotlib.pyplot as plt

def randomWalk(walkerList,tau): 
    for walker in walkerList:       
        for position in walker:          
            position+=displacement(tau)

def localEnergy(r1,r2,r12):
    return(-4 + potential(r1,r2,r12))
def potential(r1,r2,r12):
    return(1/(r12))
def displacement(tau):
    return( -2*tau + math.sqrt(tau)*random.random() )  
def importanceSampling(r1,r2,r12,tau,referenceEnergy):
   energycomparison =  tau*(referenceEnergy - localEnergy(r1,r2,r12))
   return math.exp(energycomparison) 
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
                
            variationalEnergy += localEnergy(r1,r2,r12)
        
        variationalTotalEnergy += variationalEnergy/2500
        variationalSquaredTotalEnergy+=(variationalEnergy/2500)**2
        
        totalWalker.append(r)
        
    referenceEnergy.append(variationalTotalEnergy/400)     
    standard_error = math.sqrt( ( (variationalSquaredTotalEnergy/400) - ( (variationalTotalEnergy/400)**2) ) /400)  
    print("Current variational Energy : "+str( variationalTotalEnergy / 400 ))
    print("Standard Error is " + str(standard_error))
tau_steps = [.001,.005,.01,.015,.02,.025,.03,.035,.04,.05,.06,.07,.08,.09,.1]
diffusion_energy = []
f = open('results2.txt','w')
for tau in tau_steps:
    diffusion_energy_in_each_block = []
    print("working on tau = "+str(tau))
    for block in range(100):
        totalWalker = []  
        referenceEnergy = []   
        metropolis(totalWalker,referenceEnergy)
        targetPopulation = 400
        mcheck = []
        for step in range(200):
            randomWalk(totalWalker,tau)
        energySum =0    
        for step in range(600):       
            createdWalkers =[]
            numWalkers =0
            localTotalEnergy = 0
            
            for walker in totalWalker:       
                numWalkers+=1
                
                for position in walker:          
                    position+=displacement(tau)
                    
                electron1 = math.sqrt( (walker[0]**2) + (walker[1]**2) +(walker[2]**2 ) )
                electron2 = math.sqrt( (walker[3]**2) + (walker[4]**2) +(walker[5]**2 ) ) 
                seperation = math.sqrt( ( ( walker[0] - walker[3] )**2)+ ( ( walker[1] - walker[4] )**2)+  ( ( walker[2] - walker[5] )**2) )
                
                localTotalEnergy += localEnergy(electron1,electron2,seperation)
               
                energyCheck = importanceSampling(electron1,electron2,seperation,tau,referenceEnergy[-1])
                m = math.floor(energyCheck+random.random())
                
                if m >=3:
                   for reproduce in range(3):
                        createdWalkers.append(walker)
                elif m >=1:
                    for reproduce in range(m):
                        createdWalkers.append(walker)
            
            totalWalker=createdWalkers  
#            if step % 3 == 0:
            averageEnergy = localTotalEnergy / numWalkers  
            newReferenceEnergy = averageEnergy - (math.log(numWalkers/targetPopulation) /tau)
            referenceEnergy.append(newReferenceEnergy)
            
            energySum+=averageEnergy
        diffusion_energy_in_each_block.append(energySum/600)   
        
    diffusionenergy = sum(diffusion_energy_in_each_block) / len(diffusion_energy_in_each_block)
    diffusion_energy.append(diffusionenergy)
    print ("Diffusion Energy for tau = "+str(tau)+ " is " + str( diffusionenergy)  )    
    f.write(str(tau) + "" + str(diffusionenergy) + "\n")
print("Done!") 
f.close()

plt.plot(tau_steps,diffusion_energy)
plt.xlabel('time step')
plt.ylabel('energy')
plt.show()

