# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Monte Carlo Methods in Ab initio Quantum chemistry hammond lester reynolds
import math
import random
import matplotlib.pyplot as plt
import metropolis
import dmc_functions as dmc

def randomWalk(walkerList,tau): 
    for walker in walkerList:       
        for position in walker:          
            position+=dmc.displacement(tau)

tau_steps = [.001,.005,.01,.015,.02,.025,.03,.035,.04,.05,.06,.07,.08,.09,.1]
diffusion_energy = []

f = open('results2.txt','w')
for tau in tau_steps:
    diffusion_energy_in_each_block = []
    print("working on tau = "+str(tau))
    for block in range(100):
        totalWalker = []  
        referenceEnergy = []   
        metropolis.metropolis(totalWalker,referenceEnergy)
        
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
                    position+=dmc.displacement(tau)
                    
                electron1 = math.sqrt( (walker[0]**2) + (walker[1]**2) +(walker[2]**2 ) )
                electron2 = math.sqrt( (walker[3]**2) + (walker[4]**2) +(walker[5]**2 ) ) 
                seperation = math.sqrt( ( ( walker[0] - walker[3] )**2)+ ( ( walker[1] - walker[4] )**2)+  ( ( walker[2] - walker[5] )**2) )
                
                localTotalEnergy += dmc.localEnergy(electron1,electron2,seperation)
               
                energyCheck = dmc.importanceSampling(electron1,electron2,seperation,tau,referenceEnergy[-1])
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

