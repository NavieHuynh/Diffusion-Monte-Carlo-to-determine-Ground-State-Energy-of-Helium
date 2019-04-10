from Electron import ElectronSystem
import numpy as np
import datetime

class helium_monte_carlo_system:
    def __init__(self, vmc_electron_system={}, dmc_electron_system={},counter_step=0,counter_walkers=0,  counter_energy=0, vmc_energy=0,vmc_energy_sq=0,dmc_energy=0,dmc_energy_sq=0):
        self.vmc_electron_system = vmc_electron_system
        self.vmc_energy=vmc_energy
        self.vmc_stderr=vmc_energy_sq
        self.dmc_electron_system = dmc_electron_system
        self.dmc_walker_system = dmc_electron_system
        self.dmc_reference_energy=dmc_energy
        self.dmc_stderr=vmc_energy_sq
        self.dmc_final_energy=dmc_energy
        self.captured = {}
        self.counter_step= counter_step
        self.counter_energy = counter_energy
        self.counter_walkers= counter_walkers
        self.metropolis_algorithm()

    def reset_energy_count(self):
        self.counter_energy=0
    def reset_step_count(self):
        self.counter_step=0
    def reset_walker_count(self):
        self.counter_step=0

    def metropolis_algorithm(self,metropolis_steps=100,capture_last_energy=False, capture_last_position=False):
        for i in range(metropolis_steps):
            #will create new electron system if not existing in index
            if i not in self.vmc_electron_system: 
                self.vmc_electron_system[i] = ElectronSystem(np.random.uniform(-1,1,(6,)))

            s1= 4*(self.vmc_electron_system[i].r1new - self.vmc_electron_system[i].r1)
            s2= 4*(self.vmc_electron_system[i].r1new - self.vmc_electron_system[i].r2)
            ds = s1+s2
            
            if ds <= 0:
                self.vmc_electron_system[i].update_position()
                self.counter_step+=1
            elif np.random.rand() < np.exp(-ds):
                self.vmc_electron_system[i].update_position()
                self.counter_step+=1
            
            if(capture_last_energy):
                self.counter_energy += self.vmc_electron_system[i].calc_local_energy()
            
        if(capture_last_position):
               self.captured = self.vmc_electron_system[metropolis_steps - 1]
        # message="Acceptance rate = "+ str( 100*(self.counter_step / iterations) ) + "%\n"
        # print(message, end="\r")
        self.reset_step_count()

    def calc_vmc_energy(self,vmc_steps=100,metropolis_steps=100, dmc=False):
        vmc_energy = 0
        vmc_energy_sq= 0
        for i in range(vmc_steps):
            self.metropolis_algorithm(metropolis_steps, capture_last_energy=True, capture_last_position=dmc )
            self.dmc_electron_system[i] =self.captured
            vmc_energy+= self.counter_energy / metropolis_steps
            vmc_energy_sq+=( self.counter_energy / metropolis_steps ) ** 2
            self.reset_energy_count()
            message ="vmc_energy: " + str( 100*(i / vmc_steps)) + "% done"
            print(message, end="\r")
        self.vmc_energy = vmc_energy / vmc_steps
        self.vmc_stderr = np.sqrt((vmc_energy_sq / vmc_steps) - (self.vmc_energy ** 2 ) )
        message= "[%s]vmc_energy is %s with stderr = %s" %(datetime.datetime.now(), self.vmc_energy, self.vmc_stderr)
        print(message)

    def calc_dmc_energy(self,dmc_steps=100, timestep=.007, targetpopulation=400):
        local_dmc_energy=0
        dmc_energy_sq=0
       
        def dmc_walk(x, tau =timestep ):
            return ( x + (-2*tau + np.sqrt(tau)*np.random.rand()) )
       
        def dmc_system_walk(walk_func=dmc_walk, capture_last_energy=False):
            self.reset_energy_count()
            self.reset_walker_count()
            for system in self.dmc_electron_system:
                self.dmc_electron_system[system].walk(dmc_walk)
                self.dmc_electron_system[system].update_position()
                
                if(capture_last_energy):
                    self.counter_walkers+=1
                    self.counter_energy+= self.dmc_electron_system[system].calc_local_energy()
                    m = int(self.dmc_electron_system[system].calc_dmc_cutoff_value(reference_energy=self.dmc_reference_energy,tau=timestep) )
                    # print ("\ncurrent m value is %s"% m, end="\r")
                    if m>3:
                        for create_walkers in range(3):
                            self.dmc_walker_system[create_walkers] = self.dmc_electron_system[system]
                    elif m>=1:
                        for create_walkers in range(m):
                            self.dmc_walker_system[create_walkers]= self.dmc_electron_system[system]

        for i in range(dmc_steps):
            self.calc_vmc_energy(vmc_steps = 400,metropolis_steps=2500, dmc=True)
            self.dmc_reference_energy = self.vmc_energy
            for x in range(200):
                dmc_system_walk()

            local_total = 0
            for x in range(600):
                dmc_system_walk(capture_last_energy=True)                    
                self.dmc_electron_system =self.dmc_walker_system
                self.dmc_reference_energy= (self.counter_energy / self.counter_walkers) - np.log((self.counter_walkers/targetpopulation) / timestep)  
                local_total+=self.dmc_reference_energy
                message ="local_dmc_energy:%s | %s %% done" %(self.dmc_reference_energy ,100*(x / 600) )
                print(message, end="\r")
                    
            local_dmc_energy+= local_total / 600
            dmc_energy_sq+=( local_total/600 ) ** 2
            self.dmc_stderr = np.sqrt((dmc_energy_sq / dmc_steps) - (self.dmc_final_energy ** 2 ) ) / i
            message= "[%s][%s out of %s DMC Steps] has %s total walkers with average_dmc_energy: %s and stderr: %s " %(datetime.datetime.now(),i+1,dmc_steps,len(self.dmc_electron_system), local_dmc_energy/(i+1), self.dmc_stderr)
            print(message)

        self.dmc_final_energy = local_dmc_energy / dmc_steps
        
        message= "[%s] dmc_energy is %s with stderr = %s for timestep = %s" %(datetime.datetime.now(), self.dmc_final_energy, self.dmc_stderr, timestep)
        print(message)


            

        
