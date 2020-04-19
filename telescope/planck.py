from math import exp

planck_h = 6.626e-34
light_c = 2.998e8
boltzmann_k = 1.381e-23

def planck(wavelength,temperature):
 classical = 2*planck_h*light_c**2/wavelength**5
 log_quantum = (planck_h*light_c)/(wavelength*boltzmann_k*temperature)
 # exp(enormous) gives a math range error
 if (log_quantum > 709):
  return 0 
 quantum = exp(log_quantum)

 return classical/(quantum-1)
 
# this is 'energy per unit wavelength'

sz = 1e-8
q=0
reached_tail = False
accumulated_power = 0
last_power_in_region = -7
T=40
while (not reached_tail):
 w = (q+0.5)*sz
 power_in_region = planck(w,T)*sz
 if (power_in_region < last_power_in_region and power_in_region < 1e-8):
  reached_tail = True
 accumulated_power = accumulated_power + power_in_region
 last_power_in_region = power_in_region
 q=q+1
 
q_range = q
power_in_band = 0
for q in range(800,1200):
 w = (q+0.5)*sz
 power_in_region = planck(w,T)*sz
 power_in_band = power_in_band + power_in_region
 
print power_in_band / accumulated_power
