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
 
def waveband_proportion(T,lo,hi):
 sz = 1e-8
 q=0
 reached_tail = False
 accumulated_power = 0
 last_power_in_region = -7
 while (not reached_tail):
  w = (q+0.5)*sz
  power_in_region = planck(w,T)*sz
  if (power_in_region < last_power_in_region and power_in_region < 1e-8):
   reached_tail = True
  accumulated_power = accumulated_power + power_in_region
  last_power_in_region = power_in_region
  q=q+1
 
 power_in_band = 0
 for q in range(int(lo/sz),int(hi/sz)):
  w = (q+0.5)*sz
  power_in_region = planck(w,T)*sz
  power_in_band = power_in_band + power_in_region
 
 return power_in_band / accumulated_power

for T in [40,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,3000,30000]:
 print "%s\t%s" % (T,waveband_proportion(T,8e-6,12e-6))
