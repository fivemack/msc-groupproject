asteroid_diameter = 140
telescope_diameter = 0.7
asteroid_temperature = 240
asteroid_emissivity = 0.9

asteroid_power = stefan_boltzmann * asteroid_emissivity * pi * (asteroid_diameter/2)**2 * asteroid_temperature**4
print asteroid_power

