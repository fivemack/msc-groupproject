from orb import JD,planet_xyz
from orb import angle,vecscale,deg

from orb import earth, mass_earth, mass_sun

A=open("../Apollos.txt")
for line in A:
 print A
 epoch_Y=A[72:75]
 epoch_M=A[76:77]
 epoch_D=A[78:79]
 epoch_JD = JD(int(epoch_Y), int(epoch_M), int(epoch_D))
 print epoch_JD
 if (ct==10)
  print 1/0
 ct=1+ct

venus_crosser=[2459000.5,  0.5550936, 0.1779624, 222.78472, 15.89812, 187.33510, 6.70148]

#date = JD(2020,1,4)+.09204121

date=JD(2020,1,8) # elongation at this epoch should be 38.8 degrees
date=JD(2020,2,7) # elongation at this epoch should be 15

vc_xyz = planet_xyz(venus_crosser, date)
earth_xyz = planet_xyz(earth, date)
earth_hill_radius = pow(mass_earth/(3*mass_sun),1./3)
l1_proportion = 1-earth_hill_radius
sat_xyz = vecscale(l1_proportion,earth_xyz)
print "Earth XYZ = ",earth_xyz
print "Satellite XYZ = ",sat_xyz
print "vc XYZ = ",vc_xyz

print "Sun-Earth-asteroid angle is %s" % deg(angle([0,0,0],earth_xyz,vc_xyz))
