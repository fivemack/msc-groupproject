from orb import JD,planet_xyz
from orb import angle,distance,vecscale,deg

from orb import earth, mass_earth, mass_sun

import math
from collections import defaultdict

asteroids = []

A=open("../Apollos.txt")
ct=1
for line in A:
# print line
 if (ct>3):
  skip=False
  epoch_Y=line[72:76]
  epoch_M=line[76:78]
  epoch_D=line[78:80]
  if (epoch_Y==''):
   break
  epoch_JD = JD(int(epoch_Y), int(epoch_M), int(epoch_D))
  emoid = float(line[54:62])
  a = float(line[114:120])
  M = float(line[82:87])
  w = float(line[89:94])
  N = float(line[95:100])
  i = float(line[101:106])
  e = float(line[107:112])
 
  try:
    H = float(line[65:70])
  except:
    skip=True

  if (not skip):
   obj = [epoch_JD, a, e, M, i, w, N]
   name = line[27:36]
   asteroids = asteroids + [[name, obj, H, emoid]]
 ct=1+ct

print len(asteroids)

ever_visible = [False for r in range(len(asteroids))]

date = JD(2025,1,1) # first date we are interested in
earth_xyz = planet_xyz(earth, date)
earth_hill_radius = pow(mass_earth/(3*mass_sun),1./3)
l1_proportion = 1-earth_hill_radius
sat_xyz = vecscale(l1_proportion,earth_xyz)

num=0
happy=0
too_near_sun=0
too_near_earth=0
too_faint=0

dist_hist = {}
brightness_hist = defaultdict(int)

for A in asteroids:
 obj_xyz = planet_xyz(A[1], date)
# print "Earth XYZ = ",earth_xyz
# print "Satellite XYZ = ",sat_xyz
# print "Object XYZ = ",obj_xyz

 angle_ssa = deg(angle([0,0,0],sat_xyz,obj_xyz))
 angle_esa = 180-angle_ssa # L1 geometry
 dist = distance(sat_xyz, obj_xyz)
 
 dist_bin = int(10*dist)
 if (dist_bin in dist_hist):
  dist_hist[dist_bin] = 1+dist_hist[dist_bin]
 else:
  dist_hist[dist_bin] = 1

 kludge_brightness = A[2] + 5*math.log(dist)

 print A[2], dist, kludge_brightness
 brightness_bin = int(10*kludge_brightness)
 brightness_hist[brightness_bin] = 1+brightness_hist[brightness_bin]

 if (dist > 10.0):
  # this asteroid is a comet
  print A

 ok = True
 if (angle_ssa<45.0):
  ok = False
  too_near_sun = 1+too_near_sun

 if (angle_esa<45.0):
  ok=False
  too_near_earth = 1+too_near_earth

 if (kludge_brightness > 20.0):
  ok=False
  too_faint = 1+too_faint

 if (ok==True):
  happy=1+happy
 
 num=1+num

print num, too_near_sun, too_near_earth, too_faint, happy

