from orb import JD,planet_xyz
from orb import angle,distance,vecscale,deg

from orb import earth, mass_earth, mass_sun

import math
from collections import defaultdict

asteroid_files=["../Apollos.txt","../Amors.txt","../Atens.txt"]

asteroids = []

for F in asteroid_files:
 A=open(F)
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

   # ignore intrinsically tiny objects
   if (H>25.0):
    skip=True

   if (not skip):
    obj = [epoch_JD, a, e, M, i, w, N]
    name = line[27:36]
    asteroids = asteroids + [[name, obj, H, emoid]]
  ct=1+ct

print len(asteroids)

ever_visible = [False for r in range(len(asteroids))]
weeks_visible=[0 for r in range(len(asteroids))]
minimum_distance=[999.99 for r in range(len(asteroids))]

d0 = JD(2025,1,1) # first date we are interested in
earth_hill_radius = pow(mass_earth/(3*mass_sun),1./3)
l1_proportion = 1-earth_hill_radius
J2000 = JD(1999,12,31)

venus=[J2000, 0.723332,0.006772,50.115, 3.39458, 54.884, 76.680]
venus_sat=[J2000, 0.723332,0.006772,50.115+120, 3.39458, 54.884, 76.680]
inner_sat=[J2000, 0.4, 0.01, 123, 0, 102.93, 0]
halo_sat=[J2000, 1.0,0.05, 100.47, 0, 102.93, 0]

full_dists=[[] for i in range(len(asteroids))]

# mission duration = 5 years, call it 160 megaseconds
# go through in hundred-kilosecond chunks

for kilosec in range(1600):
 date = d0 + (100000.0/86400.0)*kilosec
 earth_xyz = planet_xyz(earth, date)
 #sat_xyz = planet_xyz(venus_sat, date)
 sat_xyz = vecscale(l1_proportion,earth_xyz)
 #sat_xyz = planet_xyz(halo_sat, date)

 num=0
 happy=0
 too_near_sun=0
 too_near_earth=0
 too_faint=0

 dist_hist = {}
 brightness_hist = defaultdict(int)

 for Ai in range(len(asteroids)):
  A=asteroids[Ai]
  obj_xyz = planet_xyz(A[1], date)
# print "Earth XYZ = ",earth_xyz
# print "Satellite XYZ = ",sat_xyz
# print "Object XYZ = ",obj_xyz

  angle_ssa = deg(angle([0,0,0],sat_xyz,obj_xyz))
  angle_esa = 180-angle_ssa # L1 geometry
  dist = distance(sat_xyz, obj_xyz)
  full_dists[Ai].append(dist)
  if (dist < minimum_distance[Ai]):
   minimum_distance[Ai]=dist
 
#  dist_bin = int(10*dist)
#  if (dist_bin in dist_hist):
#   dist_hist[dist_bin] = 1+dist_hist[dist_bin]
#  else:
#   dist_hist[dist_bin] = 1

  kludge_brightness = A[2] + 5*math.log(dist)

#  print A[2], dist, kludge_brightness
#  brightness_bin = int(10*kludge_brightness)
#  brightness_hist[brightness_bin] = 1+brightness_hist[brightness_bin]

  #if (dist > 10.0):
   # this asteroid is a comet
   #print A

  ok = True
  if (angle_ssa<45.0):
   ok = False
   too_near_sun = 1+too_near_sun

  if (angle_esa<45.0):
   ok=False
   too_near_earth = 1+too_near_earth

  if (kludge_brightness > 22.0):
   ok=False
   too_faint = 1+too_faint

  if (ok==True):
   ever_visible[Ai] = True
   weeks_visible[Ai]=1+weeks_visible[Ai]
   happy=1+happy
 
  num=1+num

 count_ever_happy = sum([1 for i in ever_visible if i==True])
 print kilosec*100, num, too_near_sun, too_near_earth, too_faint, happy, count_ever_happy

 if (kilosec%100==99):
  print "Distance histogram at t=%s" % (100000*kilosec)
  hist_mindist=defaultdict(int)
  for i in minimum_distance:
   ix = int(10*i)
   hist_mindist[ix]=1+hist_mindist[ix]
  S=sorted(minimum_distance)

  print hist_mindist

  for u in range(1,20):
   qq = 0.05*u
   print ("%.2f"%qq)," ",S[int(len(S)*qq)]

hist_visibility=defaultdict(int)
for i in weeks_visible:
 hist_visibility[i]=1+hist_visibility[i]

print hist_visibility

# distance analysis per window width per asteroid

TEN=int(len(asteroids)*0.1)
FIFTY=int(len(asteroids)*0.5)
NINETY=int(len(asteroids)*0.9)

for window_width in range(2,40):
 windowed_dist = [0 for i in range(len(asteroids))]
 for anum in range(len(asteroids)):
  dspec = full_dists[anum]
  windowed_dist[anum] = min([max([dspec[t+s] for s in range(window_width)]) for t in range(len(dspec)-window_width)])
 zz = sorted(windowed_dist)
 print window_width, zz[TEN], zz[FIFTY], zz[NINETY]
