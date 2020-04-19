# bad N-body integrator

from math import sqrt
import copy

G=6.673e-11
AU=1.49597e11
masses=[1.989e30, 5.972e24]
massy=[0,0,0, 0,-(29800*masses[1]/masses[0]),0,  AU,0,0, 0,29800,0]
DRO=1.2
LEO=6371000+800000 # low Earth orbit is above Earth's surface not Earth's core
zeromass=[DRO*AU,0,0,0,24320,0,  AU,LEO,0, 7455,29800,0]

universe=[masses, massy, zeromass]

dt=60.0

def derivs(masses, massy, zeromass):
  n=len(massy)/6
  v=[0 for i in range(6*n)]
  for i in range(n):
    # x_i dot is xdot_i
    v[6*i]=massy[6*i+3]
    v[6*i+1]=massy[6*i+4]
    v[6*i+2]=massy[6*i+5]
    # xdot_i dot is GMr/r^3
    ddx=0; ddy=0; ddz=0
    for j in range(n):
     if (i!=j):
      rx=massy[6*j]-massy[6*i]
      ry=massy[6*j+1]-massy[6*i+1]
      rz=massy[6*j+2]-massy[6*i+2]
      r=sqrt(rx*rx+ry*ry+rz*rz)
      ddx=ddx+G*masses[j]*rx/(r**3)
      ddy=ddy+G*masses[j]*ry/(r**3)
      ddz=ddz+G*masses[j]*rz/(r**3)
    v[6*i+3]=ddx
    v[6*i+4]=ddy
    v[6*i+5]=ddz
  nz=len(zeromass)/6
  vz=[0 for i in range(6*nz)]
  for i in range(nz):
    for ee in [0,1,2]:
      vz[6*i+ee]=zeromass[6*i+ee+3]
    ddx=0; ddy=0; ddz=0;
    for j in range(n):
      rx=massy[6*j]-zeromass[6*i]
      ry=massy[6*j+1]-zeromass[6*i+1]
      rz=massy[6*j+2]-zeromass[6*i+2]
      r=sqrt(rx*rx+ry*ry+rz*rz)
      ddx=ddx+G*masses[j]*rx/(r**3)
      ddy=ddy+G*masses[j]*ry/(r**3)
      ddz=ddz+G*masses[j]*rz/(r**3)
    vz[6*i+3]=ddx
    vz[6*i+4]=ddy
    vz[6*i+5]=ddz
   
  return [v,vz]

f=open("pingle","w")
T = 0
for q in range(1000000):
 olduniverse=copy.deepcopy(universe)
 # distance Earth->Sun
 ug = sqrt((universe[1][6]-universe[1][0])**2+(universe[1][7]-universe[1][1])**2+(universe[1][8]-universe[1][2])**2)
 # two velocities
 v1 = sqrt(sum(t**2 for t in universe[1][3:5]))
 v2 = sqrt(sum(t**2 for t in universe[1][9:11]))
 towrite=[T,ug,v1,v2,universe[1][9],universe[1][10],universe[1][6],universe[1][7],universe[2][0],universe[2][1],universe[2][6],universe[2][7]]
 towrite=[T,sqrt((universe[2][6]-universe[1][6])**2+(universe[2][7]-universe[1][7])**2)]
 for u in towrite:
  f.write(str(u))
  f.write(" ")
 f.write("\n")

 # mid-point method: take half a step then another half-step
 midpoint=copy.deepcopy(universe)
 Ds = derivs(universe[0],universe[1],universe[2])
 D=Ds[0]; Dz=Ds[1]
 for i in range(len(D)):
  midpoint[1][i] = midpoint[1][i] + dt/2 * D[i]
 for i in range(len(Dz)):
  midpoint[2][i] = universe[2][i] + dt/2 * Dz[i]
 Ds2 = derivs(midpoint[0],midpoint[1],midpoint[2])
 D=Ds2[0]; Dz=Ds2[1]
 for i in range(len(D)):
  universe[1][i] = universe[1][i] + dt * D[i]
 for i in range(len(Dz)):
  universe[2][i] = universe[2][i] + dt * Dz[i]
 
 if (universe[1][3]>0 and olduniverse[1][3]<0):
  print "X velocity of planet changes sign at ",T
 if (universe[2][3]>0 and olduniverse[2][3]<0):
  print "X velocity of satellite changes sign at ",T
 T=T+dt
