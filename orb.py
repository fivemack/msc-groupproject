#!/usr/bin/python

from math import (tan,exp,log,sqrt,acos,asin,atan,pi,sin,cos)

def deg(theta):
  return theta*180/pi
def rad(theta):
  return theta*pi/180

def phase_integral(alpha, G):
 A=[3.332,1.862]
 B=[0.631,1.218]
 phi=[exp(-A[i]*tan(alpha/2)**B[i]) for i in range(2)]
 return (1-G)*phi[0] + G*phi[1]

# magnitude from absolute-magnitude, two distances and the phase coefficient
# dsun and dobs in AU, theta is angle sun-asteroid-observer
def magnitude(dsun,dobs,theta,H,G=0.15):
 mag = 5 + 5*log(dsun*dobs)-2.5*log(phase_integral(theta,G)) 
 return mag

# six parameters:
#  epoch / JD
#  semi-major axis (au)
#  eccentricity
#  mean anomaly / deg
#  inclination / deg
#  argument of periapsis / deg
#  longitude of ascending node / deg

bennu=[2457600.5,  1.1264,0.20375,  101.7039, 6.0349, 66.2231, 2.0609]
aethra=[2458600.5, 2.6081153,0.3897105,    272.421295, 24.994693, 255.199349, 258.392616]
earth=[2451545.0,  1.0000,0.01673,  100.47,  0.000, 102.93, 0]
mars =[2451545.0,  1.5237,0.09337,  355.43,  1.852, 336.08, 49.71]

GM_earth = 398600 # km^3 s^-2
GM_sun = 132712440000 # km^3 s^-2

def veclen(v):
  return sqrt(sum([j**2 for j in v]))

def vecdot(X,Y):
  if (len(X)!=len(Y)):
    print "X=%o and Y=%o are of different lengths\n" % (X,Y)
    return None
  return sum([X[i]*Y[i] for i in range(len(X))])

def veccross(X,Y):
  if (len(X)!=3 and len(Y)!=3):
    print "Cross-product is of 3D vectors"
    return None
  print X,Y
  return [X[1]*Y[2]-Y[1]*X[2], X[2]*Y[0]-Y[2]*X[0], X[0]*Y[1]-Y[0]*X[1]]

def vecscale(k,V):
  print "Scaling %s by %f" % (V,k)
  return [k*w for w in V]

def vecadd(X,Y):
  if (len(X)!=len(Y)):
    print "Vectors %s and %s are different lengths\n" % (X,Y)
    return None
  print "Adding %s to %s" % (X,Y)
  return [X[i]+Y[i] for i in range(len(X))]

def ecc_from_pa(perigee, apogee):
  return (apogee-perigee)/(apogee+perigee)

def period_from_pa(perigee, apogee, GM=398600):
  e = ecc_from_pa(perigee, apogee)
  h2 = perigee*GM*(1+e)
  h = sqrt(h2)
  T = 2*pi/(GM**2) * (h/sqrt(1-e**2))**3
  return T

def ecc_from_true(e, theta):
  fac = sqrt((1-e)/(1+e))
  return 2*atan(fac*tan(theta/2))

def true_from_ecc(e, theta):
  fac = sqrt((1+e)/(1-e))
  phi = 2*atan(fac*tan(theta/2))
  if (phi<0):
    phi=phi+2*pi
  return phi

def mean_from_ecc(e,M):
  m = M-e*sin(M)
  if (m<0):
    m=m+2*pi
  return m

# algorithm 3.1: E(e,M)
def ecc_from_mean(e,M):
  if (M<pi):
    E=M+e/2
  else:
    E=M-e/2
  while True:
    f=E-e*sin(E)-M
    fdash=1-e*cos(E)
    if (abs(f/fdash) < 1e-8):
      return E
    E=E-f/fdash

perigee=9600.0
apogee=21000.0
e=ecc_from_pa(perigee,apogee)
period = period_from_pa(perigee, apogee)
print "Eccentricity=%f  period=%f" % (e,period)
for theta1 in range(18):
 theta = 20*theta1
 eccles = ecc_from_true(e, rad(theta))
 mean = mean_from_ecc(e, eccles)
 print "Time to %f is %f  %f %f" %(theta, (mean/(2*pi))*period, eccles, mean)

print ecc_from_mean(e, 3.6029)
print deg(true_from_ecc(e,ecc_from_mean(e, 3.6029)))

# algorithm 4.2 in _Orbital Mechanics for Engineering Students_
def statevec_to_orbel(P,V,GM=398600):
  r = veclen(P)
  v = veclen(V)
  v_r = vecdot(P,V)/r
  print "Radial velocity is %f" % v_r
  H = veccross(P,V)
  print "Angular momentum vector is %s" % H
  h = veclen(H)
  print "Angular momentum is %f" % h
  i = deg(acos(H[2]/h))
  print "Inclination is %f" % i
  N = [-H[1],H[0],0] # perpendicular to angular momentum and to z=0
  print "Node vector is %s" % N
  n = veclen(N)
  Omega = deg(acos(N[0]/n))
  if (N[1]<0):
    Omega = 360-Omega
  print "Right ascension of ascending node is %f" % Omega
  print "v=%s" % v
  eccvec = vecadd(vecscale((v*v/GM)-(1/r),P), vecscale(-r*v_r/GM,V))
  print "Eccentricity vector is %s" % eccvec
  e = veclen(eccvec)
  print "Eccentricity is %f" % e
  omega = deg(acos(vecdot(N,eccvec) / (n*e)))
  if (eccvec[2]<0):
    omega = 360-omega
  print "Argument of perigee is %f " % omega
  true_anomaly = deg(acos(vecdot(eccvec, P)/(e*r)))
  if (v_r < 0):
    true_anomaly = 360-true_anomaly
  print "True anomaly is %f " % true_anomaly

  perigee = (h**2/GM)/(1-e)
  apogee = (h**2/GM)/(1+e)
  a = (perigee+apogee)/2
  
  period = 2*pi/sqrt(GM) * a**1.5
  
  return [a,e,i,omega,Omega,true_anomaly,period]

# example 4.3 test
P=[-6045, -3490, 2500]
V=[-3.457, 6.618, 2.533]

print statevec_to_orbel(P,V)

# orbit of Mars
au = 149597870.7 # BY DEFINITION in km
mars_in_seconds = period_from_pa(1.3814*au, 1.6660*au, GM_sun)
print mars_in_seconds/86400


# orbit of an inconvenient NEA
hermes_in_seconds = period_from_pa(0.6226*au, 2.6878*au, GM_sun)
print hermes_in_seconds
print hermes_in_seconds/86400

asclepius_elements=[2456805.5, 1.0224, 0.3570, 194.55, 4.919, 180.30, 255.30]

asclepius_in_seconds = period_from_pa(0.6574*au, 1.3874*au, GM_sun)
print asclepius_in_seconds
print asclepius_in_seconds/86400
