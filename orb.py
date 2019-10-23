#!/usr/bin/python

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

def len(v):
  return sqrt(sum([j**2 for j in v]))

def dot(X,Y):
  if (length(X)!=length(Y)):
    print "X=%o and Y=%o are of different lengths\n"
    return -6
  return sum([X[i]*Y[i] for i in range(length(X))])

# algorithm 4.2 in _Orbital Mechanics for Engineering Students_
def statevec_to_orbel(P,V,GM=398600):
  r = len(P)
  v = len(V)
  v_r = dot(P,V)/r
  print "Radial velocity is v_r"



# example 4.3 test
P=[-6045, -3490, 2500]
V=[-3.457, 6.618 2.533]

print statevec_to_orbel(P,V)
