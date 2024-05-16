#!/usr/bin/env python
#
# Model radio propagation for jicamarca
#
import ephem
import numpy as n
#import stuffr
import matplotlib.pyplot as plt
import ppigrf
import iri2016.profile as iri
#from pyglow.pyglow import Point
#from pyglow import coord as gcoord
import datetime
import coord
import scipy.constants as c
import scipy.interpolate as sint

def get_R(ne,k,Bxyz,aspect,omega=49.92e6*2.0*n.pi):
    zz=n.cross(Bxyz,k)
    t=n.cross(k,zz)
    t=t/n.sqrt(n.dot(t,t))
    # tranverse comp
    Bt=t*n.dot(Bxyz,t)
    # longitudinal comp
    Bl=k*n.dot(Bxyz,k)
    
    Bla=n.sqrt(n.dot(Bl,Bl))
    Bta=n.sqrt(n.dot(Bt,Bt))    
#    print("B %1.2g tranverse %1.2g long %1.2g"%(n.sqrt(n.dot(Bxyz,Bxyz))/1e-9,Bta/1e-9,Bla/1e-9))
 #   print("ne %1.2g"%(ne))
    X=ne*(c.e**2.0)/(c.epsilon_0*c.electron_mass*omega**2.0)
    Y_l=c.e*Bla/(c.electron_mass*omega)
    Y_t=c.e*Bta/(c.electron_mass*omega)
    # polarizations 
    rho0=1j*Y_t**2.0/(2.0*Y_l*(1.0-X)) + 1j*n.sqrt(Y_t**4.0/(4.0*(Y_l**2.0)*((1.0-X)**2.0))+ 1.0)
    rho1=1j*Y_t**2.0/(2.0*Y_l*(1.0-X)) - 1j*n.sqrt(Y_t**4.0/(4.0*(Y_l**2.0)*((1.0-X)**2.0))+ 1.0)
    n0=n.sqrt(1.0-X/(1.0-Y_t**2.0/(2.0*(1.0-X)) + n.sqrt(Y_t**4.0/(4.0*(1.0-X)**2.0) + Y_l**2.0)))
    n1=n.sqrt(1.0-X/(1.0-Y_t**2.0/(2.0*(1.0-X)) - n.sqrt(Y_t**4.0/(4.0*(1.0-X)**2.0) + Y_l**2.0)))
  #  print(n0)
   #x print(n1)
  #  print(rho0*rho1)
   # E0_x = 1.0
  #  E0_y = E0_x*rho0
  #  E1_x = 1.0
  #  E1_y = E1_x*rho1
 #   print(E0_y)
#    print(E1_y)
    return(rho0,rho1,n0,n1)


year=2015
month=10
day=22
hour=0
jro = ephem.Observer()
lon=283.125626
lat=-11.951482
alt=520.0
jro.lon="%1.4f"%(lon)
jro.lat="%1.4f"%(lat)
jro.elevation=alt

az=146.4025
el=88.4351
t = ephem.Date('%d/%02d/%02d 00:00:00'%(year,month,day))
m=ephem.Moon(jro)


dl=500.0
pos = coord.geodetic2ecef(lat, lon, alt)
k = coord.azel_ecef(lat, lon, alt, az, el)
alts=[]
ax_rats=[]
n0s=[]
n1s=[]
aspects=[]
nes=[]
pols=[]
omega=49.92e6*2.0*n.pi
frots=[]

lin_frac=[]
# circular TX
z=n.array([1.0,-1.0j],dtype=n.complex64)/n.sqrt(2.0)
zn=n.array([1.0,-1.0j],dtype=n.complex64)/n.sqrt(2.0)

# linear tx
#z=n.array([1.0,0.0],dtype=n.complex64)
#zn=n.array([1.0,0.0],dtype=n.complex64)

dists=[]
d=0.0
coeff=0.7

date=datetime.datetime(year, month, day, hour, 0)
sim = iri.IRI(date, [50,1000,1], lat, lon)
ne=n.copy(sim["ne"])
alt_km=n.copy(sim["alt_km"])
print(alt_km)
alt_km[0]=-1

alt_km[-1]=1e9
nefun=sint.interp1d(alt_km,ne)

for i in range(5000):
    #print("i")
    pos = pos + k*dl
    d+=dl/1e3
    dists.append(d)
    #print(pos)
    #print(z)
    lat_lon_h=coord.ecef2geodetic(pos[0], pos[1], pos[2])
    if lat_lon_h[2]/1e3 > 3000:
        print("done")
        break
    
    Be, Bn, Bu = ppigrf.igrf(lat_lon_h[1],lat_lon_h[0], lat_lon_h[2]/1e3, date)
#    print(Be)
    ne=nefun(lat_lon_h[2]/1e3)*coeff
#    pt = Point(, lat_lon_h[0], lat_lon_h[1], lat_lon_h[2]/1e3)
 #   pt.run_igrf()
  #  pt.run_iri()
   # ne=pt.ne*1e6*coeff
    if ne < 0.0:
        ne=0.0
    
    nes.append(ne)
    Bxyz=coord.enu2ecef(lat_lon_h[0], lat_lon_h[1], lat_lon_h[2], Be[0,0]*1e-9, Bn[0,0]*1e-9, Bu[0,0]*1e-9)
    print(Bxyz)
    aspect=180.0*n.arccos(n.dot(k,Bxyz)/(n.sqrt(n.dot(Bxyz,Bxyz))))/n.pi
    
    rho0,rho1,n0,n1=get_R(ne,k,Bxyz,aspect)
    #print("z before")
    #print(z)

    ax_rats.append(rho0)
    n0s.append(n0)
    n1s.append(n1)
    k0=(omega/c.c)*n0
    k1=(omega/c.c)*n1
    m0=n.array([1.0,rho0])
    m0=m0/n.sqrt(n.dot(m0,n.conj(m0)))
    m1=n.array([1.0,rho1])
    m1=m1/n.sqrt(n.dot(m1,n.conj(m1)))
#    print(n.dot(m0,n.conj(m1)))

    # project polarization to O and X mode
    # components
    z0=m0*n.dot(n.conj(m0),z)
    z1=m1*n.dot(n.conj(m1),z)
    
    # phase delay for each pol
    z0=z0*n.exp(1j*k0*dl)
    z1=z1*n.exp(1j*k1*dl)
    print("n0 %1.10f n1 %1.10f"%(n0,n1))
    z=z0+z1
    mag=n.sqrt(n.dot(z,n.conj(z)))

    z=z/mag
    zn=zn*n.exp(1j*(omega/c.c)*n0*dl)

    I=n.abs(z[0])**2.0+n.abs(z[1])**2.0
    Q=n.abs(z[0])**2.0-n.abs(z[1])**2.0
    U=2.0*n.real(z[0]*n.conj(z[1]))
    V=-2.0*n.imag(z[0]*n.conj(z[1]))
    L=Q+1j*U
    frots.append(0.5*n.angle(L))
    print("I %1.2f Q %1.2f U %1.2f V %1.2f theta %1.2f L %1.2f"%(I,Q,U,V,0.5*n.angle(L),n.abs(L)))
    lin_frac.append(n.abs(L))
    alts.append(lat_lon_h[2]/1e3)
    aspects.append(aspect)
    print("alt %1.2f km k0 %1.4f k1 %1.4f aspect %1.2f lat %1.2f lon %1.2f"%(lat_lon_h[2]/1e3,k0,k1,aspect,lat_lon_h[0],lat_lon_h[1]))

print(z)
for i in range(5000):
    pos = pos - k*dl
    d+=dl/1e3
    dists.append(d)
    
    lat_lon_h=coord.ecef2geodetic(pos[0], pos[1], pos[2])
    if lat_lon_h[2]/1e3 < 0:
        print("done")
        break

    Be, Bn, Bu = ppigrf.igrf(lat_lon_h[1],lat_lon_h[0], lat_lon_h[2]/1e3,date)

    ne=nefun(lat_lon_h[2]/1e3)*coeff
    
#    pt = Point(datetime.datetime(year, month, day, hour, 0), lat_lon_h[0], lat_lon_h[1], lat_lon_h[2]/1e3)
 #   pt.run_igrf()999999999
  #  pt.run_iri()
   # ne=pt.ne*1e6*coeff
    if ne < 0.0:
        ne=0.0
    
    nes.append(ne)

    Bxyz=coord.enu2ecef(lat_lon_h[0], lat_lon_h[1], lat_lon_h[2], Be[0,0]*1e-9, Bn[0,0]*1e-9, Bu[0,0]*1e-9)
    print(Bxyz)

#    Bxyz=coord.enu2ecef(lat_lon_h[0], lat_lon_h[1], lat_lon_h[2], Be, Bn, Bu)
    aspect=180.0*n.arccos(n.dot(k,Bxyz)/(n.sqrt(n.dot(Bxyz,Bxyz))))/n.pi
    
    rho0,rho1,n0,n1=get_R(ne,-k,Bxyz,aspect)
    #print("z before")
    #print(z)

    ax_rats.append(rho0)
    n0s.append(n0)
    n1s.append(n1)
    k0=(omega/c.c)*n0
    k1=(omega/c.c)*n1
    m0=n.array([1.0,rho0])
    m0=m0/n.sqrt(n.dot(m0,n.conj(m0)))
    m1=n.array([1.0,rho1])
    m1=m1/n.sqrt(n.dot(m1,n.conj(m1)))
   # print(m0)
    # project polarization to O and X mode
    z0=m0*n.dot(n.conj(m0),z)
    z1=m1*n.dot(n.conj(m1),z)
    # phase delay for each pol
    z0=z0*n.exp(1j*k0*dl)
    z1=z1*n.exp(1j*k1*dl)
    print("n0 %1.10f n1 %1.10f"%(n0,n1))
   # print("z before")
  #  print(z)
    z=z0+z1
 #   print("z after")
#    print(z)
    mag=n.sqrt(n.dot(z,n.conj(z)))
#    print(mag)
    z=z/mag
    zn=zn*n.exp(1j*(omega/c.c)*n0*dl)
 #   print(z)
    #
    I=n.abs(z[0])**2.0+n.abs(z[1])**2.0
    Q=n.abs(z[0])**2.0-n.abs(z[1])**2.0
    U=2.0*n.real(z[0]*n.conj(z[1]))
    V=-2.0*n.imag(z[0]*n.conj(z[1]))
    L=Q+1j*U
    frots.append(0.5*n.angle(L))
    print("I %1.2f Q %1.2f U %1.2f V %1.2f theta %1.2f L %1.2f"%(I,Q,U,V,0.5*n.angle(L),n.abs(L)))
    lin_frac.append(n.abs(L))
    alts.append(lat_lon_h[2]/1e3)
    aspects.append(aspect)
    print("alt %1.2f km k0 %1.4f k1 %1.4f aspect %1.2f lat %1.2f lon %1.2f"%(lat_lon_h[2]/1e3,k0,k1,aspect,lat_lon_h[0],lat_lon_h[1]))

aspects=n.array(aspects)
nes=n.array(nes)
alts=n.array(alts)
ax_rats=n.array(ax_rats)
n0s=n.array(n0s)
n1s=n.array(n1s)
plt.subplot(221)
plt.plot(dists,nes)
plt.xlabel("Ray path (km)")
plt.ylabel("N_e (1/m^3)")

plt.subplot(222)
plt.plot(dists,frots)
plt.xlabel("Ray path (km)")
plt.ylabel("Faraday rotation (deg)")

plt.subplot(223)
plt.plot(dists,lin_frac)
plt.xlabel("Ray path (km)")
plt.ylabel("Linearly polarized fraction")

plt.subplot(224)
plt.plot(dists,aspects)
plt.xlabel("Ray path (km)")
plt.ylabel("Magnetic field aspect angle (deg)")
plt.show()

plt.plot(frots)
plt.show()

plt.plot(lin_frac)
plt.title("Linearly polarized fraction (circular on transmit)")
plt.show()


plt.plot(aspects,alts,label="Center of disk")
plt.plot(aspects+0.25,alts,color="lightblue",label="Limb")
plt.plot(aspects-0.25,alts,color="lightblue")
plt.title("Aspect angle")
plt.legend()
plt.ylabel("Altitude (km)")
plt.xlabel("Aspect to B (deg)")
plt.savefig("aspect.png")
#plt.show()

plt.clf()
plt.plot(nes,alts,label="IRI N_e")
plt.title("N_e")
plt.legend()
plt.ylabel("Altitude (km)")
plt.xlabel("N_e (1/m^3")
plt.savefig("ne.png")
#plt.show()

