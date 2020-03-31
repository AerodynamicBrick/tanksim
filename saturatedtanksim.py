#Big thank you to aspirespace.org.uk and Rick Newlands
#http://www.aspirespace.org.uk/downloads/Modelling%20the%20nitrous%20run%20tank%20emptying.pdf

import math
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

def H_g(T):
    Tc=309.57
    Tr=T/Tc
    #specific enthalpy of the vapor
    #eq 4.6 thermophysical properties
    #Aplicable from -90C to 36C
    b1=-200.0
    b2=440.055
    b3=-459.701
    b4=434.081
    b5=-485.338
    return(1000*b1+b2*(1-Tr)**(1/3)+b3*(1-Tr)**(2/3)+b4*(1-Tr)+b5*(1-Tr)**(4/3))
def H_l(T):
    Tc=309.57
    Tr=T/Tc
    #specific enthalpy of the saturated liquid
    #eq 4.4 thermophysical properties
    #Aplicable from -90C to 35C
    b1=-200.0
    b2=116.043
    b3=-917.225
    b4=794.779
    b5=-589.587
    return(1000*b1+b2*(1-Tr)**(1/3)+b3*(1-Tr)**(2/3)+b4*(1-Tr)+b5*(1-Tr)**(4/3))
def H_v(T):
    Tc=309.57
    Tr=T/Tc
    #eq 4.5 thermophysical properties
    #Aplicable from -90C to 35C
    return H_g(Tr) - H_l(Tr)
def C_liquid(T):
    Tc=309.57
    Tr=T/Tc
    #isobaric specific heat capacity of the saturated liquid
    #eq 4.7 thermophysical properties
    #Aplicable from -90C to 30C
    b1=2.49973
    b2=0.023454
    b3=-3.80136
    b4=13.0945
    b5=-14.5180
    return(1000*b1*(1+b2*(1-Tr)**(-1)+b3*(1-Tr)+b4*(1-Tr)**2+b5*(1-Tr)**3))
def rho_liquid(T):
    Tc=309.57
    Tr=T/Tc
    #Density of saturated liquid
    #eq 4.2 thermophysical properties
    #Aplicable from -90C to 36C
    
    #constants
    rho_c=452.0 #critical density
    e=2.718281828459
    b1=1.72328
    b2=-.83950
    b3=.51060
    b4=-.10412
    return(rho_c*(e**(b1*(1-Tr)**(1/3)+b2*(1-Tr)**(2/3)+b3*(1-Tr)+b4*(1-Tr)**(4/3))))
def rho_vapor(T):
    Tc=309.57
    Tr=T/Tc
    #Density of saturated vapor
    #eq 4.3 thermophysical properties
    #Aplicable from -90C to 36C

    #constants
    rho_c=452.0 #critical density
    e=2.718281828459

    b1=-1.00900
    b2=-6.28792
    b3=7.50332
    b4=-7.90463
    b5=0.629427
    return(rho_c*(e**(b1*((1/Tr)-1)**(1/3)+b2*((1/Tr)-1)**(2/3)+b3*((1/Tr)-1)+b4*((1/Tr)-1)**(4/3)+b5*((1/Tr)-1)**(5/3))))
def p_vapor(T):
    Tc=309.57
    Tr=T/Tc
    #vapor pressure
    #eq 4.1 thermophysical properties
    #Aplicable from -90C to 36C
    
    #constants
    p_c=7251.0 #critical pressure
    e=2.718281828459
    
    b1=-6.71893 
    b2=1.35966
    b3=-1.3779
    b4=-4.051
    return(1000*p_c*(e**(b1*(1-Tr)+b2*(1-Tr)**(3/2)+b3*(1-Tr)**(5/2)+b4*(1-Tr)**5)))
def Dyn_visc(T):
    Tc=309.57
    Tr=T/Tc
    # Dynamic Viscosity of the Saturated Liquid
    #eq 4.9 thermophysical properties
    #Aplicable from -90C to 36C
    
    #constants
    rho_c=452.0 #critical density
    e=2.718281828459
    b1=1.6089
    b2=2.0439
    b3=5.24
    b4=0.0293423
    theta=(Tc-b3)/(T-b3)
    return(b4*exp(b1*(theta-1)^(1/3)+b2*(theta-1)^(4/3)))


#Physical constants of n2o
rho_c=452.0 #critical density
p_c=7251.0*1000 #critical pressure

#system features
chamber_pressure=350*6894.7573 #350 psi to pa
K=2 #fudge factor, generally around 2 for n2o. fit this until it looks right.
N=2 #number of injectors
A_injector = 0.00000051886844613 #area per injector in m^2
timestep = .001 #time in between loops in s
Vtank=.016

#initial conditions
N2OTemp=294.15 #assumtion
pressure=p_vapor(N2OTemp)


m_total=12 #initial mass of the liquid in the tank in kg
m_liquidNew=(Vtank-(m_total/rho_vapor(N2OTemp)))/((1/rho_liquid(N2OTemp))-(1/rho_vapor(N2OTemp))) #eq 1.12
m_vapor=m_total-m_liquidNew
m_vapor=0 #starting mass of vapour in the tank.


#program values
m_v=0 #mv is the mass vapourized.. The mv supposedly converges quickly to an actual value. 



verbose=False
##########################################################
##############Let's Begin#################################
##########################################################

if(verbose): print(N2OTemp)
if(verbose): print(rho_liquid(N2OTemp))
if(verbose): print(rho_vapor(N2OTemp))

try:
    fig, ax = ax1.subplots()
    x = np.linspace(1, 309.57, 1000)
    y = rho_vapor(x)
    y2 = rho_liquid(x)
    ax.plot(x,y,lw=2,color='#539caf',alpha=1)
    ax.plot(x,y2,lw=2,color='#9053af',alpha=1)
    ax1.show()
except:
    pass


massvaporarray=[[0,0]]
massoldarray=[[0,0]]
masstotalarray=[[0,0]]
massflowarray=[[0,0]]
pressurearray=[[0,0]]
i=0
timetot=0
warn=False



while True:
    pressure=p_vapor(N2OTemp)
    DeltaP=pressure-chamber_pressure
    if(1.2>pressure/chamber_pressure):
        warn=True
        print("WARNING, low injection pressure: "+str(pressure))
    
    timetot=timetot+timestep
    
    DeltaQ=m_v*H_v(N2OTemp) #heat removed from the liquid n2o during vaporisation
    DeltaT=-1*(DeltaQ/(m_total*C_liquid(N2OTemp)))

    if(verbose): print("N2O Temp Before: "+str(N2OTemp))
    N2OTemp=N2OTemp+DeltaT
    if(verbose): print("Delta T: "+str(DeltaT))
    if(verbose): print("N2O Temp: "+str(N2OTemp)+"\n")

    try:
        mdot_ox=math.sqrt((2*rho_liquid(N2OTemp)*DeltaP)/(K/((N*A_injector)**2))) #eq 1.7
    except ValueError:
        print("WARNING: CHAMBER PRESSURE DROPPED BELOW INJECTION PRESSURE")
        warn=True
        #
    
    if(verbose): print("mdot liquid: " + str(mdot_ox))
    if(verbose): print("rho vapor: "+str(rho_vapor(N2OTemp)))
    if(verbose): print("rho liquid: "+str(rho_liquid(N2OTemp)))
    m_liquidOld=m_total-(mdot_ox*timestep) #eq 1.8
    
    if(verbose): print("M Liquid Old: "+str(m_liquidOld))

    m_liquidNew=(Vtank-(m_total/rho_vapor(N2OTemp)))/((1/rho_liquid(N2OTemp))-(1/rho_vapor(N2OTemp))) #eq 1.12
    if(verbose): print("M Liquid New: "+str(m_liquidNew))

    m_v=m_liquidOld-m_liquidNew #eq 1.13
    if(verbose): print("Mv: "+str(m_v))

    m_vapor=m_vapor+m_v
    
    m_total=m_total-(timestep*mdot_ox)

    if(verbose): print("M total: "+str(m_total))
    
    if(m_liquidNew>m_liquidOld):
        break

    massvaporarray.append((timetot,m_vapor))
    massoldarray.append((timetot,m_liquidOld))
    masstotalarray.append((timetot,m_total))
    massflowarray.append((timetot,mdot_ox))
    pressurearray.append((timetot,pressure))

    i=i+1
    print(i)
    if(i==10000):
        #pass
        break


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

"""for ele in massvaporarray[1:]:
    ax1.plot(ele[0],ele[1], 'b.')"""
"""for ele in massoldarray[1:]:
    ax1.plot(ele[0],ele[1], 'g.')"""
"""
for ele in masstotalarray[1:]:
    ax1.plot(ele[0],ele[1], 'g.')"""

for ele in massflowarray[1:]:
    ax1.plot(ele[0],ele[1], 'g.')
    
for ele in pressurearray[1:]:
    ax2.plot(ele[0],ele[1], 'm.')
x_coordinates = [0, timetot]
y_coordinates_crit = [chamber_pressure*1.2, chamber_pressure*1.2]
y_coordinates2_fail = [chamber_pressure, chamber_pressure]


ax2.plot(x_coordinates, y_coordinates_crit, 'r')
ax2.plot(x_coordinates, y_coordinates2_fail, 'r')

plt.title("mass over time")

plt.xlabel("time (seconds)")
ax1.set_ylabel("mass flow (kg/s)")
ax1.legend()

green_patch=mpatches.Patch(color="green", label="mass flow")
#green_patch=mpatches.Patch(color="green", label="vapor mass")
red_crit_patch=mpatches.Patch(color="red", label="Critical Pressure")
red_fail_patch=mpatches.Patch(color="red", label="Failure Pressure")

ax2.set_ylabel("pressure")
mage_patch=mpatches.Patch(color="magenta", label="pressure")

ax1.legend(handles=[green_patch,red_crit_patch,red_fail_patch,mage_patch])
fig.tight_layout()
plt.show()

print(m_vapor)
if(warn):
    print("Completed with WARNINGS")
else:
    print("Completed Successfully")
    

"""
Liquid nitrous flows out of the tank
causing a drop in the level of liquid
nitrous.
"""

"""
This causes an increase in the
head space of nitrous vapour above the
liquid. The nitrous vapour pressure drops due
to this expansion.
"""

"""
Some of the liquid nitrous then
vapourises to try to raise the vapour
pressure back up.
"""

"""
The energy required to vapourise the
liquid comes from the liquid itself, and so
its temperature drops.
"""

"""
This lower temperature lowers the tank
pressure.
"""

"""
The flowrate of nitrous out of the tank
depends on the difference between the
tank and combustion chamber pressure.
"""

"""
The fuel mass flowrate eroded from the
plastic fuel grain depends on the total
flowrate of fuel plus nitrous oxidiser.
"""
