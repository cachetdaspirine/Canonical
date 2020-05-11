#!/home/hugo/anaconda3/bin/python3
from System import *
import numpy as np

#state=np.array([np.zeros(5,dtype=int) for _ in range(5)])
#state[1,1]=state[0,1]=state[4,1]=state[2,4]=1
#system=System(State=state)

system=System(Lx=5,Ly=5)
IJ=list()
for _ in range(10):
    system.AddRandParticle()
for _ in range(5):
    IJ.append(system.RemoveRandParticle())
for i in range(5):
    system.AddRandParticle(IJ[i],Radius=np.infty)
system.PrintBinary()

print("--------------------------------")
print("Print each individual cluster")
print("--------------------------------")
system.PlotPerSite()

#system.PlotPerSite()
