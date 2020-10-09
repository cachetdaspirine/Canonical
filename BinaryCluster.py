import numpy as np
import sys
import copy

TopologieDownHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieUpHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieDownTriangle = [(1,0),(-1,0),(0,1)]
TopologieUpTriangle = [(1,0),(-1,0),(0,-1)]

class BinaryCluster:
    def __init__(self,Sites,Lx,Ly,ParticleType='Triangle'):
        if ParticleType == 'Triangle':
            self.TopologieUp = TopologieUpTriangle
            self.TopologieDown = TopologieDownTriangle
        elif ParticleType=='Hexagon':
            self.TopologieUp = TopologieUpHex
            self.TopologieDown = TopologieDownHex
        # Keep track of where the sites are located in the real system
        # RealSpaceSites is a list of tuple (i,j) which represent the
        # location of each particle
        self.RealSpaceSites=Sites
        #self.ShiftSites(Lx,Ly)
        # this define the size of the box in which we are inserting this cluster
        #self.Size=max([self.RealSpaceSites.__len__()+2,5])
        self.Size=self.GetMaximumExtension()+2
        # Build an array of 0/1 in a box where there are only the one we gave as Sites
        Building=False
        while not Building:
            try :
                self.BuildArray()
            except IndexError :
                self.Size+=1
            else:
                Building=True
        self.ComputeBoundarySites()
    #------------------------------------------------------------------
    # The following function would eventually be usefull  for periodic
    # Boundary conditions
    #------------------------------------------------------------------
    # This function shift the sites in order to avoir them to cross
    # The boundary of the system
    #def ShiftSites(self,Lx,Ly):
    #    XShifted=[sites[0] for sites in self.RealSpaceSites]
    #    YShifted=[sites[1] for sites in self.RealSpaceSites]
    #    while Lx
    #Build An array of 0/1 in a smaller window
    #------------------------------------------------------------------
    def GetMaximumExtension(self):
        Xs=np.transpose(list(self.RealSpaceSites))[0]
        Ys=np.transpose(list(self.RealSpaceSites))[1]
        Xextension=max(Xs)-min(Xs)+2
        Yextension=max(Ys)-min(Ys)+2
        return max([Xextension,Yextension])
    def PrintBinary(self):
        for j in reversed(range(self.WindowArray.shape[1])):
            for i in range(self.WindowArray.shape[0]):
                print(str(self.WindowArray[i,j])+" ",end='')
            print('\n',end='')
    def BuildArray(self):
        self.ComputeCenter()
        #Build a square window of size NP*NP to be sure that the aggregate
        # fit in it. It is full of 0 for now
        self.WindowArray=np.array([np.zeros(self.Size,dtype=int)
                            for _ in range(self.Size)])
        MidX,MidY=self.Size//2, self.Size//2#middle of my aggregate
        # make sure that the aggregate is in the central sites
        # has the same orientation in the real/window space
        if(self.Xg+self.Yg)%2==1:
            MidY+=1
        self.WindowSpaceSites=set()
        for ij in self.RealSpaceSites:
            self.WindowSpaceSites.add((ij[0]-self.Xg+MidX,ij[1]-self.Yg+MidY))
        for ij in self.WindowSpaceSites:
            try:
                self.WindowArray[ij[0],ij[1]]=1
            except IndexError:
                raise

    # Given a list of indices in the real space domaine of occupied SitesIndices
    # this function compute the center of mass of the object in order to rebuild
    # the aggregate in a smaller window
    def ComputeCenter(self):
        self.Xg,self.Yg=0,0
        for ij in self.RealSpaceSites:
            self.Xg+=float(ij[0])
            self.Yg+=float(ij[1])
        self.Xg=int(self.Xg/len(self.RealSpaceSites)+1)
        self.Yg=int(self.Yg/len(self.RealSpaceSites)+1)
    def ComputeBoundarySites(self):
        #BoundarySites=set()
        self.BuildOccupiedSites()
        self.NBoundary=0
        for ij in self.OccupiedSite:
            Neighs = self.Get_Neighbors(ij,Border=True)
            self.NBoundary+=Neighs.__len__()
            #for neigh in Neighs:
            #    BoundarySites.add(neigh)
        #self.NBoundary=BoundarySites.__len__()
    def Get_Neighbors(self, ij,Occupied=False,Free=False,Border=False):
        if (ij[0]+ij[1])%2==0:
            Res = np.array(self.TopologieDown)+np.array(ij)
        else :
            Res = np.array(self.TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        if not Border:
            Resreg=np.delete(Res,np.argwhere((Res[:,0]>=self.Size) | (Res[:,0]<0) | (Res[:,1]>=self.Size) | (Res[:,1]<0)),0)
        else:
            Resreg=Res
        #Build a numpy array of tuple
        Resbis=np.empty(Resreg.__len__(),dtype=object)
        Resbis[:] = list(zip(Resreg[:,0],Resreg[:,1]))
        if Border:
            Resbis = list(Resbis)
        #check the occupancie or not
        if Border :
            for n in reversed(range(Resbis.__len__())):
                if Resbis[n][0]>=0 and Resbis[n][1]>=0 and Resbis[n][0]<self.Size and Resbis[n][1]<self.Size:
                    if self.WindowArray[Resbis[n]]!=0:
                        del Resbis[n]
        else :
            if Occupied:
                Resbis=Resbis[np.array([self.WindowArray[r]==1 for r in Resbis ])]
            elif Free:
                Resbis = Resbis[np.array([self.WindowArray[r]==0 for r in Resbis])]

        return set(Resbis)
    def BuildOccupiedSites(self):
        self.OccupiedSite=set()
        for i,line in enumerate(self.WindowArray):
            for j,site in enumerate(line):
                if site==1:
                    self.OccupiedSite.add((i,j))
