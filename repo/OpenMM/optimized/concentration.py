#######################################################
##            Concentration Gradient Calculator      ##
##                 Written By Ben Chun               ##
#######################################################

import numpy as np
from random import seed
from random import randint
from scipy import special 
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist

############################
#### Diffusion function ####
def func(r,t,D):
    rt = r/(2*np.sqrt(D*t))
    erf = special.erfc(rt)
    return erf
############################

#############################################################################################################################
### Distance of Cell from the location of Damaged Cell
### This will be used to calculate the released substance concentration from the damaged cell at the given cell location 
def DistCelltoOrigin(CellCoord,
                     Origin,
                     shape_factor,
                     state):
    '''
    [Parameter descriptions]
    CellCoord:      the coordinates of particles 
    Origin:         the origin of chemoattractant source 
    shape_factor:   = 1 for square > 1 for slab
    state:          the characteristics of diffusion of chemoattractant: steady, error, or linear 
    '''
    NoOfCell                  = np.shape(CellCoord)[0]
    CelltoOrigin_r            = np.zeros(NoOfCell)
    xOrigin, yOrigin, zOrigin = Origin
    xyOrigin = np.array([[xOrigin,yOrigin]])
    
    if shape_factor > 1:
        dist = abs(np.transpose(CellCoord)[0] - xOrigin)
        
    elif shape_factor == 1:
        dist = cdist(xyOrigin,CellCoord)
    
    return CelltoOrigin_r # array

####################################################################################################################
### Distance among cells 
### This will be used to calculate the released substance concentration from individual cell at the given cell location 
def DistCelltoCell(CellCoord):
    
    NoOfCell = np.shape(CellCoord)[0]
    CelltoCell_r = cdist(CellCoord,CellCoord)
       
    return CelltoCell_r # this is in dictionary form but it doesn't have to be 

#######################################################################
### Concentration of Chemoattractant released from the external source 
def ConcByOrigin(CelltoOrigin_r,
                 t,
                 D,
                 oriConc,
                 state,
                 Origin):
    '''
    [Parameter descriptions]
    CelltoOrigin_r: the distance of individual cell from the exteranl source  
    t:              time
    D:              Diffusion rate of chemoattractant 
    oriConc:        the max concentration of external chemoattractant constant 
    state:          the characteristic state of diffusion of chemoattractant 
    '''
    
    xOri, yOri, zOri = Origin
    if state == 'linear':
        Conc = (1-CelltoOrigin_r/xOri)*oriConc
    elif state == 'steady':
        Conc = np.ones(len(CelltoOrigin_r))*oriConc
    elif state == 'error':
        Conc = func(CelltoOrigin_r,t,D)*oriConc
    
    return Conc

##########################################################################################
### Concentration of Chemoattractant reased from the indivudal cells
def ConcByCell(ConcByOrigin,
               ConcByCell,
               CelltoCell_r,
               t,
               D,
               kd,
               auto_factor,
               cellConc):
    '''
    [Parameter description]
    ConcByOrigin: the concentration released by the external source (series of local concentration releasd from the impairment and detected by individual cell)
    ConcByCell:   the concentration released by the cell
    CelltoCell_r: the distance among cells dictionary form 
    t:            time
    D:            diffusion rate 
    kd:           constant associated with Hills coefficient for autocrinic release
    auto_factor:  a list of autocrinic factors based on the state of cells
    cellConc:     the max concentration released by cells (determining the degree of intercellular communication)
    '''

    
    
    auto_factor = np.asarray(auto_factor)
    
    NoOfCell = len(ConcByOrigin)
    ConcByCell = ConcByOrigin + ConcByCell        
    HillsCoefficient = 1/(1+((kd*auto_factor)/ConcByCell)) ## So far, there is no degree of coopertivity <---- this is where I can make adjustment for resting and activated cells 
    data = np.zeros(NoOfCell)
    
    data = cellConc*HillsCoefficient*func(CelltoCell_r,t,D).sum(axis=0)
    
    ConcByCell_new = ConcByCell + data # The individual local concentration is the collection of substance released from neighboring cells. 
    HillsCoefficient_new = 1/(1+((kd*auto_factor)/ConcByCell))
    
    return ConcByCell_new, HillsCoefficient_new

#######################################################################################################
### Generate the force of migration based on the gradient of chemoattractant along wiht Brownian Motion
def forceGen(ConcByCell,
             Origin,
             shape_factor,
             CellCoord,
             t,
             D,
             HillsCoefficient,
             searchingRange,
             highBC,
             DisplacementScaleByConc,
             cellConc,
             oriConc,
             state,
             odes,
             time_step):
    '''
    [Parameter description]
    ConcByCell:                 the concentration released by the cell
    Origin:                     the origin of chemoattractant source
    shape_factor:               = 1 for square > 1 for slab
    CellCoord:                  position of cells from pdb file 
    t:                          time
    D:                          diffusion rate 
    HillsCoefficient:           Hills coefficient for autocrinic release 
    searchingRange:             searching range around the cell
    highBC:                     [max X, max Y, max Z]
    DisplacementScaleByConc:    the degree of displacement 
    cellConc:                   the max concentration of chemoattractant released from the cell via autocrine 
    oriConc:                    the max concentration of chemoattractant released from the external source 
    state:                      the characteristics of diffusion of chemoattractant: steady, error, or linear 
    '''       
    
    
    NoOfCell = np.shape(CellCoord)[0]
    highBCx, highBCy, highBCz = highBC

    if state == 'error' or state == 'linear':
        CellCoord = np.transpose(CellCoord)
        # x
        xtest = np.array([CellCoord[0]-searchingRange,CellCoord[0],CellCoord[0]+searchingRange])

        # y
        ytest = np.array([CellCoord[1]-searchingRange,CellCoord[1],CellCoord[1]+searchingRange])
        
        # the distance from the one end where ATP is being released
        if shape_factor > 1:
            r0 = abs(Origin[0]-xtest[0]) 
            r1 = abs(Origin[0]-xtest[1]) 
            r2 = abs(Origin[0]-xtest[2]) 
        
            Conc_Origin0 = ConcByOrigin(r0,t,D,oriConc,state,Origin)
            Conc_Origin1 = ConcByOrigin(r1,t,D,oriConc,state,Origin)
            Conc_Origin2 = ConcByOrigin(r2,t,D,oriConc,state,Origin)
        
            COarray = np.array([Conc_Origin0,  #xmy
                                Conc_Origin2,  #xpy
                                Conc_Origin1,  #xym
                                Conc_Origin1,  #xyp
                                Conc_Origin1]) #xy

        elif shape_factor == 1 :
            xSource = np.array([(Origin[0]-xtest[0]),(Origin[0]-xtest[1]),(Origin[0]-xtest[2])])**2
            ySource = np.array([(Origin[1]-ytest[0]),(Origin[1]-ytest[1]),(Origin[1]-ytest[2])])**2
            rSource = np.sqrt(np.array([xSource[0]+ySource[1],   #xmy
                                        xSource[2]+ySource[1],   #xpy
                                        xSource[1]+ySource[0],   #xym
                                        xSource[0]+ySource[2],   #xyp 
                                        xSource[1]+ySource[1]])) #xy
            #                                                 
            COarray = ConcByOrigin(rSource,t,D,oriConc,state,Origin)


        xcell = CellCoord[0]
        ycell = CellCoord[1]
    
        Conc = np.zeros(NoOfCell)
        deltaX = np.zeros(NoOfCell)
        deltaY = np.zeros(NoOfCell)
    
        dirFacX1 = np.zeros(NoOfCell)
        dirFacY1 = np.zeros(NoOfCell)
        # a 
        if cellConc == 0.0:
            Concxmy = np.zeros(NoOfCell)
            Concxpy = np.zeros(NoOfCell)
            Concxym = np.zeros(NoOfCell)
            Concxyp = np.zeros(NoOfCell)
            Concxy = np.zeros(NoOfCell)
        else:    
            xcell_m = xcell - searchingRange
            xcell_p = xcell + searchingRange
            ycell_m = ycell - searchingRange
            ycell_p = ycell + searchingRange
        
            xmy = np.transpose(np.vstack([xcell_m,ycell]))
            xpy = np.transpose(np.vstack([xcell_p,ycell]))
            xym = np.transpose(np.vstack([xcell,ycell_m]))
            xyp = np.transpose(np.vstack([xcell,ycell_p]))
        
            Coord = np.transpose(CellCoord)
            distxmy = cdist(Coord,xmy)
            distxpy = cdist(Coord,xpy)
            distxym = cdist(Coord,xym)
            distxyp = cdist(Coord,xyp) 
            
            Concxmy = HillsCoefficient*cellConc*func(distxmy,t,D).sum(axis=0) + COarray[0]
            Concxpy = HillsCoefficient*cellConc*func(distxpy,t,D).sum(axis=0) + COarray[1]
            Concxym = HillsCoefficient*cellConc*func(distxym,t,D).sum(axis=0) + COarray[2]
            Concxyp = HillsCoefficient*cellConc*func(distxyp,t,D).sum(axis=0) + COarray[3]
            Conc  = ConcByCell + COarray[4]
            
            dx = (Concxpy - Concxmy)/2
            dy = (Concxyp - Concxym)/2
            
            dirFacX1 = np.random.normal(dx,abs(dx)/10)
            dirFacY1 = np.random.normal(dy,abs(dy)/10)
            
            dirFacX1[dirFacX1 < -1e-14] = -1
            dirFacX1[dirFacX1 > 1e-14] = 1
            dirFacY1[dirFacY1 < -1e-14] = -1
            dirFacY1[dirFacY1 > 1e-14] = 1
        

    elif state == 'steady':
        deltaX = np.zeros(NoOfCell)
        deltaY = np.zeros(NoOfCell)
        Conc = np.ones(NoOfCell)*oriConc
        dirFacX1 = np.ones(NoOfCell) 
        dirFacY1 = np.ones(NoOfCell)
    
    xrand = np.random.normal(0,0.5,NoOfCell)
    yrand = np.random.normal(0,0.5,NoOfCell)
    
    xrand[xrand < -1e-14] = -1
    xrand[xrand > 1e-14] = 1
    yrand[yrand < -1e-14] = -1
    yrand[yrand > 1e-14] = 1
    
    dirFacX2 = xrand
    dirFacY2 = yrand
    
    
    Ix = odes['Ix'] 
    Iy = odes['Iy'] 
    DMx = odes['DMx'] 
    DMy = odes['DMy'] 
    UMx = odes['UMx'] 
    UMy = odes['UMy'] 
    
    Inewx = Ix + ((0.01*DMx + 0.01*UMx) - (Conc*0.05 + deltaX*0.2)*Ix)*time_step 
    Inewy = Iy + ((0.01*DMy + 0.01*UMy) - (Conc*0.05 + deltaY*0.2)*Iy)*time_step
    DMnewx = DMx + (deltaX*0.2*Ix - 0.01*DMx)*time_step
    DMnewy = DMy + (deltaY*0.2*Iy - 0.01*DMy)*time_step
    UMnewx = UMx + (Conc*0.05*Ix - 0.01*UMx)*time_step
    UMnewy = UMy + (Conc*0.05*Iy - 0.01*UMy)*time_step
    FVectorX = DisplacementScaleByConc*(DMnewx*dirFacX1 + UMnewx*dirFacX2*1)
    FVectorY = DisplacementScaleByConc*(DMnewy*dirFacY1 + UMnewy*dirFacY2*1)
    FVectorZ = np.zeros(NoOfCell)
    #print('dirX',dirFacX1,'dirY',dirFacY1)
    
    odesNew = {'Ix': Inewx,
               'Iy': Inewy,
               'DMx': DMnewx,
               'DMy': DMnewy,
               'UMx': UMnewx, 
               'UMy': UMnewy}
    
    return FVectorX, FVectorY, FVectorZ, odesNew
