---
title: concentration.py File
nav_order: 1
parent: What's in Python Files?
---
## Required python package 
```python
#######################################################
##            Concentration Gradient Calculator      ##
##                 Written By Ben Chun               ##
#######################################################

import numpy as np
from random import seed
from random import randint
from scipy import special 
```

### Error function based diffusion 
- In this simulation there are three modes of calculation scheme for diffusion of chemoattractant substance in the background. One of three is error function based on the steady state diffusion in 1D. 
```python
############################
#### Diffusion function ####
def func(r,t,D):
    rt = r/(2*np.sqrt(D*t))
    erf = special.erfc(rt)
    return erf
############################
```

```python
#############################################################################################################################
### Distance of Cell from the location of Damaged Cell
### This will be used to calculate the released substance concentration from the damaged cell at the given cell location 
def DistCelltoOrigin(CellCoord,
                     Origin,
                     SourceOfOrigin,
                     state):
    '''
    [Parameter descriptions]
    CellCoord:      the coordinates of particles 
    Origin:         the origin of chemoattractant source 
    SourceOfOrigin: None = One end of simulation box in x axis, Center = the center of simulation box
    state:          the characteristics of diffusion of chemoattractant: steady, error, or linear 
    '''
    NoOfCell                  = np.shape(CellCoord)[1]
    CelltoOrigin_r            = np.zeros(NoOfCell)
    xOrigin, yOrigin, zOrigin = Origin

    for nthCell in np.arange(NoOfCell):
        
        # Differentiate the scheme of distasnce calculation between cells and the source of external chemoattractant source
        if SourceOfOrigin is None:
            xCell = CellCoord[0][nthCell]
            dist = abs(xCell-xOrigin)
        elif SourceOfOrigin == 'Center':
            xCell = CellCoord[0][nthCell]
            yCell = CellCoord[0][nthCell]
            dist = np.sqrt((xCell-xOrigin)**2 + (yCell-yOrigin)**2)
        
        CelltoOrigin_r[nthCell] = dist.copy()
        
    return CelltoOrigin_r
```

```python
####################################################################################################################
### Distance among cells 
### This will be used to calculate the released substance concentration from individual cell at the given cell location 
def DistCelltoCell(CellCoord):
    
    NoOfCell = np.shape(CellCoord)[1]
    CelltoCell_r = {}
    xcoord = CellCoord[0]
    ycoord = CellCoord[1]

    for nthCell in np.arange(NoOfCell):
        xncoord = CellCoord[0][nthCell]
        yncoord = CellCoord[1][nthCell]
        dist = np.sqrt((xcoord-xncoord)**2+(ycoord-yncoord)**2)

        CelltoCell_r[nthCell] = dist.copy() 
    
    return CelltoCell_r # this is in dictionary form but it doesn't have to be 
```

```python
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
```

```python
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
    ConcByOrigin: the concentration released by the external source 
    ConcByCell:   the concentration released by the cell
    CelltoCell_r: the distance among cells dictionary form 
    t:            time
    D:            diffusion rate 
    kd:           constant associated with Hills coefficient for autocrinic release
    auto_factor:  a list of autocrinic factors based on the state of cells
    cellConc:     the max concentration released by cells (determining the degree of intercellular communication)
    '''

    NoOfCell = len(ConcByOrigin)
    ConcByCell = ConcByOrigin + ConcByCell        
    HillsCoefficient = 1/(1+((kd*auto_factor)/ConcByCell)) ## So far, there is no degree of coopertivity <---- this is where I can make adjustment for resting and activated cells 
    data = np.zeros(NoOfCell)

    for nthCell in np.arange(NoOfCell):
        dist = CelltoCell_r[nthCell] ## array with length of NoOfCell 
        dummydata = HillsCoefficient[nthCell]*func(dist,t,D)*cellConc ## array (from specific cell to all cells)
        data = data + dummydata 
        
    ConcByCell = ConcByCell + data
    HillsCoefficient = 1/(1+((kd*auto_factor)/ConcByCell))
    
    return ConcByCell, HillsCoefficient
```

```python
#######################################################################################################
### Generate the force of migration based on the gradient of chemoattractant along wiht Brownian Motion
def forceGen(ConcByCell,
             Origin,
             SourceOfOrigin,
             CellCoord,
             t,
             D,
             HillsCoefficient,
             searchingRange,
             highBC,
             DisplacementScaleByConc,
             cellConc,
             oriConc,
             state):
    '''
    [Parameter description]
    ConcByCell:                 the concentration released by the cell
    Origin:                     the origin of chemoattractant source
    SourceOfOrigin:             None = One end of simulation box in x axis, Center = the center of simulation box 
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

    NoOfCell = np.shape(CellCoord)[1]
    highBCx, highBCy, highBCz = highBC
    ## x1,y1,z1 denotes the x,y,z coordinate of cell location
    FVectorX = np.zeros(NoOfCell)
    FVectorY = np.zeros(NoOfCell)
    FVectorZ = np.zeros(NoOfCell)
    
    # x
    x0 = CellCoord[0]-searchingRange
    x1 = CellCoord[0]
    x2 = CellCoord[0]+searchingRange

    # y 
    y0 = CellCoord[1]-searchingRange
    y1 = CellCoord[1]
    y2 = CellCoord[1]+searchingRange

    # the distance from the one end where ATP is being released
    if SourceOfOrigin is None: 
        r0 = abs(Origin[0]-x0)
        r1 = abs(Origin[0]-x1)
        r2 = abs(Origin[0]-x2)
        Conc_Origin0 = ConcByOrigin(r0,t,D,oriConc,state,Origin)
        Conc_Origin1 = ConcByOrigin(r1,t,D,oriConc,state,Origin)
        Conc_Origin2 = ConcByOrigin(r2,t,D,oriConc,state,Origin)
        CO01 = Conc_Origin0
        CO21 = Conc_Origin2
        CO10 = Conc_Origin1
        CO12 = Conc_Origin1
        CO11 = Conc_Origin1

    elif SourceOfOrigin == 'Center':
        r01 = np.sqrt((Origin[0]-x0)**2 + (Origin[1]-y1)**2)
        r21 = np.sqrt((Origin[0]-x2)**2 + (Origin[1]-y1)**2)
        r10 = np.sqrt((Origin[0]-x1)**2 + (Origin[1]-y0)**2)
        r12 = np.sqrt((Origin[0]-x0)**2 + (Origin[1]-y2)**2)
        r11 = np.sqrt((Origin[0]-x1)**2 + (Origin[1]-y1)**2)
        CO01 = ConcByOrigin(r01,t,D,oriConc,state,Origin)
        CO21 = ConcByOrigin(r21,t,D,oriConc,state,Origin)
        CO10 = ConcByOrigin(r10,t,D,oriConc,state,Origin)
        CO12 = ConcByOrigin(r12,t,D,oriConc,state,Origin)
        CO11 = ConcByOrigin(r11,t,D,oriConc,state,Origin)    


    xcell = CellCoord[0]
    ycell = CellCoord[1]
    
    xprob = np.random.normal(0.25,1,NoOfCell)
    yprob = np.random.normal(0.25,1,NoOfCell)
    
    xrand = np.random.normal(0,0.5,NoOfCell)
    yrand = np.random.normal(0,0.5,NoOfCell)
    
    for nthCell in np.arange(NoOfCell):
        
        rfromCells01 = np.sqrt((xcell-(xcell[nthCell]-searchingRange))**2 + (ycell-ycell[nthCell])**2)
        rfromCells21 = np.sqrt((xcell-(xcell[nthCell]+searchingRange))**2 + (ycell-ycell[nthCell])**2) 
        rfromCells10 = np.sqrt((xcell-xcell[nthCell])**2 + (ycell-(ycell[nthCell]-searchingRange))**2)
        rfromCells12 = np.sqrt((xcell-xcell[nthCell])**2 + (ycell-(ycell[nthCell]+searchingRange))**2)
        rfromCells11 = np.sqrt((xcell-xcell[nthCell])**2 + (ycell-ycell[nthCell])**2)

        Conc_Cell01 = sum(HillsCoefficient[nthCell]*func(rfromCells01,t,D)*cellConc)
        Conc_Cell21 = sum(HillsCoefficient[nthCell]*func(rfromCells21,t,D)*cellConc)
        Conc_Cell10 = sum(HillsCoefficient[nthCell]*func(rfromCells10,t,D)*cellConc)
        Conc_Cell12 = sum(HillsCoefficient[nthCell]*func(rfromCells12,t,D)*cellConc)
        Conc_Cell11 = sum(HillsCoefficient[nthCell]*func(rfromCells11,t,D)*cellConc)
        
        ConcTotal01 = CO01[nthCell] + Conc_Cell01
        ConcTotal21 = CO21[nthCell] + Conc_Cell21 
        ConcTotal10 = CO10[nthCell] + Conc_Cell10 
        ConcTotal12 = CO12[nthCell] + Conc_Cell12 
        ConcTotal11 = CO11[nthCell] + Conc_Cell11

        ## Calculating the concentration gradient         
        dx = (ConcTotal21 - ConcTotal01)/2
        dy = (ConcTotal12 - ConcTotal10)/2
        dz = 0
        
        # Determining in X direction force
        xp = xprob[nthCell]
        DispX = dx*xp
        
        # Determining in Y direction force
        yp = yprob[nthCell]    
        DispY = dy*yp
        
        FVectorX[nthCell] = DisplacementScaleByConc*(DispX + xrand[nthCell]*0.015)*ConcTotal11
        FVectorY[nthCell] = DisplacementScaleByConc*(DispY + yrand[nthCell]*0.015)*ConcTotal11
        FVectorZ[nthCell] = 0
    
    return FVectorX, FVectorY, FVectorZ
```