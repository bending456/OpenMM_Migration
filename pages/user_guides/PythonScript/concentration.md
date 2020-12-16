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

### Calculating distance from the source of chemoattractant source
- There are two major types of calculations based on the location of substance sources.  
    - **Center**: the source of chemoattractant is located in the center of simulation box. Therefore, the calculation of distance between the position of cell and the location of the source is based on the radial fashion. 
    - **None**: It is basically slab, the source of which is located at the end of x axis (right side).

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
    xyOrigin = np.array([xOrigin,yOrigin])

    for nthCell in np.arange(NoOfCell):
        
        # Differentiate the scheme of distasnce calculation between cells and the source of external chemoattractant source
        if SourceOfOrigin is None:
            xCell = CellCoord[0][nthCell]
            dist = abs(xCell-xOrigin)
        elif SourceOfOrigin == 'Center':
            xyCell = np.array([CellCoord[0][nthCell],CellCoord[1][nthCell]])
            dist = np.linalg.norm(xyCell-xyOrigin)
        
        CelltoOrigin_r[nthCell] = dist.copy()
        
    return CelltoOrigin_r
```

### Calculating distance from cell to cell 
- The gradient of substances in this calculation is based on *the error function*. Its source releases a certain amount, but the cells contribute the same chemical to the gradient formation (ex. ATP-mediated ATP release.)

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

### Concentration of substance as a function of cell position with respect to the source of substance. 
- There are three states: linear, steady, and error. 
    - **Linear**: a simple linear function (y=ax+b) governs the concentration gradient's estimation from the source of substances to the cell location.
    - **Steady**: the media is fully saturated with substances. 
    - **Error**: an error function is used to estimate the concentration gradient from the source of substances to the cell location. 

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

### Concentration of substance as a function of cell position with respect to other cells. 
- Their input, which is the distance of one cell from others, is a list of distance with respect to other cell positions. 
- The release of a substance from each cell is a function of the concentration of substance present at the cell's position. The local concentration is applied to the Hill's equation that generates the coefficient and controls the degree of release. 
- ***[Things to remember]*** 
```python 
CelltoCell_r[nthCell] = [ 0 1 1 1 1 1 1 1 2 3 4 5 6 7 8 9 ... nth_distance ]
``` 
- This being said, the data generated by usind this distance will contain the series of substance concentration at each individual location by the selected cell. 

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

### Force Generator based on the concentration gradient 
- This specific code will generate vector indicating the degree of force generated by the concentration gradient. 
    - Direction of vector will be against graident, meaning that from low to high concentration. 
    - ***Warning**: we have not determined the degree of migration based on the concentration gradient. Therefore, the unit is arbitrary. 

- Searching/estimating the concentration gradient
    1. Select a cell 
    2. Generate coordinates for cell's surroundings
        - x11: x coordinate for the selected cell 
        - x01: x11 - searching range 
        - x21: x11 + searching range 
        - vice versa for y coordinate. 
        - A total of 5 sets of (x,y) are generated from this process 
    3. Calculate the distance of each coordinate from the preceeding process to the location of the chemoattractant's source.
    4. Calculate the distance of each coordinate from the surroundings to other cells 
    5. Combine the concentration of substance released by other cells and the concentration of substance released from the source. 
    6. Calculate the vector in x and y direction based on the concentration differences. 
    7. Multiply each element of force by random probability with slighly bias medium 
        - The center of normal distribution is 0.25 
        - If there is no gradient then there will be no force (duh... obviously)
    8. To exaggerate the motion of cell with respect to the local concentration, we introduce another random factor obtained another normal distribution whose center is 0.0. 
        - In *Brownian* motion, the temperature should be adjust to control the motion of particle. In this case, due to the technical difficulty, it is not possible to locally reassign the temperature of each individual particle or cell. 
        - Therefore, we are simply applying a new set of force. 

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
    xtest = np.array([CellCoord[0]-searchingRange,CellCoord[0],CellCoord[0]+searchingRange])

    # y
    ytest = np.array([CellCoord[1]-searchingRange,CellCoord[1],CellCoord[1]+searchingRange])

    # the distance from the one end where ATP is being released
    if SourceOfOrigin is None:
        rSource = np.array([abs(Origin[0]-xtest[0]),abs(Origin[0]-xtest[1]),abs(Origin[0]-xtest[2])])
        Conc_Origin = ConcByOrigin(rSource,t,D,oriConc,state,Origin)
        COarray = np.array([Conc_Origin[0],Conc_Origin[2],Conc_Origin[1],Conc_Origin[1],Conc_Origin[1]])


    elif SourceOfOrigin == 'Center':
        xSource = np.array([(Origin[0]-x0),(Origin[0]-x1),(Origin[0]-x2)])**2
        ySource = np.array([(Origin[1]-y0),(Origin[1]-y1),(Origin[1]-y2)])**2
        rSource = np.sqrt(np.array([xSource[0]+ySource[1],xSource[2]+ySource[1],xSource[1]+ySource[0],xSource[0]+ySource[2],xSource[1]+ySource[1]]))
        COarray = ConcByOrigin(rSource,t,D,oriConc,state,Origin)


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
        
        CCarray = np.array([sum(HillsCoefficient*func(rfromCells01,t,D)*cellConc),
                            sum(HillsCoefficient*func(rfromCells21,t,D)*cellConc),
                            sum(HillsCoefficient*func(rfromCells10,t,D)*cellConc),
                            sum(HillsCoefficient*func(rfromCells12,t,D)*cellConc),
                            sum(HillsCoefficient*func(rfromCells11,t,D)*cellConc)])
        
        ConcTotal = np.array([COarray[0][nthCell],COarray[1][nthCell],COarray[2][nthCell],COarray[3][nthCell],COarray[4][nthCell]]) + CCarray

        ## Calculating the concentration gradient         
        dx = (ConcTotal[1] - ConcTotal[0])/2
        dy = (ConcTotal[3] - ConcTotal[2])/2
        dz = 0
        
        # Determining in X direction force
        xp = xprob[nthCell]
        DispX = dx*xp
        
        # Determining in Y direction force
        yp = yprob[nthCell]    
        DispY = dy*yp
        
        FVectorX[nthCell] = DisplacementScaleByConc*(DispX + xrand[nthCell]*0.015)*ConcTotal[4]
        FVectorY[nthCell] = DisplacementScaleByConc*(DispY + yrand[nthCell]*0.015)*ConcTotal[4]
        FVectorZ[nthCell] = 0
    
    return FVectorX, FVectorY, FVectorZ
```