#######################################################
##            Concentration Gradient Calculator      ##
##                 Written By Ben Chun               ##
#######################################################

import numpy as np
from random import seed
from random import randint
from scipy import special 

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
    NoOfCell                  = np.shape(CellCoord)[1]
    CelltoOrigin_r            = np.zeros(NoOfCell)
    xOrigin, yOrigin, zOrigin = Origin
    xyOrigin = np.array([xOrigin,yOrigin])

    for nthCell in np.arange(NoOfCell):
        
        # Differentiate the scheme of distasnce calculation between cells and the source of external chemoattractant source
        if shape_factor > 1:
            xCell = CellCoord[0][nthCell]
            dist = abs(xCell-xOrigin)
        elif shape_factor == 1:
            xyCell = np.array([CellCoord[0][nthCell],CellCoord[1][nthCell]])
            dist = np.linalg.norm(xyCell-xyOrigin)
        
        CelltoOrigin_r[nthCell] = dist.copy()
        
    return CelltoOrigin_r

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

    for nthCell in np.arange(NoOfCell):
        dist = CelltoCell_r[nthCell] ## array with length of NoOfCell
        #dummydata = HillsCoefficient*func(dist,t,D)*cellConc ## array (from specific cell to all cells)
        dummydata = HillsCoefficient[nthCell]*func(dist,t,D)*cellConc ## array (from specific cell to all cells)
        dummydata[nthCell] = 0.0 # To make sure there is no addition of cell concentration (contained within cell itself) to their local concentration. 
        data = data + dummydata 
        
    ConcByCell = ConcByCell + data # The individual local concentration is the collection of substance released from neighboring cells. 
    HillsCoefficient = 1/(1+((kd*auto_factor)/ConcByCell))
    
    return ConcByCell, HillsCoefficient

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
    if shape_factor > 1:
        #rSource = np.array([abs(Origin[0]-xtest[0]),abs(Origin[0]-xtest[1]),abs(Origin[0]-xtest[2])])
        r0 = abs(Origin[0]-xtest[0])
        r1 = abs(Origin[0]-xtest[1])
        r2 = abs(Origin[0]-xtest[2])
        
        Conc_Origin0 = ConcByOrigin(r0,t,D,oriConc,state,Origin)
        Conc_Origin1 = ConcByOrigin(r1,t,D,oriConc,state,Origin)
        Conc_Origin2 = ConcByOrigin(r2,t,D,oriConc,state,Origin)
        
        COarray = np.array([Conc_Origin0,Conc_Origin2,Conc_Origin1,Conc_Origin1,Conc_Origin1])
        #print(len(COarray))

    elif shape_factor == 1 :
        xSource = np.array([(Origin[0]-xtest[0]),(Origin[0]-xtest[1]),(Origin[0]-xtest[2])])**2
        ySource = np.array([(Origin[1]-ytest[0]),(Origin[1]-ytest[1]),(Origin[1]-ytest[2])])**2
        rSource = np.sqrt(np.array([xSource[0]+ySource[1],xSource[2]+ySource[1],xSource[1]+ySource[0],xSource[0]+ySource[2],xSource[1]+ySource[1]]))
        COarray = ConcByOrigin(rSource,t,D,oriConc,state,Origin)


    xcell = CellCoord[0]
    ycell = CellCoord[1]
    
    #xprob = np.random.normal(0.25,1,NoOfCell)
    #yprob = np.random.normal(0.25,1,NoOfCell)
    
    Conc = np.zeros(NoOfCell)
    deltaX = np.zeros(NoOfCell)
    deltaY = np.zeros(NoOfCell)
    
    probx1 = np.zeros(NoOfCell)
    proby1 = np.zeros(NoOfCell)
    
    for nthCell in np.arange(NoOfCell):
        if cellConc == 0.0:
            CCarray = np.zeros(5)
        
        else:    
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
        
        Conc[nthCell] = ConcTotal[4]
        deltaX[nthCell] = dx #np.random.normal(dx,abs(dx)*10)
        deltaY[nthCell] = dy #np.random.normal(dy,abs(dy)*10)
        
        if np.random.normal(dx,abs(dx)*2) > 0:
            probx1[nthCell] = 1
        elif np.random.normal(dx,abs(dx)*2) <= 0:
            probx1[nthCell] = -1
        
        if np.random.normal(dy,abs(dy)*2) > 0:
            proby1[nthCell] = 1
        elif np.random.normal(dy,abs(dy)*2) <= 0:
            proby1[nthCell] = -1
           
    Ix = odes['Ix'] 
    Iy = odes['Iy'] 
    DMx = odes['DMx'] 
    DMy = odes['DMy'] 
    UMx = odes['UMx'] 
    UMy = odes['UMy'] 
    
    #probx1 = np.random.normal(deltaX,abs(deltaX)*2)
    #proby1 = np.random.normal(deltaY,abs(deltaY)*2)
    if dx <= 1e-14:
        dirFacX1 = 1
    else:
        dirFacX1 = np.random.normal(dx,abs(dx)*2) #/abs(probx1)
    
    if dy <= 1e-14:
        dirFacY1 = 1
    else:
        dirFacY1 = np.random.normal(dy,abs(dy)*2) #/abs(proby1)
    
    xrand = np.random.normal(0,0.5,NoOfCell)
    yrand = np.random.normal(0,0.5,NoOfCell)
    dirFacX2 = xrand/abs(xrand)
    dirFacY2 = yrand/abs(yrand)
    
    
    Inewx = Ix + ((0.01*DMx + 0.01*UMx) - (Conc*0.1 + deltaX*0.5)*Ix)*time_step 
    Inewy = Iy + ((0.01*DMy + 0.01*UMy) - (Conc*0.1 + deltaY*0.5)*Iy)*time_step
    DMnewx = DMx + (deltaX*0.5*Ix - 0.01*DMx)*time_step
    DMnewy = DMy + (deltaY*0.5*Iy - 0.01*DMy)*time_step
    UMnewx = UMx + (Conc*0.1*Ix - (0.01 - 0.005)*UMx)*time_step
    UMnewy = UMy + (Conc*0.1*Iy - (0.01 - 0.005)*UMy)*time_step
    FVectorX = DisplacementScaleByConc*(DMnewx*dirFacX1/abs(dirFacX1) + UMnewx*dirFacX2*0.01)
    FVectorY = DisplacementScaleByConc*(DMnewy*dirFacY1/abs(dirFacY1) + UMnewy*dirFacY2*0.01)
    #print('dirX',dirFacX1,'dirY',dirFacY1)
    
    odesNew = {'Ix': Inewx,
               'Iy': Inewy,
               'DMx': DMnewx,
               'DMy': DMnewy,
               'UMx': UMnewx, 
               'UMy': UMnewy}
    
    return FVectorX, FVectorY, FVectorZ, odesNew
