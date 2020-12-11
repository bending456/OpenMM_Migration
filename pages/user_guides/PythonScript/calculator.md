---
title: calculator.py File
nav_order: 2
parent: What's in Python Files?
---

```python
########################################################
#### Sub functions that need for FeedInForce plugin ####
#### Writtne by Ben chun ###############################
########################################################

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
from sys import stdout
import numpy as np
from FeedInplugin import FeedInForce
from random import seed
from random import randint
import concentration as conc
import random
import scipy
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import yaml 
```

```python
###########################################################################
###                           PDB Generator                             ###
###########################################################################
def PDBgenNoPBC(PDBfileName,
                NoCell1,
                NoCell2,
                restingRatio1,
                restingRatio2):
    '''
    [Note]: 11/06/2020 by Ben 
    In this function, having all resting cells or activated cells will not cause any error.
    
    [Parameter Description]
    PDBfileName:    pdb file name
    NoCell1:        a total number of cells in the 1st half of box 
    NoCell2:        a total number of cells in the 2nd half of box 
    restingRatio1:  the proportion of resting cells in the 1st half of box 
    restingRatio2:  the proportion of resting cells in the 2nd half of box 

    '''    
    # --- Writing pdb file initiated ---
    structure = open(PDBfileName+".pdb","w") 
    structure.write("MODEL     1 \n")

    '''
    [------First half of box-------]
    '''
    xo1 = np.zeros(NoCell1)
    yo1 = np.zeros(NoCell1)
    # although the simulation is in 2D box, it is required to register this in PDB file
    zo1 = np.zeros(NoCell1)    
    # determining a total number of resting cells from the overall population
    ## the rest of cells are activated cells 
    restingCell1 = int(NoCell1*restingRatio1)
    
    for n in np.arange(NoCell1):
        xo1[n] = randint(0,9)
        yo1[n] = randint(0,9)
        zo1[n] = randint(0,9)
    
    for i in np.arange(NoCell1):
        x = format(xo1[i],'.3f')
        y = format(yo1[i],'.3f')
        z = format(zo1[i],'.3f')
        # Registering resting and activated cells with name ID 
        if i <= restingCell1:
            name = 'RM'
        elif i > restingCell1:
            name = 'AM'
            
        if i < 9:
            structure.write("ATOM      "+str(i+1)+"  "+name+"   "+name+"     "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 9 and i < 99:
            structure.write("ATOM     "+str(i+1)+"  "+name+"   "+name+"    "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 99 and i < 999:
            structure.write("ATOM    "+str(i+1)+"  "+name+"   "+name+"   "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 999:
            structure.write("ATOM   "+str(i+1)+"  "+name+"   "+name+"  "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")


    '''
    [------Second half of box-------]
    '''        
    xo2 = np.zeros(NoCell2)
    yo2 = np.zeros(NoCell2)
    # although the simulation is in 2D box, it is required to register this in PDB file
    zo2 = np.zeros(NoCell2)
    # determining a total number of resting cells from the overall population
    ## the rest of cells are activated cells 
    restingCell2 = int(NoCell2*restingRatio2)
    
    for n in np.arange(NoCell2):
        xo2[n] = randint(0,9)
        yo2[n] = randint(0,9)
        zo2[n] = randint(0,9)
        
    for i in np.arange(NoCell2):
        x = format(xo2[i],'.3f')
        y = format(yo2[i],'.3f')
        z = format(zo2[i],'.3f')
        # Registering resting and activated cells with name ID 
        if i <= restingCell2:
            name = 'RM'
        elif i > restingCell2:
            name = 'AM'
        
        num = i + NoCell1
        
        if num < 9:
            structure.write("ATOM      "+str(num+1)+"  "+name+"   "+name+"     "+str(num+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif num >= 9 and num < 99:
            structure.write("ATOM     "+str(num+1)+"  "+name+"   "+name+"    "+str(num+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif num >= 99 and num < 999:
            structure.write("ATOM    "+str(num+1)+"  "+name+"   "+name+"   "+str(num+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif num >= 999:
            structure.write("ATOM   "+str(num+1)+"  "+name+"   "+name+"  "+str(num+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        
    structure.write("ENDMDL")
    structure.close
                   
    return 
```


```python
#######################################################################################
###                           Cell/Particle Distributor                             ###
#######################################################################################
def genCellCoord3D(NoCell1,
                   NoCell2,
                   Rlim,
                   CentoR,
                   lowx,
                   highx,
                   lowy,
                   highy,
                   restingRatio1,
                   restingRatio2,
                   restMig,
                   actMig,
                   restAuto,
                   actAuto):
    '''
    [Parameter Description]
    NoCell1:        a total number of cells in the 1st half of box 
    NoCell2:        a total number of cells in the 2nd half of box 
    restingRatio1:  the proportion of resting cells in the 1st half of box 
    restingRatio2:  the proportion of resting cells in the 2nd half of box 
    Rlim:           the minimum distance among cells 
    CentoR:         this will determine the size of 1st section of simulation box in x axis 
    lowx:           low end of simulation box in x axis
    highx:          high end of simulation box in x axis 
    lowy:           low end of simulation box in y axis
    highy:          high end of simulation box in y axis 
    restMig:        the degree of migratory response of resting cells 
    actMig:         the degree of migratory response of activated cells 
    restAuto:       the degree of autocrinic release of resting cells 
    actAuto:        the degree of autocrinic release of activated cells 
    '''
    coord=[]
    marker=[]
    total_num_cells = NoCell1 + NoCell2
    autocrine_factor = np.zeros(total_num_cells)
    migration_factor = np.zeros(total_num_cells)
    coordtoYaml = {}
    minDist = Rlim
    bumper = 10 # this will prevent the placement beyond the actual dimension you assinged
    
    ## First box ##
    ### Configuring dimension where particles are placed
    # --- X axis ---
    if CentoR > 0:
        if CentoR > highx:
            print("The length of 1st half box should be smaller than the total length of simulation box")
            print("[CentoR < highx], Recommended to be less than half of a total length of simulation box") 
            sys.exit(0)
        elif CentoR < highx and CentoR > bumper:
            maxX1 = CentoR - bumper # in angstrom
        elif CentoR < bumper:
            print("CentoR must be bigger than the bumper size (=10)")
            sys.exit(0)
    elif CentoR == 0:
        CentoR = highx/3
        maxX1 = CentoR-bumper
        
        
    minX1 = lowx + bumper # in angstrom
    # --- Y axis ---
    maxY1 = highy - bumper 
    minY1 = lowy + bumper 

    numP1 = 0
    maxP1 = NoCell1
    maxRestingP1 = int(NoCell1*restingRatio1)
    maxIter1 = 20000*maxP1
    loopcounter1 = 1
    xo = np.zeros(NoCell1)
    yo = np.zeros(NoCell1)
    
    #print(maxP1)
    while numP1 < maxP1 and loopcounter1 < maxIter1:
        xpossible = randint(minX1,maxX1)
        ypossible = randint(minY1,maxY1)
        if numP1 == 0:
            xo[numP1] = xpossible
            yo[numP1] = ypossible
            # -- marking state of cell -- #
            marker.append('resting')
            migration_factor[numP1] = restMig
            autocrine_factor[numP1] = restAuto
            # --------------------------- #
            coord.append([xpossible,ypossible,0])
            coordtoYaml[str(numP1)] = [xpossible/10,ypossible/10,0]
            numP1 = numP1 + 1
            continue

        distance1 = np.sqrt((xo-xpossible)**2 + (yo-ypossible)**2)
        if min(distance1) >= minDist:
            xo[numP1] = xpossible
            yo[numP1] = ypossible
            coord.append([xpossible,ypossible,0])
            coordtoYaml[str(numP1)] = [xpossible/10,ypossible/10,0]
            # -- marking state of cell -- #
            if numP1 <= maxRestingP1:
                marker.append('resting')
                migration_factor[numP1] = restMig
                autocrine_factor[numP1] = restAuto
            elif numP1 > maxRestingP1:
                marker.append('activated')
                migration_factor[numP1] = actMig
                autocrine_factor[numP1] = actAuto
            # --------------------------- #
            numP1 = numP1 + 1
        loopcounter1 = loopcounter1 + 1
    
    #print(numP1)
    ## Second box ##
    ### Configuring dimension where particles are placed
    # --- X axis ---
    maxX2 = CentoR*2 - bumper # in angstrom 
    minX2 = CentoR + bumper # in angstrom
    # --- Y axis ---
    maxY2 = highy - bumper 
    minY2 = lowy + bumper 

    numP2 = 0
    maxP2 = NoCell2
    maxRestingP2 = int(NoCell2*restingRatio2)
    maxIter2 = 20000*NoCell2
    loopcounter2 = 1
    xo2 = np.zeros(NoCell2)
    yo2 = np.zeros(NoCell2)
    #print(maxP2)    
    while numP2 < maxP2 and loopcounter2 < maxIter2:
        xpossible = randint(minX2,maxX2)
        ypossible = randint(minY2,maxY2)
        if numP2 == 0:
            xo2[numP2] = xpossible
            yo2[numP2] = ypossible
            # -- marking state of cell -- #
            marker.append('resting')
            migration_factor[numP1+numP2] = restMig
            autocrine_factor[numP1+numP2] = restAuto
            # --------------------------- #
            coord.append([xpossible,ypossible,0])
            coordtoYaml[str(numP1+numP2)] = [xpossible/10,ypossible/10,0]
            numP2 = numP2 + 1
            continue

        distance1 = np.sqrt((xo2-xpossible)**2 + (yo2-ypossible)**2)
        if min(distance1) >= minDist:
            xo2[numP2] = xpossible
            yo2[numP2] = ypossible
            coord.append([xpossible,ypossible,0])
            coordtoYaml[str(numP1+numP2)] = [xpossible/10,ypossible/10,0]
            # -- marking state of cell -- #
            if numP2 <= maxRestingP2:
                marker.append('resting')
                migration_factor[numP1+numP2] = restMig
                autocrine_factor[numP1+numP2] = restAuto
            elif numP2 > maxRestingP2:
                marker.append('activated')
                migration_factor[numP1+numP2] = actMig
                autocrine_factor[numP1+numP2] = actAuto
            # --------------------------- #
            numP2 = numP2 + 1
        loopcounter2 = loopcounter2 + 1
    
    #print(numP2)
    
    #print(np.shape(coord))
    #coord = coord/10
    for i in np.arange(NoCell1+NoCell2):
        for j in np.arange(3):
            coord[i][j] = coord[i][j]/10 ## <--- angstrom to nano meter conversion
    
    ## Store coordinates in yaml file 
    with open('dummy.yml','w') as file:
          document = yaml.dump(coord,file)
                           
    return coord, marker, migration_factor, autocrine_factor 
```


```python
###############################################################################################
###                        Written by Yehven and implemented by Ben                         ###
###############################################################################################

def calcForce(positions,
              numberOfCells,
              Origin,
              SourceOfOrigin,
              oriConc,
              cellConc,
              Diff,
              t,
              kd,
              highBC,
              DisplacementScaleByConc,
              searchingRange,
              marker,
              auto_factor,
              state,
              ConcByCell):
    '''
    [Parameter Description]
    positions:               position of cells from pdb file 
    numberOfCells:           total number of cells (num of cells in 1st half + num of cells in 2nd half)
    Origin:                  the origin of chemoattractant source 
    SourceOfOrigin:          None = One end of simulation box in x axis, Center = the center of simulation box
    oriConc:                 the max concentration released from the source of chemoattractant
    cellConc:                the max concentration released by cells (determining the degree of intercellular communication)
    Diff:                    the diffusion rate of chemoattractant 
    t:                       time 
    kd:                      constant associated with Hills coefficient for autocrinic release
    highBC:                  [max X, max Y, max Z]
    DisplacementScaleByConc: Displacement Scale 
    searchingRange:          searching range around the cell 
    marker:                  a list of marker for resting and activated cells 
    auto_factor:             a list of autocrinic factors based on the state of cells 
    state:                   the characteristics of diffusion of chemoattractant: steady, error, or linear  
    ConcByCell:              the concentration accumulated by the attractant released from each individual cell
    '''
    
    # Do some magic here that will return substrate field-related force acting on each particle
    n = 0
    dummy_coord = np.zeros([3,numberOfCells])

    #ConcByCell                = np.zeros(numberOfCells) + 1e-14

    for i in enumerate(positions):           
        dummy_coord[0][n] = i[1][0].value_in_unit(nanometers)
        dummy_coord[1][n] = i[1][1].value_in_unit(nanometers)
        dummy_coord[2][n] = i[1][2].value_in_unit(nanometers)*0
        n += 1
        
    recordedPositions         = dummy_coord

    DistCelltoOrigin          = conc.DistCelltoOrigin(recordedPositions,
                                                      Origin,
                                                      SourceOfOrigin,
                                                      state)

    DistCelltoCell            = conc.DistCelltoCell(recordedPositions)

    ConcByOrigin              = conc.ConcByOrigin(DistCelltoOrigin,
                                                  t,
                                                  Diff,
                                                  oriConc,
                                                  state,
                                                  Origin)

    ConcbyCell, HC            = conc.ConcByCell(ConcByOrigin,
                                                ConcByCell,
                                                DistCelltoCell,
                                                t,
                                                Diff,
                                                kd,
                                                auto_factor,
                                                cellConc)

    fvX, fvY, fvZ             = conc.forceGen(ConcbyCell,
                                              Origin,
                                              SourceOfOrigin,
                                              recordedPositions,
                                              t,
                                              Diff,
                                              HC,
                                              searchingRange,
                                              highBC,
                                              DisplacementScaleByConc,
                                              cellConc,
                                              oriConc,
                                              state)    
       
    return fvX, fvY, fvZ, ConcbyCell
```

```python
###############################################################################################
###                        Written by Ben and implemented by Ben                            ###
###############################################################################################

def stateFunc(y,t,c):
    ## Rate of state transition 
    kf = 0.1
    kb = 0.01
    
    dydt = kf*c*(1-y) - kb*y
    
    return dydt

def calcStateVariable(numberOfCells,
                      repeat,
                      timeb,
                      time,
                      ConcbyCell,
                      stateVariable):
    
    t = scipy.linspace(timeb,time,repeat)
    newstateVar = []
    stateDepFactor = []
       
    for n in np.arange(numberOfCells):
        
        stateVar = odeint(stateFunc,stateVariable[n],t,args=(ConcbyCell[n],))
        dummy1 = np.float(stateVar[-1])
        newstateVar.append(dummy1)
        dummy2 = 1/(1+dummy1/0.50)
        stateDepFactor.append(dummy2)
        
    return newstateVar, stateDepFactor

def calcForceModified(numberOfCells,
                      fvX,
                      fvY,
                      fvZ,
                      stateDepFactor,
                      mig_factor):
    # This is where I can make an adjustment for resting and activated cells
    forces = []
    
    for i in range(numberOfCells):
    
        force_on_particle = [fvX[i]*mig_factor[i]*stateDepFactor[i], 
                             fvY[i]*mig_factor[i]*stateDepFactor[i], 
                             fvZ[i]*mig_factor[i]*stateDepFactor[i]*0]
        
        forces.append(force_on_particle)
        
    return forces
```

```python
###########################################################################
###                           Test RUN - OpenMM                         ###
###########################################################################
def testRunner(filename):
    pdb = PDBFile(filename+'.pdb')
    forcefield = ForceField('/home/bending456/Ben-Code/Modeling-Codes/Codes/OpenMM_Tutorial/Particle_in_box/Particle_Ben.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter('output.pdb', 10))
    simulation.reporters.append(StateDataReporter(stdout, 10, step=True, potentialEnergy=True, temperature=True))
    simulation.step(10)
    
    return

```