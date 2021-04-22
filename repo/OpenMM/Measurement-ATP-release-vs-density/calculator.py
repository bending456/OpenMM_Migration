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
import time as timer
from scipy.spatial.distance import cdist


###########################################################################
###                           PDB Generator                             ###
###########################################################################
def PDBgenNoPBC(PDBfileName,
                NoCell1,
                NoCell2,
                restingRatio1,
                restingRatio2,
                num_dead_cell,
                num_dead_cell_core):
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
    # determining a total number of resting cells from the overall population
    ## the rest of cells are activated cells
    
    
    restingCell1 = int(NoCell1*restingRatio1)
    restingCell2 = int(NoCell2*restingRatio2)
    TotalNoCell = NoCell1 + NoCell2 + num_dead_cell 
    refnum1 = num_dead_cell + restingCell1
    refnum2 = num_dead_cell + NoCell1
    refnum3 = num_dead_cell + NoCell1 + restingCell2
    
    ## Dead Cell first: There is no distinction between 1st and 2nd compartment for the dead cells 
    i = 0  # this is overall cell count and index for the pdb file
    for i in np.arange(TotalNoCell):
        x = format(randint(0,9),'.3f')
        y = format(randint(0,9),'.3f')
        z = format(randint(0,9),'.3f')
        
        if i <= num_dead_cell:
            name = 'DC'
        elif (i > num_dead_cell and i <= refnum1) or (i > refnum2 and i <=refnum3): 
            name = 'RM'
        else:
            name = 'AM'
               
        if i < 9:
            structure.write("ATOM      "+str(i+1)+"  "+name+"   "+name+"     "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 9 and i < 99:
            structure.write("ATOM     "+str(i+1)+"  "+name+"   "+name+"    "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 99 and i < 999:
            structure.write("ATOM    "+str(i+1)+"  "+name+"   "+name+"   "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 999 and i < 9999:
            structure.write("ATOM   "+str(i+1)+"  "+name+"   "+name+"  "+str(i+1)+"       "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
        elif i >= 9999:
            structure.write("ATOM  "+str(i+1)+"  "+name+"   "+name+"  "+str(i+1)+"      "+str(x)+"   "+str(y)+"   "+str(z)+"  1.00  0.00 \n")
            
    structure.write("ENDMDL")
    structure.close
                   
    return

############################################################################################
###                          Cell Cluster Generator                                      ###
############################################################################################
def cluster(xcenter,
            ycenter,
            xcoord,
            ycoord,
            coord,
            numOfcell):
    
    angle = 137.508
    phi = angle * ( math.pi / 180.0 )
    cspread = 1
    t = numOfcell
   
    # for loops iterate in this case from the first value until < 4, so
    for n in range (0,t):
        r = cspread * math.sqrt(n)
        theta = n * phi

        x = r * math.cos(theta) + xcenter
        y = r * math.sin(theta) + ycenter
    
        xcoord.append(x)
        ycoord.append(y)
        coord.append([x,y,0])
        #marker.append('Dead')
        #factor1.append(float(0.01))
        #factor2.append(float(100))
        
    
    return xcoord, ycoord, coord #, factor1, factor2, marker
        


#######################################################################################
###                           Cell/Particle Distributor                             ###
#######################################################################################
def genCellCoord3D(NoCell1,
                   NoCell2,
                   dead_cell,
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
                   actAuto,
                   shapefactor,
                   NoParticleZone,
                   placement):
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
    coord = []
    marker = []
    autocrine_factor = []
    migration_factor = []
    xo = []
    yo = []
    xo2 = []
    yo2 = []
    coordtoYaml = {}
    minDist = Rlim
    bumper = 1 # this will prevent the placement beyond the actual dimension you assinged
    minDistforDead = 15
    
    ### Cell Population Information 
    #### First Compartment
    numP1 = 0
    maxDeadCore = len(dead_cell)
    NoOfCellsinCluster = dead_cell
    NoOfDeadCells = np.sum(NoOfCellsinCluster)
    maxRestingP1 = int(NoCell1*restingRatio1)
    maxRestingP2 = int(NoCell2*restingRatio2)
    
    maxIter1 = 20000*NoCell1
    loopcounter1 = 1

    
    #### Second Compartment
    numP2 = 0    
    maxIter2 = 20000*NoCell2
    loopcounter2 = 1

    totalNum = 0 
    
    
    ### Configuring dimension where particles are placed
    if shapefactor > 1: # This means it is slab structure
        CentoR = highx/(shapefactor)
        modhighy = highy*(shapefactor/2)
        
        maxX1 = CentoR - bumper
        minX1 = lowx + bumper 
        maxY1 = modhighy - bumper
        minY1 = lowy + bumper
        
        maxY2 = maxY1
        minY2 = minY1 
        maxXdead = CentoR*9
        minXdead = CentoR*4
        
        # Placing cells in the FIRST compartment
        # DeadCells First
        while numP1 < maxDeadCore and loopcounter1 < maxIter1:
            xpossible = randint(minXdead,maxXdead)
            ypossible = randint(minY1,maxY1)
            
            if numP1 == 0:
                numOfCells = NoOfCellsinCluster[numP1]
                xo, yo, coord = cluster(xpossible,ypossible,
                                        xo,yo,coord,
                                        numOfCells)
                numP1 = numP1 + 1
                continue
                
            distance1 = np.sqrt((np.asarray(xo)-xpossible)**2 + (np.asarray(yo)-ypossible)**2)
            if min(distance1) >= minDistforDead:
                numOfCells = NoOfCellsinCluster[numP1]
                xo, yo, coord = cluster(xpossible,ypossible,
                                        xo,yo,coord,
                                        numOfCells)
                numP1 = numP1 + 1
        
            loopcounter1 = loopcounter1 + 1 
            
    
        # Live Cells
        numP1 = 0
        loopcounter1 = 1    
        while numP1 < NoCell1 and loopcounter1 < maxIter1:
            xpossible = randint(minX1,maxX1)
            ypossible = randint(minY2+25,maxY2-25)
            
            distance1 = np.sqrt((np.asarray(xo)-xpossible)**2 + (np.asarray(yo)-ypossible)**2)
            if min(distance1) >= minDist:
                xo.append(xpossible)
                yo.append(ypossible)
                coord.append([xpossible,ypossible,0])
                if numP1 <= maxRestingP1:
                    migration_factor.append(restMig)
                    autocrine_factor.append(restAuto)
                    marker.append('resting')
                elif numP1 > maxRestingP1:
                    migration_factor.append(actMig)
                    autocrine_factor.append(actAuto)
                    marker.append('activated')
                
                numP1 = numP1 + 1
        
            loopcounter1 = loopcounter1 + 1
            
        
        # Placing cells in the SECOND compartment      
        # Live Cells
        if placement != 'gradient':
            if placement =='near':
                maxX2 = CentoR*3 - bumper # in angstrom 
                minX2 = CentoR*2 + bumper # in angstrom
            elif placement == 'far':
                maxX2 = CentoR*9 - bumper # in angstrom 
                minX2 = CentoR*8 + bumper # in angstrom
            elif placement == 'wide':
                maxX2 = CentoR*9 - bumper # in angstrom 
                minX2 = CentoR*3 + bumper # in angstrom
                
            while numP2 < NoCell2 and loopcounter2 < maxIter2:
                xpossible = randint(minX2,maxX2)
                ypossible = randint(minY1,maxY1)
            
                distance2 = np.sqrt((np.asarray(xo2)-xpossible)**2 + (np.asarray(yo2)-ypossible)**2)
                if min(distance1) >= minDist:
                    xo2.append(xpossible)
                    yo2.append(ypossible)
                    coord.append([xpossible,ypossible,0])
                    if numP2 <= maxRestingP2:
                        migration_factor.append(restMig)
                        autocrine_factor.append(restAuto)
                        marker.append('resting')
                    elif numP2 > maxRestingP2:
                        migration_factor.append(actMig)
                        autocrine_factor.append(actAuto)
                        marker.append('activated')
                
                    numP2 = numP2 + 1
        
                loopcounter2 = loopcounter2 + 1
                       
        elif placement == 'gradient':
            maxX2 = CentoR*9 - bumper # in angstrom 
            minX2 = CentoR*8 + bumper # in angstrom
            maxX3 = CentoR*8 - bumper
            minX3 = CentoR*7 + bumper
            maxX4 = CentoR*7 - bumper
            minX4 = CentoR*6 + bumper
            maxX5 = CentoR*6 - bumper
            minX5 = CentoR*5 + bumper
            maxX6 = CentoR*5 - bumper
            minX6 = CentoR*4 + bumper
        
            maxX = [maxX2, maxX3, maxX4, maxX5, maxX6]
            minX = [minX2, minX3, minX4, minX5, minX6]
            num = [100, 80, 60, 40, 20]
            
            for n in np.arange(len(num)):
                loopcounter2 = 0
                numP2 = 0
                numPcontinue = 0
                while numP2 < num[n] and loopcounter2 < maxIter2:
                    xpossible = randint(minX[n],maxX[n])
                    ypossible = randint(minY1,maxY1)
                              
                    distance2 = np.sqrt((np.asarray(xo2)-xpossible)**2 + (np.asarray(yo2)-ypossible)**2)
                    if min(distance1) >= minDist:
                        xo2.append(xpossible)
                        yo2.append(ypossible)
                        coord.append([xpossible,ypossible,0])
                        if numPcontinue <= maxRestingP2:
                            migration_factor.append(restMig)
                            autocrine_factor.append(restAuto)
                            marker.append('resting')
                        elif numPcontinue > maxRestingP2:
                            migration_factor.append(actMig)
                            autocrine_factor.append(actAuto)
                            marker.append('activated')
                
                        numP2 += 1
                        numPcontinue += 1 
        
                    loopcounter2 += 1
        
        
                
    
    ################################################################
    elif shapefactor == 1: # This means it is square structure
        R0to1 = NoParticleZone
        R1to2 = CentoR
        R2to3 = highx
        deadCellZone = R2to3*2/3
        xori = highx/2
        yori = highy/2
        minX1 = xori - R1to2
        maxX1 = xori + R1to2
        minY1 = yori - R1to2
        maxY1 = yori + R1to2
        minX2 = xori - R2to3
        maxX2 = xori + R2to3
        minY2 = yori - R2to3
        maxY2 = yori + R2to3

        # Placing cells in the FIRST compartment
        # Dead Cell First 
        while numP1 < maxDeadCore and loopcounter1 < maxIter1:
            xpossible = randint(minX2,maxX2)
            ypossible = randint(minY2,maxY2)
            
            if numP1 == 0:
                numOfCells = NoOfCellsinCluster[numP1]
                xo, yo, coord = cluster(xpossible,ypossible,
                                        xo,yo,coord,
                                        numOfCells)
                numP1 = numP1 + 1
                continue
                
            distance1 = np.sqrt((np.asarray(xo)-xpossible)**2 + (np.asarray(yo)-ypossible)**2)
            distance2 = np.sqrt((xpossible-xori)**2 + (ypossible-yori)**2)

            #if min(distance1) >= minDist and distance2 >= R0to1 and distance2 <= R1to2:
            if min(distance1) >= minDistforDead and distance2 <= deadCellZone:
                numOfCells = NoOfCellsinCluster[numP1]
                xo, yo, coord = cluster(xpossible,ypossible,
                                        xo,yo,coord,
                                        numOfCells)
                numP1 = numP1 + 1
        
            loopcounter1 = loopcounter1 + 1 
            

            
        # Live Cells
        numP1 = 0
        loopcounter1 = 1    
        while numP1 < NoCell1 and loopcounter1 < maxIter1:
            xpossible = randint(minX1,maxX1)
            ypossible = randint(minY1,maxY1)
            
            distance1 = np.sqrt((np.asarray(xo)-xpossible)**2 + (np.asarray(yo)-ypossible)**2)
            distance2 = np.sqrt((xpossible-xori)**2 + (ypossible-yori)**2)
            if min(distance1) >= minDist and distance2 >= R0to1 and distance2 <= R1to2:
                xo.append(xpossible)
                yo.append(ypossible)
                coord.append([xpossible,ypossible,0])
                if numP1 <= maxRestingP1:
                    migration_factor.append(restMig)
                    autocrine_factor.append(restAuto)
                    marker.append('resting')
                elif numP1 > maxRestingP1:
                    migration_factor.append(actMig)
                    autocrine_factor.append(actAuto)
                    marker.append('activated')
                
                numP1 = numP1 + 1
        
            loopcounter1 = loopcounter1 + 1
            
        # Placing cells in the SECOND compartment
        # Live Cells
        while numP2 < NoCell2 and loopcounter2 < maxIter2:
            xpossible = randint(minX2,maxX2)
            ypossible = randint(minY2,maxY2)
            
            distance1 = np.sqrt((np.asarray(xo)-xpossible)**2 + (np.asarray(yo)-ypossible)**2)
            distance2 = np.sqrt((xpossible-xori)**2 + (ypossible-yori)**2)
            if min(distance1) >= minDist and distance2 >= R1to2 and distance2 <= R2to3 :
                xo.append(xpossible)
                yo.append(ypossible)
                coord.append([xpossible,ypossible,0])
                if numP2 <= maxRestingP2:
                    migration_factor.append(restMig)
                    autocrine_factor.append(restAuto)
                    marker.append('resting')
                elif numP2 > maxRestingP2:
                    migration_factor.append(actMig)
                    autocrine_factor.append(actAuto)
                    marker.append('activated')
                
                numP2 = numP2 + 1
        
            loopcounter2 = loopcounter2 + 1
            
               
    #-----------------------------------------------------------------------------
       
    for i in np.arange(len(coord)):
        for j in np.arange(3):
            coord[i][j] = coord[i][j]/10 ## <--- angstrom to nano meter conversion
    
    ## Store coordinates in yaml file 
    with open('dummy.yml','w') as file:
          document = yaml.dump(coord,file)
    
    ## Constitution of Cell Population 
    #### First Compartment
    Num_Dead = NoOfDeadCells
    Num_Resting1 = maxRestingP1
    Num_Active1 = NoCell1 - Num_Resting1
    
    #### Second Compartment
    Num_Resting2 = maxRestingP2
    Num_Active2 = NoCell2 - Num_Resting2
    
    Cell_Constitution = {'Rest1': Num_Resting1,
                         'Act1':  Num_Active1,
                         'Rest2': Num_Resting2,
                         'Act2':  Num_Active2,
                         'Dead': Num_Dead}

    return coord, marker, migration_factor, autocrine_factor, Cell_Constitution


###############################################################################################
###                        Written by Yehven and implemented by Ben                         ###
###############################################################################################

def calcForce(positions,
              numberOfCells,
              numberOfdead,
              Origin,
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
              ConcByCell,
              odes,
              step_size,
              shape_factor):
    '''
    [Parameter Description]
    positions:               position of cells from pdb file 
    numberOfCells:           total number of cells (num of cells in 1st half + num of cells in 2nd half)
    Origin:                  the origin of chemoattractant source 
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
    shape_factor:               = 1 for square > 1 for slab
    '''
    
    # Do some magic here that will return substrate field-related force acting on each particle
    n = 0
    m = 0
    dummy_coord = np.zeros([numberOfCells,2]) # since it's 2D

    for i in enumerate(positions):
        if n > numberOfdead:
            dummy_coord[m][0] = i[1][0].value_in_unit(nanometers)
            dummy_coord[m][1] = i[1][1].value_in_unit(nanometers)
            #dummy_coord[m][2] = i[1][2].value_in_unit(nanometers)*0
            m += 1
        n += 1
    recordedPositions         = dummy_coord
    
    
    # Done
    DistCelltoOrigin          = conc.DistCelltoOrigin(recordedPositions,
                                                      Origin,
                                                      shape_factor,
                                                      state)
    
    # Done
    DistCelltoCell            = conc.DistCelltoCell(recordedPositions)

    # Done
    ConcByOrigin              = conc.ConcByOrigin(DistCelltoOrigin,
                                                  t,
                                                  Diff,
                                                  oriConc,
                                                  state,
                                                  Origin)

    # Done
    ConcbyCell, HC            = conc.ConcByCell(ConcByOrigin,
                                                ConcByCell,
                                                DistCelltoCell,
                                                t,
                                                Diff,
                                                kd,
                                                auto_factor,
                                                cellConc)
       
    fvX, fvY, fvZ, odesnew    = conc.forceGen(ConcbyCell,
                                              Origin,
                                              shape_factor,
                                              recordedPositions,
                                              t,
                                              Diff,
                                              HC,
                                              searchingRange,
                                              highBC,
                                              DisplacementScaleByConc,
                                              cellConc,
                                              oriConc,
                                              state,
                                              odes,
                                              step_size)    
    
    
    if state == 'steady':
        dummy_origin = np.ones(np.array(np.shape(recordedPositions)))*[Origin[0],Origin[1]]
        DistCelltoOrigin = cdist(recordedPositions,dummy_origin)
        dummy, HC            = conc.ConcByCell(ConcByOrigin,
                                               ConcByCell,
                                               DistCelltoOrigin,
                                               t,
                                               Diff,
                                               kd,
                                               auto_factor,
                                               cellConc)
        ConcByOrigin = dummy - oriConc
    
       
    return fvX, fvY, fvZ, ConcbyCell, odesnew, ConcByOrigin

###############################################################################################
###                        Written by Ben and implemented by Ben                            ###
###############################################################################################

def calcStateVariable(numberOfCells,
                      time,
                      ConcbyCell,
                      stateVarOld):
    
    stateVarNew = stateVarOld + (np.random.normal(1e-4,1e-5)*(1-stateVarOld) - np.random.normal(1e-3,1e-4)*ConcbyCell*stateVarOld)*0.002
    stateDepFactor = 1/(1+(0.25/stateVarNew)**2)
           
    return stateVarNew, stateDepFactor

def calcForceModified(numberOfCells,
                      numberOfdead,
                      fvX,
                      fvY,
                      fvZ,
                      stateDepFactor,
                      mig_factor,
                      stateVar):
    # This is where I can make an adjustment for resting and activated cells
    forces = []
    n = 0
    
    if stateVar == 'off':
        stateDepFactor = np.ones(len(stateDepFactor))
    
    for i in range(numberOfCells):
        
        if i <= numberOfdead:
            force_on_particle = [0,0,0]
        else:
            force_on_particle = [fvX[n]*mig_factor[n]*np.random.normal(stateDepFactor[n],stateDepFactor[n]/10), 
                                 fvY[n]*mig_factor[n]*np.random.normal(stateDepFactor[n],stateDepFactor[n]/10), 
                                 0] #fvZ[i]*mig_factor[i]*0#stateDepFactor[i]*0]
            n += 1
            
        forces.append(force_on_particle)
        
    return forces

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