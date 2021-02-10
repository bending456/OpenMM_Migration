#####################################################################################
# This is an example usage of the interface to feed forces into OpenMM with Python. #
#####################################################################################

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
from sys import stdout
import numpy as np
from FeedInplugin import FeedInForce
from random import seed
from random import randint
import calculator as calc
import time as timer
import os 
import random

here=os.path.dirname(os.path.abspath(__file__))

'''
[Note]: 11/06/2020 by Ben 
This simulation is designed for 2D space but it can be also implemented in 3D
'''
## This is where parameters are being stored
def simulator(oriConc = 10,                     # the max concentration released from the source of chemoattractant
              cellConc = 1,                     # the max concentration released by cells (determining the degree of intercellular communication)
              Diff = 10,                        # the diffusion rate of chemoattractant 
              kd = 10,                          # constant associated with Hills coefficient for autocrinic release
              DispScale = 100000,               # Displacement Scale ... wrong name I guess 
              searchingRange = 0.1,             # searching range around the cell 
              numOfCells1 = 50,                 # a total number of cells in the 1st section of simulation box 
              numOfCells2 = 50,                 # a total number of cells in the 2nd section of simulation box 
              CentoR = None,                       # this will determine the size of 1st section of simulation box in x axis
              pdbFileName = 'temp',             # generated PDB file name 
              dcdfilename = 'test_output2',     # generated DCD file name 
              lowlimBoxlen =0,                  # in angstrom = 1/10 of nanometer (REMEMBER!!) ??? do I need this? 
              highlimBoxlen =1000,              # in angstrom = 1/10 of nanometer (REMEMBER!!)
              simLength = 10000,                # a total number of the simulation steps 
              dumpSize = 100,                   # dump size in the dcd generation process 
              restingRatio1 = 0.8,              # the proportion of resting cells in the overall population of cells in the 1st section 
              restingRatio2 = 0.8,              # the proportion of resting cells in the overall population of cells in the 2nd section 
              restMig = 0.5,                    # the degree of migratory response of resting cells 
              actMig = 0.1,                     # the degree of migratory response of activated cells 
              restAuto = 10,                    # the degree of autocrinic release of resting cells: the larger, the less insensitive
              actAuto = 0.1,                    # the degree of autocrinic release of activated cells
              shape = 'slab',                   # shape of simulation box: either slab or square 
              DiffState = 'steady',             # the characteristics of diffusion of chemoattractant: steady, error, or linear 
              NoParticleZone = None
              ):
    
    # This will generate the random structure that makes no sence but will be processed again
    ## Note: This will generate dummy pdb file that will pass the test run but it is irrelevant to the actual calculation. 
    calc.PDBgenNoPBC(pdbFileName,numOfCells1,numOfCells2,restingRatio1,restingRatio2) 
    pdb = PDBFile(pdbFileName+'.pdb')        
    
    # Configure Dumping frequency 
    dumpingstep = dumpSize
    dcdReporter = DCDReporter(dcdfilename+'.dcd', dumpingstep)
        
    # Simulation Box Dimension in angstrom
    ## Slab case
    if shape == 'slab':
        shape_factor = 5 # What where this goes to.
        dummy = (np.int(highlimBoxlen*shape_factor))/(10)
        OriginOfExternalSignal = [dummy,0,0] # in nanometer 
        # This indicates the far end of X-axis
    ## Square case
    elif shape == 'square':
        shape_factor = 1
        dummy = (np.int(highlimBoxlen))/(10)
        OriginOfExternalSignal = [dummy/2,dummy/2,0]
        # This indicates the center of simulation box
   
   
    # Defining the dimension of simulation box 
    x = [lowlimBoxlen,highlimBoxlen*shape_factor]
    y = [lowlimBoxlen,highlimBoxlen]
    z = [0,0]
    highBC = np.array([x[1],y[1],z[1]])
    
    # Configure Simulation Length 
    Repeat = simLength
    stepFreq = 1 

    # Parameters determining the dynamics of cells 
    MassOfCell = 200
    RepulsiveScale = 0.05 # this can determine the distance among cells themselves
    ### Simulation parameters
    temperature = 1.0 #0.0    # K   <---- set to 0 to get rid of wiggling for now
    frictionCoeff = 1.0  # 1/ps
    step_size = 0.002    # ps
    
    ##########################################################
    ###              OpenMM Simulation Control             ###
    ##########################################################
    # Create an OpenMM System object
    system = System()

    # Create CustomNonbondedForce. We use it to implement repulsive potential between particles.
    # Only keep repulsive part of L-J potential
    nonbond = CustomNonbondedForce("(sigma/r)^12; sigma=0.5*(sigma1+sigma2)")
    nonbond.addPerParticleParameter("sigma")
    # Here we don't use cutoff for repulsive potential, but adding cutoff might speed things up a bit as interaction
    # strength decays very fast.
    nonbond.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    # Add force to the System
    system.addForce(nonbond)
    
    # FeedInForce is an OpenMM Force object that serves as the interface to feed forces into the simulation. It stores
    # additional forces that we want to introduce internally and has a method
    # .updateForceInContext(OpenMM::Context context, vector< vector<double> > in_forces). in_forces should be a vector of
    # 3d vectors - one 3d vector per particle and each 3d vector should contain X,Y,Z components of the force.
    # updateForceInContext() will copy forces from in_forces into internal variable that will keep them
    # until updateForceInContext() is called again. Every simulation step OpenMM will add forces stored in FeedInForce to
    # the simulation. You don't need to call updateForceInContext() every step if the forces didn't change - OpenMM will
    # use the forces from the last call of updateForceInContext().
    in_force = FeedInForce()
    # Add in_force to the system
    system.addForce(in_force)

    num_particles = len(pdb.getPositions())
    
    ### Make it 2D 
    energy_expression = 'k * (z^2)'
    force = openmm.CustomExternalForce(energy_expression)
    force.addGlobalParameter('k', 999)
    for particle_index in range(num_particles):
        force.addParticle(particle_index, [])
    system.addForce(force)
    ## This will generate more reasonable structure but still needs to take number of particle from the arbitrary structure 
    ## The reason I choose to do so is generating PDB file is little tricky due to its sensitivity to line arrangement 
    min_dist_among_cells = 2.5
    xlow_end = x[0]
    xhigh_end = x[1]
    ylow_end = y[0]
    yhigh_end = y[1]
    '''
    [Outcomes of genCellCoord3D]
    initPosInNm: coordnates of cells in array format (num_particles by [x,y,z]) 
    marker: a list of markers for resting and activated cells (length of num_particles) 
    migration_factor: a list of constants for migration degree for resting and activated (length of num_particles)
    autocrine_factor: a list of constants for autocrinic degree for resting and activated (lengh of num_particles)
    '''
    initPosInNm, marker, mig_factor, auto_factor = calc.genCellCoord3D(numOfCells1,
                                                                       numOfCells2,
                                                                       min_dist_among_cells,
                                                                       CentoR,
                                                                       xlow_end,
                                                                       xhigh_end,
                                                                       ylow_end,
                                                                       yhigh_end,
                                                                       restingRatio1,
                                                                       restingRatio2,
                                                                       restMig,
                                                                       actMig,
                                                                       restAuto,
                                                                       actAuto,
                                                                       shape_factor,
                                                                       NoParticleZone)

    for i in range(num_particles):
        # Populate the system with particles.
        # system.addParticle(mass in atomic mass units)
        system.addParticle(MassOfCell)  # 100.0 is a particle mass
        # Add particles to CustomNonbondedForce. Here we define repulsion strength sigma for each particle. If particles
        # have different value of sigma, combination rule sigma=0.5*(sigma1+sigma2) will be applied as defined in
        # CustomNonbondedForce force expression.
        sigma = RepulsiveScale
        nonbond.addParticle([sigma])

    # Create integrator
    integrator = BrownianIntegrator(temperature, frictionCoeff, step_size)

    # Create platform
    platform = Platform.getPlatformByName('CUDA')

    # Create simulation
    simulation = Simulation(pdb.topology, system, integrator, platform)
    print('Simulation is initiated')
    print('REMARK  Using OpenMM platform %s' % simulation.context.getPlatform().getName())

    # Set initial positions for the particles
    simulation.context.setPositions(initPosInNm)

    # Simulate
    # We can add dcd or other reporters to the simulation to get desired output with
    ## Note: Ben's version 
    simulation.reporters.append(dcdReporter)

    time = 0.001 # <----- can't start at 0; Otherwise, the calculation blows up 
    
    
    ## Set the initial stateVariable: Starting from 0 for now. 
    #odeiter = 5
    stateVariable = np.ones([num_particles])

    # storing initial position 
    rest_cell_no = int(numOfCells1*restingRatio1)+int(numOfCells2*restingRatio2)+2
    act_cell_no = num_particles - rest_cell_no
    resting_start = np.zeros([2,rest_cell_no])
    activated_start = np.zeros([2,act_cell_no])

    positions = simulation.context.getState(getPositions=True).getPositions()

    rest_cnt = 0
    act_cnt = 0
    cell_cnt = 0

    for pos in enumerate(positions):
        if marker[cell_cnt] == 'resting' :
            resting_start[0][rest_cnt] = pos[1][0].value_in_unit(angstroms)
            resting_start[1][rest_cnt] = pos[1][1].value_in_unit(angstroms)
            rest_cnt += 1

        elif marker[cell_cnt] == 'activated':
            activated_start[0][act_cnt] = pos[1][0].value_in_unit(angstroms)
            activated_start[1][act_cnt] = pos[1][1].value_in_unit(angstroms)
            act_cnt += 1

        cell_cnt += 1

    print("|------ calibration initiated -----|")
    for num_iter in range(1, 5000): ### What is this? 
        positions = simulation.context.getState(getPositions=True).getPositions()
        simulation.step(stepFreq)
        state = simulation.context.getState(getEnergy=True, getForces=True)
        time += stepFreq
    
    print("|----- simulation initiated -----|")
    
    time_state = 1e-14
    
    Ix = np.ones(num_particles)
    Iy = np.ones(num_particles)
    DMx = np.zeros(num_particles)
    DMy = np.zeros(num_particles)
    UMx = np.zeros(num_particles)
    UMy = np.zeros(num_particles)
    
    odes = {'Ix': Ix,
            'Iy': Iy,
            'DMx': DMx,
            'DMy': DMy,
            'UMx': UMx, 
            'UMy': UMy}
    
    for num_iter in range(1, Repeat): ### What is this? 
        positions = simulation.context.getState(getPositions=True).getPositions()
        # forces_vec will contain external forces that we want to feed into OpenMM.
        # Get current positions of the particles for concentration field/forces calculation
        ##################################################################################################################
        ###                                    Where All the Magical Things Happen                                     ###
        ##################################################################################################################
        if num_iter == 1:
            ConcByCell = np.zeros(num_particles) + 1e-14

        fvX, fvY, fvZ, ConcbyCell, odes = calc.calcForce(positions,
                                                         num_particles,
                                                         OriginOfExternalSignal,
                                                         oriConc,
                                                         cellConc,
                                                         Diff,
                                                         time_state,
                                                         kd,
                                                         highBC,
                                                         DispScale,
                                                         searchingRange,
                                                         marker,
                                                         auto_factor,
                                                         DiffState,
                                                         ConcByCell,
                                                         odes,
                                                         step_size,
                                                         shape_factor)
        
        stateVariable, stateDepFactor = calc.calcStateVariable(num_particles, 
                                                               time_state, 
                                                               ConcbyCell, 
                                                               stateVariable) 
            
        forces_vec                    = calc.calcForceModified(num_particles, 
                                                               fvX, 
                                                               fvY, 
                                                               fvZ, 
                                                               stateVariable, 
                                                               mig_factor)
        
        #print(forces_vec[1])
        
        # Feed external forces into OpenMM
        in_force.updateForceInContext(simulation.context, forces_vec)
        # Advance simulation for 1 steps
        simulation.step(stepFreq)
        state = simulation.context.getState(getEnergy=True, getForces=True)
        time += stepFreq
        time_state += stepFreq
    
    

#!/usr/bin/env python
import sys
import numpy as np
##################################
#
# Revisions
#       10.08.10 inception
#
##################################
#
# Message printed when program run without arguments
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose:

Usage:
"""
  msg+="   -numOfCells [integer number]:  a number of cells in the box \n" 
  msg+="   -lowlimBoxlen [integer number]: low boundary of simulation box in angstrom \n"
  msg+="   -highlimBoxlen [integer number]: high boundary of simulation box in angstrom \n" 
  msg+="   -t [integer number]: a total step of simulation \n" 
  msg+="   -ds [integer number]: dumpting size \n" 
  msg+="   -name [fileName]: output file name \n" 
  msg+="   -input [fileName]: inputfile name- no need to change it \n" 
  msg+="   -CentoR [float]: radius of area where no microglia present \n"
  msg+="   -oriConc [float]: max concentration of substance from origin \n" 
  msg+="   -cellConc [float]: max concentration of substance from cell \n"
  msg+="   -Diff [float]:diffusion coefficient of substance \n"  
  msg+="   -kd [float]: \n"
  msg+="   -DiffScale [float]: diffusion rate of cells \n" 
  msg+="   -searchingRange [float]: searching range in angstrom \n" 
  msg+="   -restingRatio [float]: ratio of resting cells to activated cells \n"
  msg+="   -restMig [float]: resting cells migration rate factor \n"
  msg+="   -actMig [float]: activated cells migration rate factor \n"
  msg+="   -restAuto [float]: resting cells autocrine factor \n"
  msg+="   -actAuto [float]: activated cells autocrine factor \n"
  msg+="   -shape [string]: slap or square (default=None)\n"
  msg+="   -DiffState [string]: steady or linear or erf (default=None) \n"
  msg+="""


Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)
  
  numOfCells1 = 50
  numOfCells2 = 50
  lowlimBoxlen = 0 # in angstrom
  highlimBoxlen =1000 # in angstrom
  simLength = 10000
  dumpSize = 10
  CentoR = None
  pdbFileName = 'temp'
  dcdfilename = 'test_output'
  oriConc = 1 
  cellConc = 0.1  
  Diff = 20 
  kd = 0.5  
  DispScale = 100000
  searchingRange = 0.1
  restingRatio1 = 0.8
  restingRatio2 = 0.8
  restMig = 10
  actMig = 1
  restAuto = 5
  actAuto = 0.01
  shape = 'slap'
  DiffState = 'steady'
  NoParticleZone = None
    
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if arg=="-numOfCells1":
        numOfCells1 = np.int(sys.argv[i+1])
        
    if arg=="-numOfCells2":
        numOfCells2 = np.int(sys.argv[i+1])
        
    if arg=="-lowlimBoxlen":
        lowlimBoxlen = np.int(sys.argv[i+1])
        
    if arg=="-highlimBoxlen":
        highlimBoxlen = np.int(sys.argv[i+1])
    
    if arg=="-t":
        simLength = np.int(sys.argv[i+1])
        
    if arg=="-ds":
        dumpsize = np.int(sys.argv[i+1])
        
    if arg=="-name":
        dcdfilename = sys.argv[i+1]
    
    if arg=="-input":
        pdbFileName = sys.argv[i+1]
        
    if arg=="-oriConc":
        oriConc = np.float(sys.argv[i+1])
    
    if arg=="-cellConc":
        cellConc = np.float(sys.argv[i+1])

    if arg=="-Diff":
        Diff = np.float(sys.argv[i+1])
    
    if arg=="-kd":
        kd = np.float(sys.argv[i+1])
    
    if arg=='-DispScale':
        DispScale = np.float(sys.argv[i+1])
    
    if arg=="-searchingRange":
        searchingRange = np.float(sys.argv[i+1])
        
    if arg=="-restingRatio1":
        restingRatio1 = np.float(sys.argv[i+1])
        
    if arg=="-restingRatio2":
        restingRatio2 = np.float(sys.argv[i+1])

    if arg=="-CentoR":
        CentoR = np.float(sys.argv[i+1])

    if arg=="-restMig":
        restMig = np.float(sys.argv[i+1])

    if arg=="-actMig":
        actMig = np.float(sys.argv[i+1])

    if arg=="-restAuto":
        restAuto = np.float(sys.argv[i+1])

    if arg=="-actAuto":
        actAuto = np.float(sys.argv[i+1])  
    
    if arg=="-shape":
        shape = sys.argv[i+1]
        
    if arg=="-DiffState":
        DiffState = sys.argv[i+1]
        
    if arg=="-NoParticleZone":
        NoParticleZone = np.float(sys.argv[i+1])

simulator(oriConc = oriConc, # diffusion rate of substance from the origin
          cellConc = cellConc,   # diffusion rate of substance from the cell
          Diff = Diff,
          kd = kd,  
          DispScale = DispScale,
          searchingRange = searchingRange,
          numOfCells1 = numOfCells1,
          numOfCells2 = numOfCells2,
          CentoR = CentoR,
          pdbFileName = pdbFileName,
          dcdfilename = dcdfilename,
          lowlimBoxlen = lowlimBoxlen, # in angstrom
          highlimBoxlen = highlimBoxlen, # in angstrom
          simLength = simLength,
          dumpSize = dumpSize,
          restingRatio1 = restingRatio1,
          restingRatio2 = restingRatio2,
          restMig = restMig,
          actMig = actMig,
          restAuto = restAuto,
          actAuto = actAuto,
          shape = shape,
          DiffState = DiffState,
          NoParticleZone = NoParticleZone
          )
