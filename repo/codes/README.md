## File Structure
- operator.py
    |
    - runner.py  [input = example value]           [description of input variables)]
        |        oriConc = 10,                     # the max concentration released from the source of chemoattractant
        |        cellConc = 1,                     # the max concentration released by cells (determining the degree of intercellular communication)
        |        Diff = 10,                        # the diffusion rate of chemoattractant 
        |        kd = 10,                          # constant associated with Hills coefficient for autocrinic release
        |        DispScale = 100000,               # Displacement Scale ... wrong name I guess 
        |        searchingRange = 0.1,             # searching range around the cell 
        |        numOfCells1 = 50,                 # a total number of cells in the 1st section of simulation box 
        |        numOfCells2 = 50,                 # a total number of cells in the 2nd section of simulation box 
        |        CentoR = 0,                       # this will determine the size of 1st section of simulation box in x axis
        |        pdbFileName = 'temp',             # generated PDB file name 
        |        dcdfilename = 'test_output2',     # generated DCD file name 
        |        lowlimBoxlen =0,                  # in angstrom = 1/10 of nanometer (REMEMBER!!) ??? do I need this? 
        |        highlimBoxlen =1000,              # in angstrom = 1/10 of nanometer (REMEMBER!!)
        |        SourceOfOrigin = None,            # None = One end of simulation box in x axis, Center = the center of simulation box 
        |        simLength = 10000,                # a total number of the simulation steps 
        |        dumpSize = 100,                   # dump size in the dcd generation process 
        |        restingRatio1 = 0.8,              # the proportion of resting cells in the overall population of cells in the 1st section 
        |        restingRatio2 = 0.8,              # the proportion of resting cells in the overall population of cells in the 2nd section 
        |        restMig = 0.5,                    # the degree of migratory response of resting cells 
        |        actMig = 0.1,                     # the degree of migratory response of activated cells 
        |        restAuto = 10,                    # the degree of autocrinic release of resting cells: the larger, the less insensitive
        |        actAuto = 0.1,                    # the degree of autocrinic release of activated cells
        |        shape = 'slab',                   # shape of simulation box: either slab or square 
        |        DiffState = 'steady'              # the characteristics of diffusion of chemoattractant: steady, error, or linear 
        |
        - calculator.py
            |            PDBgenNoPBC
            |            genCellCoord3D
            |            calcForce
            |            stateFunc
            |            calcStateVariable
            |            calcForceModified
            |
            - concentration.py
                |                DistCelltoOrigin
                |                DistCelltoCell
                |                ConcByOrigin
                |                ConcByCell
                |                forceGen
                
            
            
       