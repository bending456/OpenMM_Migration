0. Install CUDA (https://medium.com/@exesse/cuda-10-1-installation-on-ubuntu-18-04-lts-d04f89287130)
1. Install OpenMM from source (http://docs.openmm.org/7.0.0/userguide/library.html#compiling-openmm-from-source-code)

Download the package from https://github.com/openmm/openmm.git

    git clone https://github.com/openmm/openmm.git

2. Install FeedInForce developed by Yehven 

Download the package from https://github.com/YevChern/FeedInForce_openmm

    git clone https://github.com/YevChern/FeedInForce_openmm

go to the plugin directory 

    ccmake CMakeLists.txt 
    
    cmake CMakeLists.txt

then edit your openMM directory section and cmake install prefix section with /path/to/OpenMM
Also, make sure your python version. 
If python2 and 3 are installed in your system, then /usr/bin/python for python2 and /usr/bin/python3 for python3

or 

    cmake CMakeLists.txt -DOPENMM_DIR='/path/to/OpenMM' -DCMAKE_INSTALL_PREFIX='/path/to/OpenMM' -DPYTHON_EXECUTABLE='/usr/bin/python' or '/usr/bin/python3'

execute make install

    make install 

if you are using python installed over system 

    sudo make PythonInstall 

if not 

    make PythonInstall 
    
If you run into error by saying the following:

    compilation terminated.
    error: command 'x86_64-linux-gnu-gcc' failed with exit status 1
    python/CMakeFiles/PythonInstall.dir/build.make:66: recipe for target 'PythonInstall' failed
    make[3]: *** [PythonInstall] Error 1
    CMakeFiles/Makefile2:348: recipe for target 'python/CMakeFiles/PythonInstall.dir/all' failed
    make[2]: *** [python/CMakeFiles/PythonInstall.dir/all] Error 2
    CMakeFiles/Makefile2:355: recipe for target 'python/CMakeFiles/PythonInstall.dir/rule' failed
    make[1]: *** [python/CMakeFiles/PythonInstall.dir/rule] Error 2
    Makefile:240: recipe for target 'PythonInstall' failed
    make: *** [PythonInstall] Error 2
    
The best bet is something went wrong in your setup.py. In general, it should be automatically corrected after executing ccmake and cmake install. 
To fix this, go to FeedInForce_openmm/plugin/python and correct openmm_dir in setup.py as 

    openmm_dir = '/path/to/your/openmm'

3. Running exmaples

Run example at FeedInForce_openmm/example 

    python run.py

or 

    python3 run.py 

if you are running into error for libFeedInForce.so isn't available 

    export LD_LIBRARY_PATH='$LD_LIBRARY_PATH:/path/to/OpenMM/lib:/path/to/OpenMM'

or 

Add the following in .bashrc and source

    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/OpenMM/lib:/path/to/OpenMM

Exit to terminal and 

    source .bashrc