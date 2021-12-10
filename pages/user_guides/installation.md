---
title: Installation of OpenMM and FeedInForce
nav_order: 2
---

## Before you start... 
This specific tutorial is meant to instructe those who are familiar with Linux/Unix environment with Python background. 
Possibly more than my experience and skillset. 

### 0. Install Prerequisites (CUDA, FLEX, BISON, SWIG, CMAKE, Doxygen, OpenCL, CCMAKE...)

#### a. Installing CUDA: [Instruction Link - CUDA](https://medium.com/@exesse/cuda-10-1-installation-on-ubuntu-18-04-lts-d04f89287130)

After install, execute the following in the terminal:

    nvidia-smi

If the following error message is given

    Failed to initialize NVML: Driver/library version mismatch

Then, simply reboot the system 

If it works just fine, then you should expect the following outcome in the terminal:

    Thu Dec  9 16:06:49 2021
    +-----------------------------------------------------------------------------+
    | NVIDIA-SMI 495.29.05    Driver Version: 495.29.05    CUDA Version: 11.5     |
    |-------------------------------+----------------------+----------------------+
    | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
    |                               |                      |               MIG M. |
    |===============================+======================+======================|
    |   0  NVIDIA GeForce ...  On   | 00000000:01:00.0  On |                  N/A |
    | 36%   34C    P0    32W / 130W |    501MiB /  2000MiB |      3%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+
    +-----------------------------------------------------------------------------+
    | Processes:                                                                  |
    |  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
    |        ID   ID                                                   Usage      |
    |=============================================================================|
    |    0   N/A  N/A      1025      G   /usr/lib/xorg/Xorg                148MiB |
    |    0   N/A  N/A      1562      G   /usr/lib/xorg/Xorg                163MiB |
    |    0   N/A  N/A      2029      G   /usr/bin/gnome-shell              154MiB |
    |    0   N/A  N/A      4281      G   ...AAAAAAAAA= --shared-files        2MiB |
    +-----------------------------------------------------------------------------+ 

#### b. Installing FLEX

    sudo apt-get update
    sudo apt-get install flex

#### c. Installing BISON

    sudo apt-get install bison

#### d. Installing Doxygen: [Instruction Link - Doxygen](https://www.doxygen.nl/download.html)

***Note*** recommend to install 1.8.X version rather than 1.9.X [Error Example](https://github.com/openmm/openmm/issues/3317)

The following command should install 1.8.17

    sudo apt install doxygen

#### e. Installing OpenSSL

    sudo apt-get install libssl-dev

#### f. Installing CMAKE: [Instruction Link - CMAKE](https://cmake.org/install/)

#### g. Installing CCMAKE:
    
    sudo apt-get install cmake-curses-gui

#### h. Installing SWIG

    sudo apt install swig

#### i. Installing OpenCL

    sudo apt install ocl-icd-opencl-dev



### 1. Install OpenMM from source 
[Official Link](http://docs.openmm.org/7.0.0/userguide/library.html#compiling-openmm-from-source-code)

Download the package from [Git Link](https://github.com/openmm/openmm.git)

    git clone https://github.com/openmm/openmm.git

[Download the version it's been working with the current implementation](https://github.com/openmm/openmm/archive/refs/tags/7.4.1.tar.gz)

[Compiling Instruction](http://docs.openmm.org/latest/userguide/library/02_compiling.html)

### 2. Install FeedInForce developed by Yehven 

Download the package from [Git Link](https://github.com/YevChern/FeedInForce_openmm)

    git clone https://github.com/YevChern/FeedInForce_openmm

go to the plugin directory 

    ccmake CMakeLists.txt 
    
    cmake CMakeLists.txt

then edit your openMM directory section and cmake install prefix section with */path/to/OpenMM*

For example:

    CMAKE_BUILD_TYPE
    CMAKE_INSTALL_PREFIX             /usr/local/openmm
    CUDA_HOST_COMPILER               /usr/bin/cc
    CUDA_SDK_ROOT_DIR                CUDA_SDK_ROOT_DIR-NOTFOUND
    CUDA_TOOLKIT_ROOT_DIR            /usr/local/cuda-11
    CUDA_USE_STATIC_CUDA_RUNTIME     ON
    CUDA_rt_LIBRARY                  /usr/lib/x86_64-linux-gnu/librt.so
    FeedIn_BUILD_CUDA_LIB            ON
    FeedIn_BUILD_PYTHON_WRAPPERS     ON
    OPENMM_DIR                       /usr/local/openmm
    PYTHON_EXECUTABLE                PYTHON_EXECUTABLE-NOTFOUND
    SWIG_EXECUTABLE                  /usr/bin/swig

*If FeedIn_BUILD_PYTHON_WRAPPERS is not on, then the next step won't work.*


Also, make sure your python version. 
If python2 and 3 are installed in your system, then */usr/bin/python* for *python2* and */usr/bin/python3* for *python3*

or 

    cmake CMakeLists.txt -DOPENMM_DIR='/path/to/OpenMM' -DCMAKE_INSTALL_PREFIX='/path/to/OpenMM' -DPYTHON_EXECUTABLE='/usr/bin/python' or '/usr/bin/python3'

execute *make install* by the following 

    sudo make install 

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
    
The best bet is something went wrong in your *setup.py*. 
In general, it should be automatically corrected after executing ccmake and cmake install. 
To fix this, go to *FeedInForce_openmm/plugin/python* and correct openmm_dir in *setup.py* as 

    openmm_dir = '/path/to/your/openmm'

### 3. Running exmaples

Run example at *FeedInForce_openmm/example* 

    python run.py

or 

    python3 run.py 

if you are running into error for *libFeedInForce.so* isn't available 

    export LD_LIBRARY_PATH='$LD_LIBRARY_PATH:/path/to/OpenMM/lib:/path/to/OpenMM'

or 

Add the following in *.bashrc* and *source*

    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/OpenMM/lib:/path/to/OpenMM

Exit to terminal and 

    source .bashrc