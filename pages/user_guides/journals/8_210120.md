---
title: 01/20/2021
nav_order: 8
parent: Work Log
---

#### Today's Date
Jan-20-2021

--------------------------------------------------------------------------------
## Main Task

### 1. What was the major task #1
#### Task Detail 
- Ran into weird error saying cuda platform error, which indicates the version is not compatible with OpenMM 
```
[Wed Jan 20 11:44 AM:bending456@kant]
< /home/shared/OpenMM_Migration/repo/withODE > python -m simtk.testInstallation

OpenMM Version: 7.5
Git Revision: 7bf58d51830e5dc108232c92e0d69a9608854983

There are 3 Platforms available:

1 Reference - Successfully computed forces
2 CPU - Successfully computed forces
3 CUDA - Error computing forces with CUDA platform

CUDA platform error: Error initializing CUDA: CUDA_ERROR_COMPAT_NOT_SUPPORTED_ON_DEVICE (804) at /tmp/openmm/platforms/cuda/src/CudaContext.cpp:138

Median difference in forces between platforms:

Reference vs. CPU: 1.96076e-06

All differences are within tolerance.
```
- It turns out after rebooting it, everything works fine

- the latest task I did and issue I ran into was the migration in diagonal direction with ODE setup. 
 - Done 
 
- New trial 
 - Very low Diffusion coefficient for chemokine diffusion 
  - current: 0.5 with 20,000 steps - 4 mins
  - new: 0.1 with 50,000 steps - expected to be 9 mins? - 10 mins
  
 - Displacement Scale = 300 -> 10

----------------------------------------------------------
[Note]
