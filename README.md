# MATINS

**MATINS** is a 3D code for the MAgneto-Thermal evolution of Isolated Neutron Stars, incorporating realistic microphysics and a fully coupled temperature–magnetic field evolution. The code employs a finite-volume scheme discretized on a cubed-sphere coordinate system. The magnetic field formalism and grid implementation are described in Dehman et al. (2023, *MNRAS*, 518, 1222–1242), while the thermal evolution framework is detailed in Ascenzi et al. (2024, *MNRAS*, 533, 201–224). The code is optimized with OpenMP and performs best with Intel compilers.

# Main Developers of the MATINS code

- Clara Dehman  
- Daniele Viganò  
- Stefano Ascenzi  

# Requirements

These instructions explain how to obtain a copy of the project and set it up on your local machine.

## 1. Clone the repository
```commandline
git clone https://github.com/ice-csic-astroexotic/MATINS
```

## 2. LAPACK Installation and Parallel Usage

### Linux
Install LAPACK with:
```commandline
sudo apt-get install liblapacke-dev checkinstall
```

### Mac
Install LAPACK via Homebrew or download from [Netlib](http://www.netlib.org/lapack/):

```commandline
brew install lapack
```

### Parallel Usage 
For parallel usage of lapack, it must be linked with OpenBlas. To check if the linked version uses OpenBlas:<br>
1 - Check lapack library: ldd mt3d

	The result should be something like that:
    ```commandline
	linux-vdso.so.1 (0x00007ffc543bd000)
	liblapack.so.3 => /usr/lib/x86_64-linux-gnu/liblapack.so.3 (0x00007fe470171000)
	libgfortran.so.5 => /usr/lib/x86_64-linux-gnu/libgfortran.so.5 (0x00007fe46fea9000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fe46fd5a000)
	libmvec.so.1 => /lib/x86_64-linux-gnu/libmvec.so.1 (0x00007fe46fd2c000)
	libgomp.so.1 => /usr/lib/x86_64-linux-gnu/libgomp.so.1 (0x00007fe46fcea000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fe46fccf000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fe46fadd000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fe46faba000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fe4722ce000)
	libquadmath.so.0 => /usr/lib/x86_64-linux-gnu/libquadmath.so.0 (0x00007fe46fa70000)
	libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007fe46fa68000)
    ```
	
2 - Follow links of lapack: ls -la /usr/lib/x86_64-linux-gnu/liblapack.so.3

	This is a posible result:

	/usr/lib/x86_64-linux-gnu/liblapack.so.3 -> /etc/alternatives/liblapack.so.3-x86_64-linux-gnu

	So follow doing ls -la to the result until it reaches the end.

	A library like /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3 should be found.

If linked lapack is not the openblas version, one can install openblas and link with -lopenblas instead of -llapack

## 3. Usage

### To compile, simply:
```commandline
make (or make -j8 for faster compilation).
```
### To run, prepare the in/input.dat file and type:
```commandline
./mt3d
```

### OpenBlas
When already compiled and linked with OpenBlas' lapack, the number of openmp processors must be set typing: 
```commandline
export OMP_NUM_THREADS=X
```
Being X the desired number. Beware not to use more processes than real cores. Current machines show double of real cores due to multithreading techniques.

Number of real cores can be get with: 
```commandline
grep -m 1 'cpu cores' /proc/cpuinfo
```
Once OMP_NUM_THREADS is set in the terminal, the executable can be run normally. One can check the usage of processors with top or any other system monitor. When using openmp, the percentage of usage will exceed 100%.

The export of OMP_NUM_THREADS can be set in bashrc to make it persistent in following executions.

### Intel
Using intel as compiler gives an important boost in performance. <br>
Currently, there is a free intel software toolset:

https://www.intel.com/content/www/us/en/developer/articles/news/free-intel-software-developer-tools.html

Base and HPC kits must be installed. To apply intel variables and paths run: 
```commandline
source /opt/intel/oneapi/setvars.sh
```
To activate openmp in mkl libraries, the following environment variable must be set:
```commandline
export MKL_NUM_THREADS=X
```
Intel provides its own version for compiler, lapack and openmpi. 

By default, intel's flag -O3 applies performance to floating point calculations that are considered "unsafe". The variations implied by “unsafe” are usually very tiny; however, their impact on the final result of a longer calculation may be amplified if the algorithm involves cancellations (small differences of large numbers).
To make floating point calculations safe and consistent, both flags "-fp-model precise -fp-model source" must be added to compilation. The calculations will be safer, but slower. The user should think which option fits better in its simulation.
This web resource may be helpful to take the decision:
https://www.intel.com/content/dam/develop/external/us/en/documents/fp-consistency-102511-326704.pdf

## 4. Visualization (in paraview)

### For the magnetic field vector field
Right click on the file <br>
Add filter (in previous releases it can be found in Filters menu) <br>
Alphabetical <br>
Glyph <br>
Glyph Type -> Arrow <br>
Orientation Array -> B <br>
Scale Array -> B , o bien, No scale array <br>
Adjust Scale Factor at will <br>
Coloring -> B <br>
Apply 

### For the temperature:
Representation -> Surface <br>
Coloring -> Temperature

## 5. VTUNE usage for profiling

Vtune is used to do the profiling of memory or time (Hotspots), for Intel compilers. It needs the -g -tracebak options.
You need to launch from the terminal where the code is:

```commandline
vtune-gui
```

and select which analysis you want to do, in configuration. We have tried Memory performance and Hotspots. You get the results when the code is done.
Use then Callers/Callees to see the % of usage of all subroutines.

## 6. STACK SIZE (Segmentation fault)

Using dynamic arrays throw a segmentation fault when compiled with intel.
To solve this problem, either compile using flag -heap-arrays (https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/advanced-optimization-options/heap-arrays.html)
Or set a high limit for stack size in the OS: ulimit -s unlimited BETTER PERFORMANCE

The flag allows a parameter to set a maximum size for arrays to be created on stack.

Performance difference:
```commandline
FLAG: 	6m47
ULIMIT:	5m45
```

## 7. NUMA information

Relevant for memory affinity
```commandline
numactl -H
lscpu
```


