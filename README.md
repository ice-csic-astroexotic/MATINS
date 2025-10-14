# MATINS
3D code for magneto-thermal evolution, with realistic microphysics, and coupled temperature-magnetic field evolution.
The magnetic evolution and details about the cubed-sphere grid is described in Dehman et al., 2023, MNRAS, Volume 518, Issue 1, pp.1222-1242. The technical thermal evolution is described in Ascenzi et al. 2024, MNRAS, Volume 533, Issue 1, pp. 201-224.
The code currently uses OPENMP. It works better with Intel compilers. 

# Main developer of MATINS
Clara Dehman,Daniele Viganò, Stefano Ascenzi

# Requirements:
Lapack libraries
# sudo apt-get install liblapacke-dev checkinstall

# To install Lapack library for Mac: 
LAPACK releases are [available on netlib](http://www.netlib.org/lapack/).

# To compile, simply:
make
# To run, prepare the in/input.dat file and type:
./mt3d


# In paraview
# For the temperature:
Representation -> Surface
Coloring -> Temperature

# For the magnetic field vector field
# Right click on the file
Add filter   # In previous releases it can be found in Filters menu
Alphabetical
Glyph
Glyph Type -> Arrow
Orientation Array -> B
Scale Array -> B , o bien, No scale array
Adjust Scale Factor at will
Coloring -> B
Apply
