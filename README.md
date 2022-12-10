# 3d-grid
Cubed sphere 3D code
This version has the possibility of using realistic microphysics, and coupled temperature-magnetic field evolution.
It is currently under debugging and testing.

# Requirements:
Lapack libraries
# sudo apt-get install liblapacke-dev checkinstall

# To install Lapack library for Mac: 
LAPACK releases are [available on netlib](http://www.netlib.org/lapack/).

# To compile, simply:
make

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

