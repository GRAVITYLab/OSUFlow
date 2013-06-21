== COMPILE with fPIC ==
To be able to compile the plugin as a shared library, MPI and diy should be configured with --with-fpic
Next, set CC="mpicc -fPIC" and CXX="mpicxx -fPIC" as environment variable to use cmake when configuring Paraview and this entire OSUFlow library
