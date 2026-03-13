# Fever.jl

## Useful commands

### generate mesh
gmsh -2 cable2d_axi.geo -format msh2 -o cable2d_axi.msh

### example for adding a new package to the project
julia --project=. -e 'using Pkg; Pkg.add("CSV")'

### run outside repl 
julia --project=. fever.jl cases/joint_LIMES/case_joint.toml

### initialization on new machine
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'
