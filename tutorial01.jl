include("Define_Lattice.jl");
include("Define_Lattice_indexing.jl");
include("Define_Hamiltonian.jl");

xLen = 10; 
yLen = 10;
sublattice = 2;
intFreedom = 2;


mylat = 
  make_lattice( xLen, yLen; N = subLattice, m = intFreedom)