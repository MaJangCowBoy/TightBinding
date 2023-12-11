"""
  Lattice structure to be used in the tight binding model.
  This structure is used to store the Hamiltonian matrix.
  Basic information of the lattice is stored in the fields 
  "sublattice", "spin_orbital", and "lattice_dimension".
  The Hamiltonian matrix is stored in the field "Hamiltonian" of this structure.
"""

struct LatticeSystem
  sublattice::Int64
  spin_orbital::Int64
  lattice_dimension::Array{Int64}
  Hamiltonian::Array{Float64,2}
  Matdim::Int64
end


"""
  Structure could be constructed by the following function.

  sublattice number: N
  sublattice (spin/orbital) degree of freedom: m
  lattice dimension: x * y * z
  Hamiltonian matrix dimension: x * y * z * N * m
"""

function make_lattice(;x::Int64, y::Int64, z::Int64 = 1, N = 1, m = 1)

  LatSize = [x,y,z];
    # system size
  MatDim = x * y * z * N * m;
    # Hamiltonian matrix dimension
  Hamiltonian = zeros(Float64, MatDim, MatDim);
    # Hamiltonian matrix
  return LatticeSystem(N,m, LatSize, Hamiltonian, MatDim);
end