"""
  This function returns the position of the lattice site.
  The position is represented by a 4-element array. (2D lattice case is only implemented.)
  The first element is the x-coordinate, the second is the y-coordinate,
  the third is the sublattice number, and the fourth is the sublattice degree of freedom.
"""

function pos2index(mylat::LatticeSystem, pos::Array{Int64})
  
  if length(pos) != 5
    error("pos must be a 5-element array");
  end

  x = mylat.lattice_dimension[1];
  y = mylat.lattice_dimension[2];
  z = mylat.lattice_dimension[3];
  N = mylat.sublattice;
  m = mylat.spin_orbital;

  if pos[1] > x || pos[2] > y || pos[3] > z || pos[4] > N || pos[5] > m
    error("position is out of range");
  elseif pos[1] < 1 || pos[2] < 1 || pos[3] < 1 || pos[4] < 1 || pos[5] < 1
    error("position is out of range");
  end

  idx = pos[5] + 
        m * (pos[4] - 1) + 
        m * N * (pos[3] - 1) + 
        m * N * z * (pos[2] - 1) + 
        m * N * z * y * (pos[1] - 1);

  return idx;
end


"""
  inverprocess of pos2index, this function returns the position of the lattice site.
  This function is mainly used for debugging, or to check the Hamiltonian matrix.
"""

function index2pos(mylat::LatticeSystem, idx::Int64)
  if length(idx) != 1
    error("idx must be a 1-element Int64");
  end

  x = mylat.lattice_dimension[1];
  y = mylat.lattice_dimension[2];
  z = mylat.lattice_dimension[3];
  N = mylat.sublattice;
  m = mylat.spin_orbital;

  pos = zeros(Int64, 5);

  pos[1] = div(idx - 1, m * N * z * y) + 1;
  pos[2] = div(mod1(idx, m * N * z * y) - 1, m * N * z) + 1;
  pos[3] = div(mod1(idx, m * N * z) - 1, m * N) + 1;
  pos[4] = div(mod1(idx, m * N) - 1, m) + 1;
  pos[5] = mod1(idx, m);

  return pos;
end

"""
  this function returns the position of the lattice site.
  The position is represented by a 4-element array. (2D lattice case is only implemented.)
  The first element is the x-coordinate, the second is the y-coordinate,
  Periodic boundary indexing is implemented.
"""

function neighboring(mylat::LatticeSystem, R₀::Array{Int64}, dR::Array{Int64})

  if length(R₀) != 3 || length(dR) != 3
    error("R₀ and dR must be a 3-element array");
  end

  x = mylat.lattice_dimension[1];
  y = mylat.lattice_dimension[2];
  z = mylat.lattice_dimension[3];

  if R₀[1] > x || R₀[2] > y || R₀[3] > z
    error("position is out of range");
  elseif R₀[1] < 1 || R₀[2] < 1 || R₀[3] < 1
    error("position is out of range");
  end

  R₁ = R₀ + dR;  R₁ = mod1.(R₁, (x,y,z));
  R₂ = R₀ - dR;  R₂ = mod1.(R₂, (x,y,z));

  return R₁, R₂;

end
