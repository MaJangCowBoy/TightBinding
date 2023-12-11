"""
  This function returns the position of the lattice site.
  The position is represented by a 4-element array. (2D lattice case is only implemented.)
  The first element is the x-coordinate, the second is the y-coordinate,
  the third is the sublattice number, and the fourth is the sublattice degree of freedom.
"""

function pos2index(mylat::LatticeSystem, pos::Array{Int64})
  
  if length(pos) != 4
    error("pos must be a 4-element array");
  end

  x = mylat.lattice_dimension[1];   y = mylat.lattice_dimension[2];
  N = mylat.sublattice;  m = mylat.spin_orbital;

  if pos[1] > x || pos[2] > y || pos[3] > N || pos[4] > m
    error("position is out of range");
  elseif pos[1] < 1 || pos[2] < 1 || pos[3] < 1 || pos[4] < 1
    error("position is out of range");
  end

  idx = pos[4] + (pos[3] - 1) * m + (pos[2] - 1) * N * m + (pos[1] - 1) * y * N * m;

  return idx;
end

"""
  this function returns the position of the lattice site.
  The position is represented by a 4-element array. (2D lattice case is only implemented.)
  The first element is the x-coordinate, the second is the y-coordinate,
  Periodic boundary indexing is implemented.
"""

function neighboring(mylat::LatticeSystem, R₀::Array{Int64}, dR::Array{Int64})

  if length(R₀) != 2 || length(dR) != 2
    error("r must be a 2-element array");
  end

  x = mylat.lattice_dimension[1];   y = mylat.lattice_dimension[2];

  if R₀[1] > x || R₀[2] > y
    error("position is out of range");
  elseif R₀[1] < 1 || R₀[2] < 1
    error("position is out of range");
  end

  R₁ = R₀ + dR;  R₁ = mod1.(R₁, (x,y));
  R₂ = R₀ - dR;  R₂ = mod1.(R₂, (x,y));

  return R₁, R₂;

end
