"""
  function make_hamiltonian is constructed in a following way.
  Input: mylat::LatticeSystem, t::Float64, dr::Array{Int64}
  mylat is the lattice structure, t is the hopping parameter, 
  and dr is the relative position of the hopping site. 

  note that the hopping site is determined 
  by the relative position of the hopping site and the position of the lattice site.

  Output: mylat::LatticeSystem

  This function returns the lattice structure with the Hamiltonian matrix.
  The Hamiltonian matrix is stored in the field "Hamiltonian" of the lattice structure.
"""

function make_hamiltonian(mylat::LatticeSystem, t::Float64, dR::Array{Int64}; N₀ = 1, m₀ = 1, N₁ = 1, m₁ = 1)

  if length(dR) != 2
    error("dr must be a 2-element array");
  end

  x = mylat.lattice_dimension[1];   y = mylat.lattice_dimension[2];

  for i = 1:x, j = 1:y
    R₀ = [i,j,N₀,m₀];
    R₁, R₂ = neighboring(mylat, R₀[1:2], dR);
    append!(R₁, [N₁, m₁]);  append!(R₂, [N₁, m₁]);

    idxR₀ = pos2index(mylat, R₀);  idxR₁ = pos2index(mylat, R₁);  idxR₂ = pos2index(mylat, R₂);

    mylat.Hamiltonian[idxR₀,idxR₁] = -t;
    mylat.Hamiltonian[idxR₀,idxR₂] = -t;

  end

  return mylat;
end