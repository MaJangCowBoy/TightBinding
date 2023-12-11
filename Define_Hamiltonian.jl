"""
  function make_interaction is constructed in a following way.
  Input: mylat::LatticeSystem, t::Float64, dr::Array{Int64}
  mylat is the lattice structure, t is the hopping parameter, 
  and dr is the relative position of the hopping site. 

  note that the hopping site is determined 
  by the relative position of the hopping site and the position of the lattice site.

  Output: mylat::LatticeSystem

  This function returns the lattice structure with the Hamiltonian matrix.
  The Hamiltonian matrix is stored in the field "Hamiltonian" of the lattice structure.

  Note that, even you give t, as positive value,
  this leads to -t hopping term in the Hamiltonian matrix.
  Kinda confusing, but this is because of the definition of the Hamiltonian matrix.

"""

function make_hopping(mylat::LatticeSystem, t::Float64, dR::Array{Int64}; N₀ = 1, m₀ = 1, N₁ = 1, m₁ = 1)

  if length(dR) != 3
    error("dR must be a 3-element array");
  end

  x = mylat.lattice_dimension[1];
  y = mylat.lattice_dimension[2];
  z = mylat.lattice_dimension[3];

  for i = 1:x, j = 1:y, k = 1:z
    R₀ = [i,j,k,N₀,m₀];
    R₁ = neighboring(mylat, R₀[1:end-2], dR);
    append!(R₁, [N₁, m₁]);

    idxR₀ = pos2index(mylat, R₀);  idxR₁ = pos2index(mylat, R₁);

    if idxR₀ == idxR₁ && imag(t) != 0
      error("t must be real if the hopping site is the same as the lattice site");
    end

    mylat.Hamiltonian[idxR₀,idxR₁] = -t;
      # hopping from lattice site [i,j,k]   , sublattice N₀, spin/orbital m₀ 
      #           to lattice site [i,j,k]+dR, sublattice N₁, spin/orbital m₁ 

    if idxR₀ != idxR₁
      mylat.Hamiltonian[idxR₁,idxR₀] = -conj(t);
        # hopping from lattice site [i,j,k]   , sublattice N₀, spin/orbital m₀
        #           to lattice site [i,j,k]+dR, sublattice N₁, spin/orbital m₁
    end

  end

  return mylat;
end
