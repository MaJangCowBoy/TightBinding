using LinearAlgebra

"""
  function make_hopping is constructed in a following way.
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

"""
  function make_SOC_hopping is constructed in a following way.
  this makes an interaction between spin and orbital degree of freedom.
  spin-flip hopping is introduced.
  
  H(SOC hop from n to m)
    =  [ Cₘ⬆† ]ᵀ * [ H⬆⬆ H⬆⬇ ] * [cₙ⬆]
       [ Cₘ⬇† ]    [ H⬇⬆ H⬇⬇ ]   [Cₙ⬇] , where,

  (Hₙₘ)ᵢⱼ is defined as, (n,m = ↑ or ↓ ),

  (Hₙₘ)ᵢⱼ = iλ * (nᵢⱼ[1] * σ₁ + nᵢⱼ[2] * σ₂ + nᵢⱼ[3] * σ₃)

  where, σ₁, σ₂, σ₃ are Pauli matrices, and nᵢⱼ is a unit vector.
  detailed feature of nᵢⱼ depends on the system.

"""

function make_SOC_hopping(mylat::LatticeSystem, λ::Float64, nᵢⱼ::Vector{Float64}, dR::Array{Int64}; N₀ = 1, N₁ = 1)
  
  if length(dR) != 3
    error("dR must be a 3-element array");
  elseif length(λ) != 1
    error("λ must be a Float64 value");
  elseif length(nᵢⱼ) != 3
    error("nᵢⱼ must be a 3-element array");
  elseif N₀ == N₁ && dR == [0,0,0]
    error("R₀ and R₁ must be different");
  end

  x = mylat.lattice_dimension[1];
  y = mylat.lattice_dimension[2];
  z = mylat.lattice_dimension[3];

  σ₁ = [ 0   1; 
         1   0];
  σ₂ = [ 0 -im;
        im   0];
  σ₃ = [ 1   0;
         0  -1];

  nᵢⱼ = normalize(nᵢⱼ);

  for i = 1:x, j = 1:y, k = 1:z
    R₀⬆ = [i,j,k,N₀,1];
    R₀⬇ = [i,j,k,N₀,2];
    R₁  = neighboring(mylat, R₀⬆[1:end-2], dR);
    R₁⬆ = vcat(R₁, [N₁, 1]);
    R₁⬇ = vcat(R₁, [N₁, 2]);

    tSOC = im * λ * (nᵢⱼ[1] * σ₁ + nᵢⱼ[2] * σ₂ + nᵢⱼ[3] * σ₃);

    idxR₀⬆ = pos2index(mylat, R₀⬆);  idxR₀⬇ = pos2index(mylat, R₀⬇);
    idxR₁⬆ = pos2index(mylat, R₁⬆);  idxR₁⬇ = pos2index(mylat, R₁⬇);
    
    mylat.Hamiltonian[idxR₁⬆,idxR₀⬆] += tSOC[1,1];  mylat.Hamiltonian[idxR₁⬆,idxR₀⬇] += tSOC[1,2];
    mylat.Hamiltonian[idxR₁⬇,idxR₀⬆] += tSOC[2,1];  mylat.Hamiltonian[idxR₁⬇,idxR₀⬇] += tSOC[2,2];
    # remind the form of the Hamiltonian matrix,
    # H(SOC hop from n to m)
    #   = [ Cₘ⬆† ]ᵀ * im * λ [ H⬆⬆ H⬆⬇ ] * [cₙ⬆]
    #     [ Cₘ⬇† ]           [ H⬇⬆ H⬇⬇ ]   [Cₙ⬇]

    mylat.Hamiltonian[idxR₀⬆,idxR₁⬆] += conj(tSOC[1,1]);  mylat.Hamiltonian[idxR₀⬆,idxR₁⬇] += conj(tSOC[2,1]);
    mylat.Hamiltonian[idxR₀⬇,idxR₁⬆] += conj(tSOC[1,2]);  mylat.Hamiltonian[idxR₀⬇,idxR₁⬇] += conj(tSOC[2,2]);
    # since the Hamiltonian matrix is Hermitian,
    # H(SOC hop from n to m) = H(SOC hop from m to n)†

  end

  return mylat;
end