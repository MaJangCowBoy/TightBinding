include("Define_Lattice.jl");
include("Define_Lattice_indexing.jl");
include("Define_Hamiltonian.jl");

xLen = 10; 
yLen = 10;
sublattice = 2;
intFreedom = 2;

mylat = 
  make_lattice(x = xLen, y = yLen, z = zLen, N = subLattice, m = intFreedom)

dR1 = [ 0, 0, 0];
dR2 = [+1, 0, 0];
dR3 = [ 0,-1, 0];

t = 1.0;

mylat = 
  make_hopping(mylat, t, dR1; N₀ = 1, m₀ = 1, N₁ = 2, m₁ = 1)

mylat = 
  make_hopping(mylat, t, dR2; N₀ = 1, m₀ = 1, N₁ = 2, m₁ = 1)

mylat = 
  make_hopping(mylat, t, dR3; N₀ = 1, m₀ = 1, N₁ = 2, m₁ = 1)


λ = 0.1;  nᵢⱼ = [0.0, 0.0, 1.0];

mylat =
  make_SOC_hopping(mylat, λ, nᵢⱼ, dR1; N₀ = 1, N₁ = 2)

mylat =
  make_SOC_hopping(mylat, λ, nᵢⱼ, dR2; N₀ = 1, N₁ = 2)
  
mylat =
  make_SOC_hopping(mylat, λ, nᵢⱼ, dR3; N₀ = 1, N₁ = 2)