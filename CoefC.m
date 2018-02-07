function C = CoefC(S0,S)

% retroune directement le vecteur C avec tout les coef.
% > Article de C. Chen, P. Berini et al. "Efficient and accurate numreical
%   analysis of multilayer planar optical waveguides in lossy anisotropic
%   media", optics express, vol7(8), 2000.
% > Version dérecursivée (plus rapide)
% > utile pout FindZerosm.m

C = zeros(1, S0);
C(S0+1) = 1;

for n=1:S0
    kp = S0 - n +1;
    C(kp) = (1/(kp-1-S0)) *( S(1:n) * C( (kp+1):(S0+1) ).' ); 
end
