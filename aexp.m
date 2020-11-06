% Alexander Athougies
% Adiabatic Ideal Gas Expansion/Compression
%
% Assumes ideal gas of constant gamma

function P = aexp(V,Po,Vo,gamma)

C = Po*Vo^gamma;

P = C/(V^gamma);

end