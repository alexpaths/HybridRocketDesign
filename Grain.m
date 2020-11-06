% Alexander Athougies
% CPSS Hybrid Rocket Motor
% Initial Sizing
% 
% Vfuel = in^3
% L = in
% a,n = unitless
% dmox = lbm/sec
% OtoF = unitless
function [rinp, rin] = Grain(Vfuel, rout, L, a, n, dmox, OtoF)
    rhoHTPB = .0325; % lbm / in^3 - http://www.braeunig.us/space/propel.htm
    dmf = -dmox/OtoF;
    
    rinp = sqrt(rout^2 - (Vfuel/(pi*L)));
    
    rin = ((-1/2/L)*(1/a/rhoHTPB)*(dmf/(dmox^n)))^(1/(1-2*n));
end