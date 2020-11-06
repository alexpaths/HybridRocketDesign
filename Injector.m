% Alexander Athougies
% CPSS
% 
% Calculates Nox mass flow through an injector
% 
% INPUTS:
% 1) Ninj = number of injectors
% 2) dinj = injector diameter (in)
% 3) Kinj = head loss coefficient
% 4) Psource = source pressure (psia)
% 5) Pdump = pressure of dump point (psia)
% 6) rho = liquid density (lbm/in^3)
% 
% OUT:
% 1) dm = mass flow (lbm/sec)

function dm = Injector(Ninj, dinj, Kinj, Psource, Pdump, rho)

rho = rho * (12^3); % lbm / ft^3

if Psource < Pdump
    dm = -(Ninj*pi/4*dinj^2)*(2.238*Kinj/rho/(Pdump - Psource))^(-1/2);
else
    dm = (Ninj*pi/4*dinj^2)*(2.238*Kinj/rho/(Psource - Pdump))^(-1/2);
end

end