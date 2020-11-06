% Alexander Athougies
% CPSS
% 
% Calculates Nozzle Thrust and Isp given Dimensions and Mass Flow
% 
% INPUTS:
% 1) dm = mass flow expected(lbm/sec)
% 2) P0 = chamber pressure (psia)
% 3) T0 = chamber temperature (R)
% 4) At = throat area (in^2)
% 5) Ae = exit Area (in^2)
% 6) Aet = intake Area (in^2)
% 7) Patm = Atmospheric Pressure (psia)
% 8) Re = gas constant (in-lbf/lbm-R)
% 9) gamma = ratio of specific heats
% 
% OUT:
% 1) F = force (lbf)
% 2) Isp = specific impulse (sec)
% 3) dme = real mass flow (lbm/sec)

function [F, Isp,dme] = nozzlecalc(dm, P0, T0, At, Ae, Aet, Patm, Re, gamma)

Knoz = 80;
dme = 0;
a = 12*sqrt(gamma*Re*T0/12*32.2);
rho = sqrt(P0/Re/T0);
dt = 2*sqrt(At/pi());

% Find M1
M1 = dm/(a*rho*Aet);

% Find Theoretical Compression Ratio
E1 = 1 / M1*((2/(gamma+1))*(1+(gamma-1)/2*M1^2))^((gamma + 1) / 2 / (gamma - 1));  % 4.29 - Zucrow

% Find Real E
E = Aet / At;
E2r = Ae / At;

% check choking
if (E1*0.90) <= E
    dme = At*P0/sqrt(gamma*Re*T0)*(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)))); % 4.38 to 4.39 Zucrow
    Me = 1;
    E2 = 0;
    while E2 < E2r
       E2 = 1 / Me*((2/(gamma+1))*(1+(gamma-1)/2*Me^2))^((gamma + 1) / 2 / (gamma - 1));  % 4.29 - Zucrow
       Me = Me + .001; 
    end
    ve = Me*a; % in / sec
    Pe = P0*(1+(gamma - 1)/2*Me^2)^(-gamma / (gamma - 1)); % 4.26 Zucrow
%     while Pe < Patm
%         Pe = P0*(1+(gamma - 1)/2*Me^2)^(-gamma / (gamma - 1)); % 4.26 Zucrow
%         Me = Me - .001;
%     end
elseif E1 > E
    % no choked
    dme = Injector(1, dt, Knoz, P0, Patm, rho);
%     dme = At*P0/sqrt(gamma*Re*T0)*(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)))); % 4.38 to 4.39 Zucrow
    ve = dme/At/rho;
    Me = ve / a;
    Pe = P0*(1+(gamma - 1)/2*Me^2)^(-gamma / (gamma - 1));
end

% Find Thrust
if dm > 0 && dme > 0
F = dme*ve/32.2/12 + (Pe - Patm)*Ae; % lbf
Isp = F / dme; % sec
else
    F = 0;
    Isp = 0;
    dme = Injector(1, dt, 1, P0, Patm, rho);
end

end
