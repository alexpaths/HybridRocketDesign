% Alexander Athougies
% CPSS Hybrid
% 
% Structural Envelopes
% 
% INPUTS:
%   1) p = pressure (psi)
%   2) r = inner radius (in)
%   3) t = casing thickness (in)
function hybridstruct(p,r,t)

%% Given
% Aluminum
sigmamax = 35; % Max Tension (ksi) - Yield
sigmamin = 0; % Max Compression or 0 (ksi)
pr = .3; % Poission's Ratio

%% Burst Factors / MOS'
BF = 3;

%% Solve Casing
% Hoop
sigma1 = p*BF*r/t;
sigma2 = sigma1/2;

t_min = p*BF*r/sigmamax/1000;
tsigma = sigmamax/2;

% End Caps
Area = pi*r^2;
Force = p*BF * Area;
sigma3 = Force / (pi*((r+t)^2 - r^2));

% sigma2 = sigma2 + sigma3;

%% Plot Mohr's Circles
figure('Name','Structural Analysis')
% Aluminum
Mcirc(sigmamax, -sigmamin, 'r')
% Casing
Mcirc(sigma1/1000, sigma2/1000, 'g')
% Minimum thickness
Mcirc(sigmamax,tsigma,'b')

% Plot Properties
axis equal
legend('6061-T6 Aluminum',['Casing: BF = ',mat2str(BF)],['Minimum Thickness: ', mat2str(t_min)])
xlabel('Axial Stress (ksi)')
ylabel('Shear Stress (ksi)')
title({'Mohr''s Circles: Hybrid Rocket Casing',['Max Chamber Pressure: ',mat2str(round(p)),' psi']},'FontWeight','bold')

end