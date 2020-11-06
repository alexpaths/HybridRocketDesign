% Alexander Athougies
% CPSS Hybrid Rocket Motor Project
% Initial Sizing
% 
% INPUTS
%   tburn = sec
%   Pch = psia
%   Pt = psia
%   dPinj = psia
%   Vt = in^3
%   Vl = in^3
%   OtoF = NA
%   Isp = sec
%   g = NA
function [Vox,Vf,mox,mf,dmox,dmf,T,I] = Hybrid(tburn,Pch,Pt,dPinj,Vt,Vl,OtoF,Isp,g)

%% Given
% g0 = 32.2; % Gravity at Sea Level - ft/sec^2
rhoNOS = .025; % Density of Liquid Nitrous Oxide (lbm/in^3) - Matheson Data Sheet
rhoHTPB = .0325; % lbm / in^3 - http://www.braeunig.us/space/propel.htm
% Tt = 75; % Tank Temperature = STP (degF)

% % Assume Self Pressurized Tank
% Pt = 812; % psia - Vapor Pressure (NOX)

%% Calculations
% Find Tank maximum mass holding capability
mox = rhoNOS * (Vt + Vl) * (1 - ((Pch + dPinj)/Pt)^(1/g)); % lbm

% Find Fual Mass
mf = mox / OtoF; % lbm

% Find AVG Mass Flux
dmox = mox / tburn; % lbm/sec
dmf = mf / tburn;

% Find Volumes
Vox = mox / rhoNOS; % in^3
Vf = mf / rhoHTPB; 

% Thrust-to-Weight and Impulse
dmt = dmox + dmf;
Ft = Isp * dmt*32.2; %lbf
W = (mf + mox)*32.2; % lbf
T = Ft / W;  % Thrust to Weight - lbf/lbf;

I = Ft * tburn; % lbf-sec

end