% Alexander Athougies
% CPSS Hybrid Rocket Motor Project
% Initial Sizing
% 
% Design Trends
%
%
function HybridTrends()

clc
close all

%% Set Inputs
minTburn = 7; %sec
maxTburn = 10;
Pt = 600; %psia
minPc = 300;
maxPc = 500;
dPinj = 25;
Isp = 200;
OtoF = 5;

%% Givens/Fixed
Vt = 295; % in^3
Vl = 0;
g = 1.3; % gamma
N = 100; %# of Points

%% Calculate
A = linspace(minTburn,maxTburn,N);
B = linspace(minPc,maxPc,N);
[tburn,Pch] = meshgrid(A,B);
mox = zeros(N,N);
mf = zeros(N,N);
T = zeros(N,N);
I = zeros(N,N);
dmox = zeros(N,N);
dmf = zeros(N,N);

for i = 1:N % varying tburn - x
    for j = 1:N % varying Pch - y
        [Vox,Vf,mox(i,j),mf(i,j),dmox(i,j),dmf(i,j),T(i,j),I(i,j)] = Hybrid(tburn(i,j),Pch(i,j),Pt,dPinj,Vt,Vl,OtoF,Isp,g);
        progressbar((i + j/N)/(N+1),'Hybrid Trends')
    end
end

%% PLOT
figure('Name','Impulse')
contourf(tburn,Pch,I)
colorbar
xlabel('Burn Time (sec)')
ylabel('Chamber Pressure (psia)')
title('Total Impulse (Est.) (lb-sec)')
figure('Name','Thrust to Weight')
contourf(tburn,Pch,T)
colorbar
xlabel('Burn Time (sec)')
ylabel('Chamber Pressure (psia)')
title('Thrust to Weight')
figure('Name','Oxidizer Mass Available')
contourf(tburn,Pch,mox)
colorbar
xlabel('Burn Time (sec)')
ylabel('Chamber Pressure (psia)')
title('Oxidizer Mass (lbm)')


% figure('Name','Test')
% contourf(tburn,Pch,(mox/max(max(mox))+I/max(max(I))+T/max(max(T)) + dmox/max(max(dmox))))
% colorbar
% xlabel('Burn Time (sec)')
% ylabel('Chamber Pressure (psia)')

figure('Name','Ox Mass Flux')
contourf(tburn,Pch,dmox)
colorbar
xlabel('Burn Time (sec)')
ylabel('Chamber Pressure (psia)')
title('AVG Oxidizer Mass Flux Available (lbm/sec)')
figure('Name','Fuel Mass Flux')
contourf(tburn,Pch,dmf)
colorbar
xlabel('Burn Time (sec)')
ylabel('Chamber Pressure (psia)')
title('AVG Fuel Mass Flux (lbm/sec)')

% figure
% plot(Pch(:,50),I(:,50))

%% Grain Length
N = 500;
L = linspace(5,15,N); % in
rout = 3.5/2; % in
OtoF = 5;
a = .195;
n = .325;
rinp = zeros(N,1);
rin = zeros(N,1);

[Vox,Vf,mox,mf,dmox,dmf,T,I] = Hybrid(8,300,812,25,295,0,OtoF,200,1.3);

for i = 1:N
 [rinp(i), rin(i)] = Grain(Vf, rout, L(i), a, n, dmox, OtoF);
end

figure('Name','Inner Grain port Radius')
hold on
plot(L,rin,'r')
plot(L,rinp,'b')
legend('Mass Flow','Physical')
xlabel('Length (in)')
ylabel('Port Radius (in)')
grid on


end