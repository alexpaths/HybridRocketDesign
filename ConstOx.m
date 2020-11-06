% Alexander Athougies
% CPSS Hybrid Rocket Motor
% Initial Sizing
% 
% Assuming a constant mass oxidizer flow calc fuel mass flow over time
% 
function ConstOx()
close all
clc

%% Grain Length
rhoHTPB = .0325; % lbm / in^3 - http://www.braeunig.us/space/propel.htm
N = 500;
L = linspace(5,15,N); % in
rout = 3.5/2; % in
OtoF =6;
a = .2184; 
n = .33358;
Pt = 812; % psia
rt = .425; % in
At = pi*(rt^2); % in^2
rinp = zeros(N,1);
rin = zeros(N,1);

tburn = 7; % sec
minPc = 100; % psia
dPinj = 25; % psia
Vt = 295; % in^3
Vl = 3; % in^3
g = 1.2;


[Vox,Vf,mox,mf,dmox,dmf,T,I] = Hybrid(tburn,minPc,Pt,dPinj,Vt,Vl,OtoF,200,g);

V = linspace(Vt-Vox+Vl,Vt+Vl,N);
P = zeros(N,1);
for i = 1:N
    P(i) = aexp(V(i),Pt,Vt-Vox+Vl,g);
end

figure('Name',['Tank Pressure: Adiabatic Expansion... mox = ', mat2str(mox), ' lbm'])
hold on
plot(V,P-dPinj,'--')
plot(V,P)
legend('Chamber Pressure','Tank Pressure')
xlabel('Gas Volume in Tank (in^3)')
ylabel('Pressure (psia)')
grid on
hline(Pt,'r','Initial Tank Pressure')
hline(minPc,'r','Minimum Operating Pressure: Design')
vline(Vt+Vl,'r','Tank Volume')


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
hline(rout,'k','Maximum Radius')

index = find(rin > rout);
if ~isempty(index)
    vline(L(index(length(index))),'k','Ox Rich Burn')
    A = L(index(length(index)));
end
index = find(rin < rinp);
if ~isempty(index)
    vline(L(index(1)),'k','Fuel Rich Burn')
    B = L(index(1));
end

%% mass flow
if ~isempty(index)
 L = (A+B)/2;
[rinp, rin] = Grain(Vf, rout, L, a, n, dmox, OtoF);
Noztrends(Pt,dPinj,minPc,At,(dmf + dmox))

r = linspace(rinp,rout,N);
dmf = zeros(N,1);

for i = 1:N
    dmf(i) = rhoHTPB*2*pi*L*a*((dmox/pi)^n)*(r(i)^(1-2*n));
end

figure('Name',['Mass Flow through burn: L = ',mat2str(L),' in'])
plot(r,dmf)
xlabel('Port Radius (in)')
ylabel('Fuel Mass Flow (lbm/sec)')
grid on
vline(rin,'r','AVG Mass Flux')
end

% %% Actual Port Radius
% ODETIME = linspace(0,tburn,N);
% ODEINPUT = [rinp,14.7];
% R = 53.3533*12; % in-lbf / lbm-R
% T = 2000; % R
% g = 1.2;
% g0 = 32.2;
% Patm = 14.7;
% [t,y] = ode45(@odefcn,ODETIME,ODEINPUT);
% rin = y(:,1);
% Pch = y(:,2);
% dmf = zeros(size(t));
% ms = zeros(size(t));
% 
% for i = 1:length(t)
%     if rin(i) <= rout
%     dmf(i) = rhoHTPB*2*pi*L*a*((dmox/pi)^n)*(rin(i)^(1-2*n));
%     else
%         dmf(i) = 0;
%     end
%     ms(i) = Pch(i) * At / sqrt(g*R*T/g0) * sqrt((2*g^2)/(g-1)*((Patm/Pch(i))^(2/g))*(1-(Patm/Pch(i))^((g-1)/g)));
% end
% 
% figure('Name','Actual Fuel Mass Flow')
% plot(t,dmf)
% grid on
% xlabel('Time (sec)')
% ylabel('Fuel Mass Flow (lbm/sec)')
% figure('Name','Actual Total Mass Flow')
% hold on
% plot(t,dmf+dmox,'--')
% plot(t,ms)
% legend('Total Mass Flow','Nozzle Mass Flow')
% grid on
% xlabel('Time (sec)')
% ylabel('Mass Flow (lbm/sec)')
% figure('Name','Port Radius')
% plot(t,rin)
% grid on
% xlabel('Time (sec)')
% ylabel('Port Radius (in)')
% figure('Name','Chamber Pressure')
% plot(t,Pch)
% grid on
% xlabel('Time (sec)')
% ylabel('Chamber Pressure (psia)')
% 
%     function dy = odefcn(t,y)
%         dy = zeros(2,1);
%         A = pi * y(1)^2;
%         if y(1) >= rout
%             dy(1) = 0;
%             y(1) = rout;
%             dmf = 0;
%         else
%             dy(1) = a*(dmox/A)^n;
%             dmf = rhoHTPB*2*pi*L*a*((dmox/pi)^n)*(y(1)^(1-2*n));
%         end
%        Vt = pi*L*(rout^2) - (pi*L*(rout^2 - y(1)^2));
%        ms = y(2) * At / sqrt(g*R*T/g0) * sqrt((2*g^2)/(g-1)*((Patm/y(2))^(2/g))*(1-(Patm/y(2))^((g-1)/g)));
%        dy(2) = R*T*(dmox+dmf-ms)/Vt;
%     end
end
