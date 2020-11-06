% Alexander Athougies
% CPSS Hybrid Rocket
%
% Simulator (Rev. 8)...
%
% Assumptions: Isentropic, Adiabatic, Perfect Gas, Incompressible,
% Constant Burn Temps, Combustions stops with fuel/ox (aka. no problems)
%
% function [t,Pt,Pc,rin,mox,mf,F,I,type,avgF] = HybridSim()
function HybridSimRev8()

%% Housekeeping
clc
close all
clear all

%% Inputs
% % CPSS Hybrid
% Simulation
tburn = 10; % sec - simulation time

% Inital Conditions
mox = 6.8; % lbm - Initial Oxidizer Mass - In Tank
Vt = 295 + 24*pi/4*(.5^2); %in^3 - Tank Internal Volume + Line
x = .05; % Quality of Mixture ( % vapor )

% Grain Geometry
L = 13.465; % in - Grain Length
rout = 3.49/2; % in - Outer Radius of Grain
rinp = 1.8/2; % in - Inner Radius of Grain

% Nozzle Geometry
rt = .7/2; % in - Radius of Nozzle Throat
Vfree = 0; % in^3 - free space inside combustion chamber
E2 = 5; % Nozzle Area Expansion Ratio

% Atmospheric Conditions
Pc = 14.7; %psia - Initial Chamber Pressure
Patm = Pc; %psia - Atmospheric Pressure
Tatm = 505; % atmospheric Temp (R)
g0 = 32.2; % ft/sec^2 - Gravity

% Grain Burn Characteristics
a = .28; % Regression Variables - See Excel File
n = .3;

% Injector Geometry
dinj = 3/32; % in - diameter of injector(s)
N = 12; % number of injectors
Kinj = 44.25; % Injector Head Loss Coefficient

% Chemical Properties
Dome = load('vapordome_nox_v2.mat');
% % Loaded From File
% All Nitrous Vapor Dome Data under structure 'Dome'
g1 = Dome.g1; % Nitrous gamma
g2 = 1.2; % Exhaust gamma
Rnos = 378.5862; % in-lbf / lbm-R (.1889 J / g-K)
% Rex = function of O/F
rhoNOSl = table(Dome.Temperature_E, Dome.density_E_l, Tatm); % Density of Liquid Nitrous Oxide (lbm/in^3)
rhoNOSg = table(Dome.Temperature_E, Dome.density_E_g, Tatm); % Density of Gas Nitrous Oxide (lbm/in^3)
rhoHTPB = .0325; % lbm / in^3 - http://www.braeunig.us/space/propel.htm
rhoAIR = 4.34E-5; % lbm / in^3 - Density of Air at T = 20 degC
Pt = table(Dome.Temperature_E, Dome.Pressure_E, Tatm); % Tank Pressure (psia)

Dome.Af = 1.3;

%% Givens
% Nozzle Properties
Aet = pi*rout^2;
At = pi*rt^2;
Ae = At * E2;

%% Calculate Initial Conditions
Vox = (1 - x) * mox / rhoNOSl; % Oxidizer Volume (in^3)
Vf = (rout^2 - rinp^2)*pi*L; % Fuel Mass Volume (in^3)
mf = Vf * rhoHTPB; % lbm
Vu = Vt - Vox; % Initial Ullage Volume (in^3)
mair = (pi*L*rinp^2 + Vfree) * rhoAIR; % Mass of Air in Chamber (lbm)
moxg = x * mox; % lbm
moxl = (1-x) * mox; % lbm
ChamberTemp = Tatm; % R
Ptc = Pt; % Psia
dmoxg = 0;

%% Initial Internal Energy
hg = table(Dome.Temperature_E, Dome.enthalpy_E_g, Tatm); % BTU / lbm
hl = table(Dome.Temperature_E, Dome.enthalpy_E_l, Tatm); % BTU / lbm
ug = hg - Pt / rhoNOSg * 32.2 * 32.2 / 12 / 25037; % BTU / lbm
ul = hl - Pt / rhoNOSl * 32.2 * 32.2 / 12 / 25037; % BTU / lbm
Ut = moxg * ug + moxl * ul; % BTU

% %% In Progress Display
% figure('Name','Display','NumberTitle','off')
% hold on
% display = plot(0,Pt,'rx','MarkerSize',10);
% display2 = plot(0,Pc,'bx','MarkerSize',10);
% legend('Tank','Chamber')
% grid on
% axis([0 tburn 0 800])

%% Main Loop
% ODETIME = [0,tburn];
ODETIME = linspace(0,tburn,10000);
ODEINPUTS = [Ut, Pc, rinp, moxl, moxg, mf, mair];
ODEOPTIONS = odeset('RelTol',1E-3,'Events',@event);
[t,y] = ode45(@odefcn,ODETIME,ODEINPUTS,ODEOPTIONS);

%% Convert Outputs
Ut = real(y(:,1));
Pc = real(y(:,2));
rin = real(y(:,3));
moxl = real(y(:,4));
moxg = real(y(:,5));
mf = real(y(:,6));
mch = real(y(:,7));

%% Misc Calculations
F = zeros(size(t));
Isp = zeros(size(t));
dmox = zeros(size(t));
dmf = zeros(size(t));
T = zeros(size(t));
Ptc = zeros(size(t));

Ptc(1) = Pt;
rhoNOSl = table(Dome.Temperature_E, Dome.density_E_l, Tatm); % Density of Liquid Nitrous Oxide (lbm/in^3)
rhoNOSg = table(Dome.Temperature_E, Dome.density_E_g, Tatm); % Density of Gas Nitrous Oxide (lbm/in^3)

for i = 1:length(t)
    
    Voxc = moxl(i)/ rhoNOSl;
    
    % OxMass Flows
    Ttc = Ptc(i)*(Vt - Voxc)/(Pt*(Vt-Vox)/Tatm);
    moxc = moxl(i) + moxg(i);
    if moxl <= 0
        dmox(i) = Injector(N, dinj, Kinj, Ptc(i), Pc(i), rhoNOSg);
    else
        dmox(i) = Injector(N, dinj, Kinj, Ptc(i), Pc(i), rhoNOSl);
    end
        
    % Port Radius
    if rin(i) < rout
        Ap = pi * rin(i)^2;
        if dmox(i) > 0
            Gox = dmox(i) / 32.2 / Ap; % slugs/in^2-sec
            rp = a*Gox^n; % in/sec
        else
            rp = 0;
        end
    else
        rp = 0;
    end
    
    % Fuel Mass Flow
    dmf(i) = (2*pi*L * rin(i) * rp) * (rhoHTPB); % By rp
    
    OF = dmox(i) / dmf(i);
    
    if dmf(i) > 0
        T(i) = GetTemp(OF);
        ChamberTemp = T(i);
    elseif dmf(i) <= 0 && dmox(i) > 0
        T(i) = Ttc;
    else
        T(i) = ChamberTemp;
    end
    
    % Calculate Force / Isp Outside Loop
    [F(i), Isp(i), dmn] = nozzlecalc(dmf(i) + dmox(i), Pc(i), T(i), At, Ae, Aet, Patm, Rex, g2);
    
    % Calculate next Tank Pressure step
    if i < length(t)
        if moxl(i) > 0
            Vl = moxl(i) / rhoNOSl;
            Vu = Vt - Vl;
            rhoNOSg = moxg(i) / Vu; %lbm/in^3
            rhoNOSl = table(Dome.density_E_g, Dome.density_E_l, rhoNOSg); % lbm / in^3
            hg = table(Dome.density_E_g, Dome.enthalpy_E_g, rhoNOSg); % BTU/lbm
            hl = table(Dome.density_E_g, Dome.enthalpy_E_l, rhoNOSg); % BTU/lbm
            H = moxg(i)*hg + moxl(i)*hl;
            Ptc(i+1) = (H - Ut(i)) / Vt * (25037*12/32.2/32.2);
        else
            rhoNOSg = moxg(i) / Vt;
            hg = table(Dome.density_E_g, Dome.enthalpy_E_g, rhoNOSg); % BTU/lbm
            H = moxg(i) * hg;
            Ptc(i+1) = (H - Ut(i)) / Vt * (25037*12/32.2/32.2);
        end
    else
        Ptc(i) = Ptc(i-1);
    end

    progressbar(i / length(t),'Calculate Thrust')
end

%% PLOTS
figure('Name','Pressures')
hold on
plot(t, Ptc,'--')
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
title('Tank and Chamber Pressures','FontSize',12,'FontWeight','bold')
grid on
plot(t,Pc,'r')
ylabel('Pressure (psia)','FontSize',10,'FontWeight','bold')
legend('Tank','Chamber')

figure('Name','Grain Geometry')
plot(t,rin)
title('Port Radius','FontSize',12,'FontWeight','bold')
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Radius (in)','FontSize',10,'FontWeight','bold')
hline(rout,'r','Max Radius')
grid on

figure('Name','Masses')
hold on
grid on
plot(t, moxl)
plot(t, moxg,'r--')
plot(t, mf,'g-.')
legend('Ox Liquid','Ox Gas','Fuel')
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Mass (lbm)','FontSize',10,'FontWeight','bold')

figure('Name','Thrust')
subplot(2,1,1)
plot(t,F)
hold on
title('Thrust')
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Force (lbf)','FontSize',10,'FontWeight','bold')
grid on
plot(Dome.T1(:,1)-.235,Dome.T1(:,2),'r--')
plot(Dome.T2(:,1)-.235,Dome.T2(:,2),'g-.')
legend('Model','05/25/2010','06/01/2010')

subplot(2,1,2)
plot(t,Isp)
title('Isp')
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Isp (sec)','FontSize',10,'FontWeight','bold')
grid on

%% Outputs
save('output.mat')
clc
PeakThrust = real(max(F))
index = find(t<5);
AvgThrust = real(mean(F(index)))
Impulse = real(trapz(t(index),F(index)))
SpecImpulse = real(mean(Isp(index)))
OxMassFlow = real(mean(dmox(index)))
FuelMassFlow = real(mean(dmf(index)))

%% Embedded Fcns
    function dy = odefcn(t,y)
        %         global Ptc
        dy = zeros(7,1);
        
        %% Get Current #'s
        Utc = y(1);
        Pcc = y(2);
        rinc = y(3);
        moxl = y(4);
        moxg = y(5);
        mfc = y(6);
        mcc = y(7);
        Voxc = moxl / rhoNOSl;
        moxc = moxl + moxg;
        
        %% OxMass Flows
        Ttc = Ptc*(Vt - Voxc)/(Pt*(Vt-Vox)/Tatm);
        if moxl <= 0
            dmox = Injector(N, dinj, 8, Ptc, Pcc, rhoNOSg);
        else
            dmox = Injector(N, dinj, Kinj, Ptc, Pcc, rhoNOSl);
        end
        
        %% Port Radius
        Ap = pi * rinc^2;
        if rinc < rout
            if dmox > 0
                Gox = dmox / 32.2 / Ap; % slugs/in^2-sec
                rp = a*Gox^n; % in/sec
            else
                rp = 0;
            end
        else
            rp = 0;
        end
        dy(3) = rp;
        
        %% Fuel Mass Flow
        dmf = (2*pi*L * rinc * rp) * (rhoHTPB); % By rp
        dy(6) = -dmf;
        
        %% Get Combustion Temp
        OF = dmox / dmf;
        
        %         disp(mat2str(OF))
        if dmf > 0
            T = GetTemp(OF);
            ChamberTemp = T;
        elseif dmf <= 0 && dmox > 0
            T = Ttc;
        else
            T = ChamberTemp;
        end
        
        %% Nozzle Mass Flow
        if dmf == 0
            Rex = Rnos;
        else
            Rex = ExitR(OF); % in-lbf/lbm-R
        end
        [F, I, dmn] = nozzlecalc(dmf + dmox, Pcc, T, At, Ae, Aet, Patm, Rex, g2); % [lbf, sec, lbm/sec]
        
        %% Mass of Gas in Chamber
        if moxc <=0
            dmt = -dmn;
        elseif rinc >= rout && moxc > 0
            dmt = dmox - dmn;
        else
            dmt = dmox + dmf - dmn;
        end
        dy(7) = dmt;
        
        %% Chamber Pressure
        if moxc <= 0
            Vch = pi*L*rinc^2; % in^3
            Rair = 53.3533*12; % in-lbf / lbm-R
            dPcc = (Rair) * T / Vch * (-dmn); % psia/sec
        elseif rinc >= rout && moxc > 0
            Vch = pi*L*rinc^2;
            dPcc = (Rnos) * T / Vch * dmt;
        else
            Vc = pi*L*rinc^2 + Vfree;
            dVc = (2*pi*L * rinc * rp); 
            dPcc = Rex*T *(dmt*Vc - mcc*dVc) / (Vc^2);
%             dPcc = Rex*T*(dmt/Vc + mcc*dVc/(Vc^2));
        end
        dy(2) = dPcc;
              
        %% Tank Pressure
        if Ptc <= Patm
            Ptc = Patm;
            dUt = 0;
            dmoxg = 0;
            dmoxl = 0;
        elseif moxl > 0
            dmoxg = Dome.Af * dmox / rhoNOSl / (rhoNOSg/(rhoNOSg^2) - rhoNOSl/(rhoNOSl^2));
            dmoxl = -1*(dmoxg + dmox); % lbm/sec
            Vl = moxl / rhoNOSl;
            Vu = Vt - Vl;
            rhoNOSg = moxg / Vu; %lbm/in^3
            rhoNOSl = table(Dome.density_E_g, Dome.density_E_l, rhoNOSg); % lbm / in^3
            hg = table(Dome.density_E_g, Dome.enthalpy_E_g, rhoNOSg); % BTU/lbm
            hl = table(Dome.density_E_g, Dome.enthalpy_E_l, rhoNOSg); % BTU/lbm
            H = moxg*hg + moxl*hl;
            Ptc = (H - Utc) / Vt * (25037*12/32.2/32.2);
            dW = Ptc * dmox / rhoNOSl / 25037 * 32.2 / 12; % BTU / sec
%             dW = 0;
            dUt = -(dW + dmox*hl);  % BTU / sec
        else
            rhoNOSg = moxg / Vt;
            dmoxg = -dmox;
            dmoxl = 0;
            hg = table(Dome.density_E_g, Dome.enthalpy_E_g, rhoNOSg); % BTU/lbm
            H = moxg * hg;
            Ptc = (H - Utc) / Vt * (25037*12/32.2/32.2);
            dUt = -dmox*hg;
        end
        
        dy(1) = dUt;
        dy(4) = dmoxl;
        dy(5) = dmoxg; % lbm/sec
        
        %% Progress
        clc
        t
        dmn
        F
        I
        Pcc
        Ptc
        rp
        Rex
%         pause(1/32)
        progressbar(t/tburn,'Hybrid Simulator')
%         set(display,'XData',t,'YData',Ptc)
%         set(display2,'XData',t,'YData',Pcc)
    end

% Non-ODE Related Functions
    function T = GetTemp(OtoF)
        if OtoF <= .728
            T = 1768.8*OtoF + 510.33;
        elseif OtoF > 30
            T = ChamberTemp;
        else
            T = -7.3514*OtoF^2 - 311.15*OtoF + 5927.1;  % see excel: ISP.xls
        end
%         if OtoF <= 2
%             T = 643.86*OtoF+510.33;
%         elseif OtoF > 30
%             T = ChamberTemp;
%         else
%             T = -33.102*OtoF+1828.2;  % see excel: ISP.xls
%         end
    end

    function Rex = ExitR(OtoF)
       Rex = 8.314472 / (OtoF*44.013 + 54.09) * (OtoF + 1); % kJ / g-K
%         Rex = 1000*0.0945*OtoF^(-0.737); % kJ / g-K
        Rex = Rex * 0.5 / 4.448 * 39.37 * 453.59237; % in-lb / lbm-R
    end

    function [value,isterminal,direction] = event(t,r)
       value = Ptc - Patm;
       isterminal = 1;
       direction = 0;
    end
end