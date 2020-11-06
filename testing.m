
close all

%% Inputs
OtoF = 6;
dm = 1.5;
At = pi/4*(.7^2);
Ae = pi/4*(1.65^2);
Aet = pi/4*(3.5^2);
Patm = 14.7;
gamma = 1.33;
P0 = linspace(0,1000,100);
T0 = linspace(400,6000,100);

%% Convert
Rex = 1000*0.0945*OtoF^(-0.737); % kJ / g-K
Rex2 = 8.314472 / (OtoF*44.013 + 54.09) * (OtoF + 1) * 1000; % kJ / g-K
Re = Rex * 0.5 / 4.448 * 39.37 * 453.59237; % in-lb / lbm-R
[P0, T0] = meshgrid(P0, T0);

%% Loop
F = zeros(size(P0));
Isp = zeros(size(P0));

for i = 1:100
    for j = 1:100
        [F(i,j), Isp(i,j),dme] = nozzlecalc(dm, P0(i,j), T0(i,j), At, Ae, Aet, Patm, Re, gamma);
        
    end
    progressbar(i/100,'Testing')
end

%% Plot
figure('Name','Force')
contour(P0, T0, F)
colorbar
xlabel('Chamber Pressure (psia)')
ylabel('Chamber Temperature (R)')

figure('Name','Isp')
contour(P0,T0,Isp)
colorbar
xlabel('Chamber Pressure (psia)')
ylabel('Chamber Temperature (R)')