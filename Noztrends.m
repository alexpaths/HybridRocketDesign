% Alexander Athougies
% CPSS Hybrid 
% Initial Sizing
% 
% Nozzle Flow Trends

function Noztrends(Pt,dPinj,Pc,At,dm)
Pt = Pt - dPinj; %psia
% Pc = 300; % psia (minimum)

N = 100;
Patm = 14.7; %psia
g = 1.2;
T = 2000; % R
R = 53.3533; % ft lbf/ lbm-R
g0 = 32.2; % ft/sec^2

Rconv = R*12; %in lbf/lbm-R

Pch = linspace(Pt,Pc,N);
% At = 1; % in^2

% [Pch,At] = meshgrid(Pch,At);
ms = zeros(N,1);

for i = 1:N
%     for j = 1:N
        ms(i) = Pch(i) * At / sqrt(g*Rconv*T/g0) * sqrt((2*g^2)/(g-1)*((Patm/Pch(i))^(2/g))*(1-(Patm/Pch(i))^((g-1)/g)));
%     end
end

figure('Name','Nozzle Mass Flow Variation Due to Changing Chamber Pressure')
plot(Pch,ms)
legend('Actual MassFlow')
title(['Throat Area = ',mat2str(At),' in^2'])
grid on
xlabel('Chamber Pressure (psia)')
ylabel('Mass Flow (lbm/sec)')
hline(dm,'k','Designed Mass Flow')

end