% Alexander Athougies
% CPSS Hybrid
% Optimization
% 
% Monte Carlo Optimizer (for Rev 5 of Hybrid Sim)
% 
% OUTPUT = all species
% EOUTPUT = all species (optimized for Excel export)
% OPT = optimized cases
% EOPT = optimized cases (optimized for Excel export)
function [OUTPUT, EOUTPUT, OPT, EOPT] = MCHybridRev5()

%% Housekeeping
clc
close all
clear

%% OUTPUT Files
c = clock;
version = 11; % file version
filename = ['MCHybridRev5_Run',mat2str(version),'_RunLog_',mat2str(c(2)),mat2str(c(3)),mat2str(c(1)),'.log'];
of = fopen(filename,'w');
fprintf(of,['CPSS Hybrid Rocket Motor\nMonte Carlo Optimizer\nRun Log\nMCHybridRev4\nc. Alexander Athougies 2010\n\nLog Version: ',mat2str(c(2)),'.',mat2str(c(3)),'.',mat2str(c(1)),'  ',mat2str(c(4)),':',mat2str(c(5)),':',mat2str(c(6)),'\n\n']);
fprintf(['CPSS Hybrid Rocket Motor\nMonte Carlo Optimizer\nMCHybridRev4:\t',filename,'\nc. Alexander Athougies 2010\n\nRun: ',mat2str(c(2)),'.',mat2str(c(3)),'.',mat2str(c(1)),'  ',mat2str(c(4)),':',mat2str(c(5)),':',mat2str(c(6)),'\n\n']);


%% Simulation Variables
progressbar(0,'Hybrid Monte Carlo Optimizer')
nmax = 100; % number of cases to run
N = 5000; %points
perctol = .90; % Percentage Tolerance
maxplot = 10; % Maximum number of thrust curves to plot

%% Fixed
simtime = 10; % sec
rout = 3.66/2; % in
Pt = 750; % psia
Vt = 295; % in^3
Patm = 14.7; % psia
a = .2184; % Regression Variables - See Excel File
n = .33358;
NOSgamma = 1.2;
EXgamma = 1.3;
T = 2000; % deg R
Tatm = 530; % atmospheric Temp (R)

%% File Inputs
fprintf(of,'Grain Outer Radius: %3.3d in\nTank Pressure: %3.1d psia\nTank Volume: %3.1d in^3\nAtmospheric Pressure: %3.1d psia\nRegression Variables:\n\ta:\t%1.5d \n\tn:\t%1.5d\nRatios of Specific Heat:\n\tNOS:\t%1.1d\n\tExhaust:\t%1.1d\nChamber Temperature: %4.1d R\n\nAtmospheric Temperature:  %4.1d R\n',rout,Pt,Vt,Patm,a,n,NOSgamma,EXgamma,T,Tatm);
fprintf('Grain Outer Radius: %3.3d in\nTank Pressure: %3.1d psia\nTank Volume: %3.1d in^3\nAtmospheric Pressure: %3.1d psia\nRegression Variables:\n\ta:\t%1.5d \n\tn:\t%1.5d\nRatios of Specific Heat:\n\tNOS:\t%1.1d\n\tExhaust:\t%1.1d\nChamber Temperature: %4.1d R\nAtmospheric Temperature:  %4.1d R\n',rout,Pt,Vt,Patm,a,n,NOSgamma,EXgamma,T, Tatm);


%% Variables
mox = linspace(1,6,N); %lbm
OtoF = linspace(3,7,N);
L = linspace(8,15,N); % in
rt = linspace(.4,1,N); % in
% dinj are standard ANSI Drill Sizes (inches)
dinj = [1/16,1/8,3/32,1/2,1/4,.1285,.1360,.1405,9/64,.1440,.1470,.1495,.1520,.1540,.1562,.1570,.1590,.1610,.1660,.1695,.1719,.1730,.1770,.1800,.1820,.1850,.1875,.1890,.1910,.1935,.1960,.1990,.2010,.2031,.2040,.2055,.2090,.2130,.2187,.2210,.2280,.2420,.2460,.250,.2570,.2610,.2656,.2660,.2720,.2770,.2811,.2812,.29,.2950,.2968,.3020,.3125,.3160,.3230,.3261,.3320,.3390,.3437,.3480,.3580,.3594,.3680,.3750,.3770,.3860,.3906,.3970,.404,.4062,.4130,.4219,.4375,.4531,.4687,.4844,.5,.5156,.5312,.5469,.5625,.5781,.5937,.6094,.6250,.6406,.6562,.6719,.6875,.7031,.7187,.7344,.75,.7656,.7812,.7969,.8125,.8281,.8437,.8594,.8750,.8906,.9062,.9219,.9375,.9531,.9687,.9844,1]; % in
Ninj = linspace(1,20,N);
Kinj = linspace(.2,1,N);
E = linspace(2,7,N);

%% Preallocate
OUTPUT = cell(nmax,24);
EOUTPUT = cell(nmax,22);
progressbar(1/(nmax+1),'Hybrid Monte Carlo Optimizer')

%% Calculations
for index = 1:nmax
    %% Print Case to Screen
    fprintf(['\n\n\nCase #: ',mat2str(index),' of ',mat2str(nmax),'\n']);
    fprintf(of,['\n\n\nCase #: ',mat2str(index),' of ',mat2str(nmax),'\n']);
    
    %% Get Random
    moxc = mox(GetRand(N));
    OtoFc = OtoF(GetRand(N));
    Lc = L(GetRand(N));
    rtc = rt(GetRand(N));
    dinjc = dinj(GetRand(length(dinj)));
    Ninjc = round(Ninj(GetRand(N)));
    Kinjc = Kinj(GetRand(N));
    Ec = E(GetRand(N));
    
    %% Log Randoms
    fprintf(of,'\nOxidizer Mass:\t%2.3d lbm\nOxidizer to Fuel Ratio:\t%2.3d\nGrain Length:\t%2.3d in\nThroat Radius:\t%1.3d in\nInjector(s) Diameter:\t%1.3d in\nNumber of Injectors:\t%1.1d\nHead Loss Coefficient:\t%1.3d\nNozzle Expansion Ratio:\t%2.3d\n',moxc,OtoFc,Lc,rtc,dinjc,Ninjc,Kinjc,Ec);
    
    %% Calc
    [actburn, t, F, I, type, avgF, Isp] = HybridSimfRev5(simtime, moxc, OtoFc, Lc, rout, rtc, Pt, Vt, Patm, a, n, dinjc, Ninjc, Kinjc, NOSgamma, EXgamma, T, Ec, Tatm, of);
    
    %% Store Fixed
    fprintf(of,'\nActual Burn Time: %3.1d sec\n',actburn);
    
    OUTPUT{index,1} = rout;
    OUTPUT{index,2} = Pt;
    OUTPUT{index,3} = Vt;
    OUTPUT{index,4} = Patm;
    OUTPUT{index,5} = a;
    OUTPUT{index,6} = n;
    OUTPUT{index,7} = NOSgamma;
    OUTPUT{index,8} = EXgamma;
    OUTPUT{index,9} = T;
        
    EOUTPUT{index,1} = rout;
    EOUTPUT{index,2} = Pt;
    EOUTPUT{index,3} = Vt;
    EOUTPUT{index,4} = Patm;
    EOUTPUT{index,5} = a;
    EOUTPUT{index,6} = n;
    EOUTPUT{index,7} = NOSgamma;
    EOUTPUT{index,8} = EXgamma;
    EOUTPUT{index,9} = T;
    
    %% Store Moving
    OUTPUT{index,10} = moxc;
    OUTPUT{index,11} = OtoFc;
    OUTPUT{index,12} = Lc;
    OUTPUT{index,13} = rtc;
    OUTPUT{index,14} = dinjc;
    OUTPUT{index,15} = Ninjc;
    OUTPUT{index,16} = Kinjc;
    OUTPUT{index,17} = Ec;
    OUTPUT{index,18} = actburn;
    OUTPUT{index,19} = I;
    OUTPUT{index,20} = type;
    OUTPUT{index,21} = avgF;
    OUTPUT{index,22} = Isp;
    OUTPUT{index,23} = t;
    OUTPUT{index,24} = F;

    EOUTPUT{index,10} = moxc;
    EOUTPUT{index,11} = OtoFc;
    EOUTPUT{index,12} = Lc;
    EOUTPUT{index,13} = rtc;
    EOUTPUT{index,14} = dinjc;
    EOUTPUT{index,15} = Ninjc;
    EOUTPUT{index,16} = Kinjc;
    EOUTPUT{index,17} = Ec;
    EOUTPUT{index,18} = actburn;
    EOUTPUT{index,19} = I;
    EOUTPUT{index,20} = type;
    EOUTPUT{index,21} = avgF;
    EOUTPUT{index,22} = Isp;

    %% Progress
    progressbar((index+1)/(nmax+1),'Hybrid Monte Carlo Optimizer')
end

%% Find Best Impulse
dummy = zeros(1,nmax);
for i = 1:nmax
    dummy(i) = OUTPUT{i,19};
end
tolerance = (max(dummy)*perctol); % set tolerance here
index = find(dummy > tolerance);
if ~isempty(index)
OPT = cell(length(index),24);
EOPT = cell(length(index),22);
fprintf(of,'\nOptimal Case Numbers:\n');
for i = 1:length(index)
    for j = 1:24
        OPT{i, j} = OUTPUT{index(i), j};
    end
    for j = 1:22
        EOPT{i,j} = EOUTPUT{index(i),j};
    end
    fprintf(of,['\t',mat2str(index(i)),'\n']);
end


%% Plot thrust curves
if length(index) <= maxplot
figure('Name','Thrust Curves: Impulse')
hold on
M = cell(length(index),1);
for i = 1:length(index)
    c = GetColor();
    A = OPT{i,23};
    B = OPT{i,24};
    plot(A,B,'Color',c)
    M{i,1} = ['Case #',mat2str(index(i)),': ',OPT{i,20},mat2str(OPT{i,21}),' (Isp = ',mat2str(OPT{i,22}),')'];
end 
legend(M)
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Force (lbf)','FontSize',10,'FontWeight','bold')
c = clock;
title({'Thrust Curves of Best Impulse Species';['MCHybridRev5.Run: ',mat2str(version),' Log: ',mat2str(round(c(2))),'.',mat2str(round(c(3))),'.',mat2str(round(c(1)))]},'FontSize',12,'FontWeight','bold')
grid on
else
   fprintf(['Too Many Cases to Plot: Try increasing tolerance... # of Cases: ',mat2str(length(index)),'\n\n']); 
end

else
    OPT = 0;
    EOPT = 0;
end

%% Find Max Average Force
dummy = zeros(1,nmax);
for i = 1:nmax
    dummy(i) = OUTPUT{i,21};
end
tolerance = (max(dummy)*perctol); % set tolerance here
index = find(dummy > tolerance);
if ~isempty(index)
OPT2 = cell(length(index),24);
EOPT2 = cell(length(index),22);
for i = 1:length(index)
    for j = 1:24
        OPT2{i, j} = OUTPUT{index(i), j};
    end
    for j = 1:22
        EOPT2{i,j} = EOUTPUT{index(i),j};
    end
end


%% Plot thrust curves
if length(index) <= maxplot
figure('Name','Thrust Curves: Force')
hold on
M = cell(length(index),1);
for i = 1:size(index,2)
    c = GetColor();
    A = OPT2{i,23};
    B = OPT2{i,24};
    plot(A,B,'Color',c)
    M{i,1} = ['Case #',mat2str(index(i)),': ',OPT2{i,20},mat2str(OPT2{i,21}),' (Isp = ',mat2str(OPT2{i,22}),')'];
end
legend(M)
xlabel('Time (sec)','FontSize',10,'FontWeight','bold')
ylabel('Force (lbf)','FontSize',10,'FontWeight','bold')
c = clock;
title({'Thrust Curves of Best Average Force Species';['MCHybridRev5.Run: ',mat2str(version),' Log: ',mat2str(round(c(2))),'.',mat2str(round(c(3))),'.',mat2str(round(c(1)))]},'FontSize',12,'FontWeight','bold')
grid on
else
   fprintf(['Too Many Cases to Plot: Try increasing tolerance... # of Cases: ',mat2str(length(index)),'\n\n']); 
end

clear OPT2 EOPT2

end


%% Close Files
fclose('all');

%% Embedded Functions
    function OUT = GetRand(N)
        OUT = ceil(rand*N);
        if OUT > N
            OUT = N;
        elseif OUT < 1
            OUT = 1;
        end
    end

    function OUT = GetColor()
        OUT = [rand, rand, rand];
    end
end