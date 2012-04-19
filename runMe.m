clear
clc
close all
global modelAtm trackCov MassforPlot MassforPlotAlpha Struct trackMode
global Tagg CpT DpT  CvT
global trackCondSink

trackCondSink = zeros(1,11);
trackMode = zeros(1,4);

%modelAtm.NumBins=4;
modelAtm.NumBins = 9;
modelAtm.TotalPop = 10;
n = modelAtm.NumBins;
TotalPop = modelAtm.TotalPop;

%Set simulation options: 1=yes, 0 = no ----------------------------
modelAtm.AgingYN = 0; 

modelAtm.AmmonSeedYN = 0;
modelAtm.AmmonSeedYNPop2 = 1;
modelAtm.PulseYN = 1;
modelAtm.KelvinYN = 0; %!!!!!!! turn off to increase speed
modelAtm.KelvinYN = 1; %
modelAtm.AlphaPYN = 1; %Alpha pinene organic = 1; single-bin organic = 0;
modelAtm.AlphaPYNPop2 = 0;
modelAtm.AlphaPYNVap = 0;
modelAtm.CoagYN = 0;

%Set simulation parameters ----------------------------
TimeVector = [0 7]*3600;

EmitTime = 3*3600;
modelAtm.EmitTime = EmitTime;
modelAtm.DecayConstant = 3;     %hr-1
modelAtm.AgingResTime = 5;      %hr  

modelAtm.BGBin = 3; %Specify C* bin of initial particles; 1 is C*=1e-2, 2 is C* =1e-1, 3 is C*=1e0, 4 is C*=1\  /X\( .. )/X\ 
modelAtm.EmitBin = 1;

SulfMassConc = 0; %Sulf concentration in vapor (held constant for now)
modelAtm.ROG = 1e-55;% reacting organic rate [=] ug/m3-s
%LoadSOAProps(ROG, EmitTime, EmitBin);
LoadAtmos;
LoadSOAProps;
%modelAtm.Injection = 1287; %120307.... 120202
modelAtm.Injection = 429; %120326... 120409
modelAtm.V_small = 2;
[Caer_part, Cvap_part] = Partition(modelAtm.Injection/modelAtm.V_small)

for i = 1:length(Caer_part)
    Caer_frac(i) = Caer_part(i)/sum(Caer_part);
end

modelAtm.Caer_part = Caer_part;
modelAtm.Caer_frac = Caer_frac;
modelAtm.Cvap_part = Cvap_part;
modelAtm.Caer_tot_initial = sum(Caer_part);
modelAtm.SOA.InitialComp = modelAtm.Caer_frac;
modelAtm.NumBins = length(modelAtm.SOA.InitialComp);
n = modelAtm.NumBins;
modelAtm.SOA.lambda = 100e-9;
LoadSulfProps(SulfMassConc);

%TotalPop = 20;

DF = 150;
DF = 10; %120409

modelAtm.DF = DF;

%[Dp_bin_vect, N_tot_bin] = LoadData120202(TotalPop);
%[Dp_bin_vect, N_tot_bin] = LoadData120307(TotalPop);
%[Dp_bin_vect, N_tot_bin] = LoadData120326(TotalPop);
[Dp_bin_vect, N_tot_bin] = LoadData120409(TotalPop);
for i = 1:length(Dp_bin_vect)
Mass_bin_vect(i) = N_tot_bin(i)*pi/6*modelAtm.SOA.rho*(Dp_bin_vect(i)*1e-9)^3*1e6*1/DF;
end
modelAtm.sumMass_bin = sum(Mass_bin_vect);

DiamSusp = Dp_bin_vect;
TotalSusp = N_tot_bin*1/DF*1e6;

NumConc0 = TotalSusp;
DiamSusp_m = DiamSusp*1e-9;
modelAtm.Pop = TotalPop;

for j = 1:TotalPop
    [BGCp_1bin MSulf_1bin] = LoadPopulation(TotalSusp(j),DiamSusp_m(j),j);
     
    for i = 1:modelAtm.NumBins
     Cp0((j-1)*modelAtm.NumBins+i) = BGCp_1bin(i);  % only for one-bin organics!
     Cp0_stack(j,i) = BGCp_1bin(i);
    end
     MSulf0(j) = MSulf_1bin;
       
    PopString = int2str(j);
    eval(['modelAtm.Pop' PopString '.Dp0 = DiamSusp_m(j);']) 
    eval(['modelAtm.Pop' PopString '.NumConc0 = TotalSusp(j);'])
    
    %Cv0_contrib(j,:) = LoadEqmVapors2(j);
    LoadPopulationProps(j);
    
end

modelAtm.sumCp0_stack = sum(Cp0_stack);

%DF = 1; % for checking

%modelAtm.Cv0 = Cv0_contrib(14,:)*1/DF; % just choose size bin with most particles for now
%modelAtm.Cv0 = Cv0_contrib(ceil(modelAtm.Pop/2+1),:)*1/DF; % just choose size bin with most particles for now
modelAtm.Cv0 = Cvap_part*1/DF;
Cv0 = modelAtm.Cv0;
Dp0 = DiamSusp;

% if (TotalPop == 20)
% Dp0=[modelAtm.Pop1.Dp0 modelAtm.Pop2.Dp0 modelAtm.Pop3.Dp0 modelAtm.Pop4.Dp0 modelAtm.Pop5.Dp0 modelAtm.Pop6.Dp0 modelAtm.Pop7.Dp0 modelAtm.Pop8.Dp0 ...
%     modelAtm.Pop9.Dp0 modelAtm.Pop10.Dp0 modelAtm.Pop11.Dp0 modelAtm.Pop12.Dp0 modelAtm.Pop13.Dp0 modelAtm.Pop14.Dp0 modelAtm.Pop15.Dp0 modelAtm.Pop16.Dp0 ...
%     modelAtm.Pop17.Dp0 modelAtm.Pop18.Dp0 modelAtm.Pop19.Dp0 modelAtm.Pop20.Dp0]*1e9; %meters to nanometers%
% MSulf0=[modelAtm.Pop1.MSulf0 modelAtm.Pop2.MSulf0 modelAtm.Pop3.MSulf0 modelAtm.Pop4.MSulf0 modelAtm.Pop5.MSulf0 modelAtm.Pop6.MSulf0 modelAtm.Pop7.MSulf0 modelAtm.Pop8.MSulf0 ...
%     modelAtm.Pop9.MSulf0 modelAtm.Pop10.MSulf0 modelAtm.Pop11.MSulf0 modelAtm.Pop12.MSulf0 modelAtm.Pop13.MSulf0 modelAtm.Pop14.MSulf0 modelAtm.Pop15.MSulf0 modelAtm.Pop16.MSulf0 ...
%     modelAtm.Pop17.MSulf0 modelAtm.Pop18.MSulf0 modelAtm.Pop19.MSulf0 modelAtm.Pop20.MSulf0];
% Cp0=[modelAtm.Pop1.Cp0 modelAtm.Pop2.Cp0 modelAtm.Pop3.Cp0 modelAtm.Pop4.Cp0 modelAtm.Pop5.Cp0 modelAtm.Pop6.Cp0 modelAtm.Pop7.Cp0 modelAtm.Pop8.Cp0 ... 
%     modelAtm.Pop9.Cp0 modelAtm.Pop10.Cp0 modelAtm.Pop11.Cp0 modelAtm.Pop12.Cp0 modelAtm.Pop13.Cp0 modelAtm.Pop14.Cp0 modelAtm.Pop15.Cp0 modelAtm.Pop16.Cp0 ...
%     modelAtm.Pop17.Cp0 modelAtm.Pop18.Cp0 modelAtm.Pop19.Cp0 modelAtm.Pop20.Cp0];
% elseif (TotalPop==10)
%     Dp0=[modelAtm.Pop1.Dp0 modelAtm.Pop2.Dp0 modelAtm.Pop3.Dp0 modelAtm.Pop4.Dp0 modelAtm.Pop5.Dp0 modelAtm.Pop6.Dp0 modelAtm.Pop7.Dp0 modelAtm.Pop8.Dp0 ...
%     modelAtm.Pop9.Dp0 modelAtm.Pop10.Dp0]*1e9;
% MSulf0=[modelAtm.Pop1.MSulf0 modelAtm.Pop2.MSulf0 modelAtm.Pop3.MSulf0 modelAtm.Pop4.MSulf0 modelAtm.Pop5.MSulf0 modelAtm.Pop6.MSulf0 modelAtm.Pop7.MSulf0 modelAtm.Pop8.MSulf0 ...
%     modelAtm.Pop9.MSulf0 modelAtm.Pop10.MSulf0];
% Cp0=[modelAtm.Pop1.Cp0 modelAtm.Pop2.Cp0 modelAtm.Pop3.Cp0 modelAtm.Pop4.Cp0 modelAtm.Pop5.Cp0 modelAtm.Pop6.Cp0 modelAtm.Pop7.Cp0 modelAtm.Pop8.Cp0 ... 
%     modelAtm.Pop9.Cp0 modelAtm.Pop10.Cp0];
% end

%Initialize -----------------------------------

tspan = TimeVector;
y0=[Cv0 Cp0 Dp0 MSulf0 NumConc0];
[T,Y] = Main3(tspan,y0);


Yagg = Y;
Tagg = T;

Struct.Yagg = Yagg;
Struct.Tagg = Tagg;
%Struct.EndIndex = EndIndex;

ProcessData;
%CreateFigS;
save('firstfile.mat')

