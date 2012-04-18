function [TotalCS BGCp TotalSulfMass] = LoadPopBin(TotalN, ModeM, FracOrg)

%Mode = ModeM*1e9; %convert mode diameter from meters to nanometers
Dp = ModeM;
Dp_nm = ModeM*1e9;

%Dp = ModeM;
global modelAtm

rho_SOA = modelAtm.SOA.rho;
rho_Sulf = modelAtm.Sulf.rho;
rho = rho_SOA*FracOrg + rho_Sulf*(1-FracOrg);

lambda = modelAtm.SOA.lambda; %m
Diff = modelAtm.SOA.Diff; %m2/s

Kn = 2*lambda*1e9/Dp;
Beta = FuchsC(Kn);
   
TotalV = TotalN*pi/6*Dp^3; % m3/m3
TotalCS = TotalN*2*pi*Dp*Diff*Beta; %[=] 1/s

TotalOrgMass = FracOrg*TotalV*rho_SOA;  % [=] ug/m3

TotalSulfMass = (1-FracOrg)*TotalV*rho_Sulf; % [=] ug/m3
TotalSulfMass = (1-FracOrg)*TotalV*rho;

CStarBasis = modelAtm.CStarBasis;

BGCp = zeros(1,modelAtm.NumBins);   %BGCp [=] ug/m3
BGCp = TotalOrgMass*modelAtm.SOA.InitialComp;





