
clear;
clc;

T = 298;                % K
P = 1.013e5;            % Pa
R = 8.314;              % J/mol-K or Pa/m3
NAv = 6.022e23;         % #/mol

M = P/(R*T)*1e-6*NAv;   % #/cm3

rho_SOA = 1e12;         % ug/m3

C_O3 = 50*1e-9*M;       % #/cm3
C_OH = 1e6;             % #/cm3
C_A = 0.2*1e-9*M;       % #/cm3
C_B = 0.2*1e-9*M;

k_aO3 = 1.01e-15*exp(-732/T);
k_aOH = 1.21e-11*exp(444/T);

k_bO3 = 1.74e-15*exp(-1297/T);
k_bOH = 2.38e-11*exp(357/T);

r_aO3 = k_aO3*C_O3*C_A;
r_aOH = k_aOH*C_OH*C_A;

Rate = r_aO3 + r_aOH;

MolRate = Rate/NAv;
MassRate = MolRate*152*1e6*1e6  % ug/m3-s

VolRate = MassRate/rho_SOA;     % m3/m3-s

DpRate = 2/(pi*(150*1e-9)^2)*VolRate*3600
