function LoadAtmos

global modelAtm

  


%Assign general atmos. and air properties
modelAtm.Temp = 300;
modelAtm.Press = 1.0135*10^5;  %Pa = N/m2 = kg/m-s2
modelAtm.air.MW = 28*0.79+32*0.21; % g/mol
modelAtm.air.NumConc = modelAtm.Press/(8.314*modelAtm.Temp)*6.022e23; %#/m3 
modelAtm.air.NumConc = modelAtm.Press/(8.314*modelAtm.Temp)*6.022e23; %#/m3 
modelAtm.air.viscosity = 1.8e-5*(modelAtm.Temp/298)^0.85;
modelAtm.air.lambda = 67e-9;

%modelAtm.CStarBasis = [1e-2 1e-1 1e0 1e1];
%modelAtm.CStarBasis = [1e-2 1e-1 1e0 1e1 1e2 1e3 1e4];
modelAtm.CStarBasis = [1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6];
%modelAtm.CStarBasis = [0.34]*10;
%modelAtm.CStarBasis = [7.85E-2];
%modelAtm.alphaProd = [1 0 0 0];

%%Chamber info
modelAtm.EddyDiffCoeff = 0.3; %1/s
modelAtm.StoV = 0.365; %1/m
modelAtm.Area = 29; %m2
modelAtm.ChamberR=1.4; %m FOR WALL LOSS CONST
%modelAtm.EddyDiffCoeff = 0.1; %1/s FOR WALL LOST CONST

modelAtm.kwall = 5e-5; %1/s  %For now, use a constant value for kwall
modelAtm.kwall = 1e-4; %Neil suggested this value, like Albert's experiments

%modelAtm.kwall = 0;
modelAtm.V_small = 0.1; % m3
modelAtm.V_big = 10; %m3

modelAtm.AgingResTime = 5;