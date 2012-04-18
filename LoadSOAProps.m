function LoadSOAProps()

global modelAtm

%Assign SOA Properties-----------------------------------------
modelAtm.SOA.MW = 152; %[=]g/mol
modelAtm.SOA.MW = 215;
modelAtm.SOA.MW = 200; % PAPER BASECASE D_K
%modelAtm.SOA.MW = 250; % PAPER MAX D_K
%modelAtm.SOA.MW = 100; % PAPER MIN D_K

modelAtm.SOA.rho = 1e12;              %[=]ug/m3, liquid
modelAtm.SOA.rho = 1.4e12; % PAPER BASECASE D_K
%modelAtm.SOA.rho = 1e12; %PAPER MAX D_K
%modelAtm.SOA.rho = 2e12; % PAPER MIN D_K

%%modelAtm.SOA.MassConc = 0.015;       %[=]ug/m3
modelAtm.SOA.CollDiam = 5.34e-8;     %[=]cm, with air
modelAtm.SOA.alpha = 1;              %sticking coeff
%%modelAtm.SOA.NumConc = MassConc_SOA*1e-12*1/(modelAtm.SOA.MW/6.023e23); %[=]#/cm3
z_SOA=modelAtm.SOA.MW/modelAtm.air.MW        %ratio of molecular masses, S&P p401
%modelAtm.SOA.lambda2=(1/(pi*(1+z_SOA)^(0.5)*modelAtm.air.NumConc*modelAtm.SOA.CollDiam^2))*(1/100);  %[=]m  S&P eq9.11 p401
modelAtm.SOA.Diff =(1e-4)*3/(8*pi)*(pi*(1.38065e-23)^3*modelAtm.Temp^3*(1+z_SOA)/...
    (2*(modelAtm.SOA.MW*1e-3/6.022e23)))^(1/2)/(modelAtm.Press*(modelAtm.SOA.CollDiam*1e-2)^2);
        %Diffusivity   [=]m2/s
modelAtm.SOA.Diff =3/(8*pi)*(pi*(1.38065e-23)^3*modelAtm.Temp^3*(1+z_SOA)/...
    (2*(modelAtm.SOA.MW*1e-3/6.022e23)))^(1/2)/(modelAtm.Press*(modelAtm.SOA.CollDiam*1e-2)^2);
modelAtm.SOA.va = (8*1.38e-23*298/(pi*modelAtm.SOA.MW*1e-3/6.022e23))^(1/2); %m/s
%modelAtm.SOA.va = 171

modelAtm.SOA.lambda = 3*modelAtm.SOA.Diff/modelAtm.SOA.va; %m
%modelAtm.SOA.lambda = 1e-7;
%  Use Fuchs & Sutugin version

modelAtm.SOA.sigma_guess = 0.05; % N/m &PAPERBASE CASE D_K
%%%%%modelAtm.SOA.sigma_guess = 0.03;
%modelAtm.SOA.sigma_guess = 0.1; % PAPER MAX D_K
%modelAtm.SOA.sigma_guess = 0.01; % PAPER MIN D_K



%modelAtm.SOA.alphaProd = [0 0 0 0];
%modelAtm.SOA.alphaProd(modelAtm.EmitBin) = 1;
%modelAtm.SOA.alphaProdAlphaP = [0.0284 0 0.3617 0.6099];


modelAtm.UptakeFact = 0.005;
modelAtm.UptakeFact = 0.1;
%modelAtm.UptakeFact = 1;

%modelAtm.SOA.InitialComp = [0.0090 0 0.1143 0.1919 0.2557 0.2695 0.1440]; % based on Ellis's data and partitioning calc
modelAtm.SOA.InitialComp = [0.0090 0 0.1143 0.1919 0.2557 0.2695 0.1440];
modelAtm.SOA.InitialComp = [0.0059 0 0.0753 0.1269 0.1756 0.2482 0.3126 0.0522];
modelAtm.SOA.InitialComp = [0.0753 0.1269 0.1756 0.2482 0.3126 0.0522]; %1 to 100,000
modelAtm.SOA.InitialComp = [1];
modelAtm.NumBins = length(modelAtm.CStarBasis)

%if modelAtm.UptakeFact >= 0.1
%    modelAtm.SOA.InitialComp(7) = 0;
%    modelAtm.SOA.InitialComp(6) = 0;
%    modelAtm.SOA.InitialComp(5) = 0; %C2
%end

%if modelAtm.UptakeFact >= 0.01
%    modelAtm.SOA.InitialComp(7) = 0;
%    modelAtm.SOA.InitialComp(6) = 0; %C2
%end

%if modelAtm.UptakeFact >= 0.001
%    modelAtm.SOA.InitialComp(7) = 0; %C2
%end

%if modelAtm.AlphaPYNVap == 1
    %modelAtm.SOA.alphaProd = [0.0284 0 0.3617 0.6099]; %This is the alpha-pinene distribution (considering only first 4 bins)
%    modelAtm.SOA.alphaProdAlphaP = [0.0627 0 0.5873 0.2899 0.0501 0.0078]; %based on Ellis's data and partitioning calc
    
    
    %modelAtm.SOA.alphaProd = [0 0 0.8 0.1]; %For testing!
%end

modelAtm.SOA.alphaProd = [0.004 0.000 0.051 0.086 0.120 0.183 0.400 0.350 0.200];
%modelAtm.SOA.alphaProd = [1];

end

