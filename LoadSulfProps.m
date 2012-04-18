function LoadSulfProps(SulfConc)

global modelAtm

%Assign Sulfate Properties-----------------------------------------
modelAtm.Sulf.MW = 132.14;             %[=]g/mol
modelAtm.Sulf.rho = 1.77e12;  %[=]ug/m3, liquid
%modelAtm.Sulf.rho = 1e12;
%modelAtm.Sulf.MassConc = 8e-3;       %[=]ug/m3   calculated in Excel
%modelAtm.Sulf.MassConc = eps;
modelAtm.Sulf.MassConc = SulfConc;
modelAtm.Sulf.CollDiam = 4.20e-8;    %[=]cm, with air
modelAtm.Sulf.alpha = 1;             %sticking coeff
modelAtm.Sulf.NumConc = modelAtm.Sulf.MassConc*1e-12*1/(modelAtm.SOA.MW/6.023e23); %[=]#/cm3
z_Sulf=modelAtm.Sulf.MW/modelAtm.air.MW;      %ratio of molecular masses, S&P p401
%modelAtm.Sulf.lambda=(1/(pi*(1+z_Sulf)^(0.5)*modelAtm.air.NumConc*modelAtm.Sulf.CollDiam^2))*(1/100);  %[m]  S&P eq9.11 p401
modelAtm.Sulf.Diff = (1e-4)*3/(8*pi)*(pi*(1.38065e-16)^3*modelAtm.Temp^3*(1+z_Sulf)/...
    (2*(modelAtm.Sulf.MW/6.022e23)))^(1/2)/(modelAtm.Press*(modelAtm.Sulf.CollDiam)^2);
        %Diffusivity   [=]m2/s
modelAtm.Sulf.va = (8*1.38e-23*298/(pi*modelAtm.Sulf.MW*1e-3/6.022e23))^(1/2);

modelAtm.Sulf.lambda = 3*modelAtm.Sulf.Diff/modelAtm.Sulf.va;
%  Use Fuchs & Sutugin version