function LoadSuspendedPopulation(Total, Tag,DiamSusp)
%Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
%Tag =1, specify total number concentration

global modelAtm

Dp0 = DiamSusp*1e-9;
modelAtm.Pop1.Dp0 = Dp0; %Mode of # distribution, (in m)

%Define population parameters

if modelAtm.AmmonSeedYN == 1
    FracOrg = 0; 
    rho = modelAtm.Sulf.rho;
else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
   % FracOrg = 0.5;
   % rho = modelAtm.SOA.rho*FracOrg + modelAtm.Sulf.rho*(1-FracOrg);
end

if Tag==1
    modelAtm.Pop1.NumConc0 = Total;
elseif Tag ==2
    modelAtm.Pop1.NumConc0 = Total/rho*6/pi*1/Dp0^3;
else 
    modelAtm.Pop1.NumConc0 = 5*1e9;
end

%modelAtm.Pop1.NumConc = NumConc0;
 
 PopString = int2str(1);

 NumConc = modelAtm.Pop1.NumConc0;
 Diam = modelAtm.Pop1.Dp0;

    [CondSinkBG BGCp MSulf] = LoadPopMono(1,NumConc,Diam,FracOrg,modelAtm.BGBin,1);

    modelAtm.Pop1.CondSink0=CondSinkBG;

    modelAtm.Pop1.Cp0 = BGCp;

    modelAtm.Pop1.MSulf0 = MSulf;
    modelAtm.Pop1.coverage = 0;

CpTot = zeros(1,modelAtm.NumBins);
for i = 1:modelAtm.NumBins 
    CpTot(i) = modelAtm.Pop1.Cp0(i);  
end