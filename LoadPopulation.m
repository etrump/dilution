function [BGCp MSulf] = LoadPopulation(NumConc,DiamSusp, j)
%Diameter passed in in [m]
%# per m3


global modelAtm

%Define population parameters
%if modelAtm.AmmonSeedYN == 1%
%    FracOrg = 0; 
%    rho = modelAtm.Sulf.rho;
%else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
   % FracOrg = 0.5;
   % rho = modelAtm.SOA.rho*FracOrg + modelAtm.Sulf.rho*(1-FracOrg);
%end



%modelAtm.Pop1.NumConc = NumConc0;
 
 PopString = int2str(j);
 
 [CondSinkBG BGCp MSulf] = LoadPopBin(NumConc,DiamSusp,FracOrg);
 
 eval(['modelAtm.Pop' PopString '.CondSink0 = CondSinkBG;'])
 eval(['modelAtm.Pop' PopString '.Cp0 = BGCp;']) 
 eval(['modelAtm.Pop' PopString '.MSulf0 = MSulf;']) 
 eval(['modelAtm.Pop' PopString '.coverage = 0;']) 
 

%CpTot = zeros(1,modelAtm.NumBins);
%%for i = 1:modelAtm.NumBins 
%    CpTot(i) = modelAtm.Pop1.Cp0(i);  
%end