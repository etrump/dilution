function [BGCv] = LoadEqmVapors2(CurrPop)

global modelAtm

PopString = int2str(CurrPop);
eval(['BGCp = modelAtm.Pop' PopString '.Cp0']);
eval(['BGDp = modelAtm.Pop' PopString '.Dp0']);
CStarBasis = modelAtm.CStarBasis;
n = modelAtm.NumBins;
modelAtm.Pop1.Kelvin0 = 1;


for i=1:n
   if modelAtm.KelvinYN==1
        Kelvin(i) = exp(4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/(8.314*298*modelAtm.SOA.rho*1e-6*BGDp));
        modelAtm.Pop1.Kelvin0 = Kelvin(i);
    else
        Kelvin(i) = 1;            
   end
    
   x(i) = BGCp(i)/sum(BGCp);
   BGCv(i) = x(i)*Kelvin(i)*CStarBasis(i);
   
  % if CurrPop==1 && modelAtm.AmmonSeedYN==1
  %     BGCv = eps*eps*ones(1,4);
   %    BGCv = zeros(1,4);
  % end
   
  % if CurrPop==2 && modelAtm.AmmonSeedYNPop2==1
  %     BGCv = eps*eps*ones(1,4);
  %     BGCv = zeros(1,4);
  % end
  % BGCv(i) = Kelvin(i)*CStarBasis(i);
  % BGCv(i) = 0;
end