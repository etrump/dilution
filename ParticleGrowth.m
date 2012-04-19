function [dCv dCp] = ParticleGrowth(Cv, Cp, CondSink, j, KelvinYN,ROG,Dp)

global modelAtm

%j is the population number
n=modelAtm.NumBins;
KelvinYN = modelAtm.KelvinYN;
CStarBasis = modelAtm.CStarBasis;
alpha = modelAtm.SOA.alphaProd;

%CurPop = int2str(j);

%eval(['Dp = modelAtm.Pop' CurPop '.Dp;'])
UptakeFact = modelAtm.UptakeFact;
%UptakeFact = 1; % NO MASS XFER LIMITATIONS
    if KelvinYN==1
        Kelvin = KelvinTerm(j,Dp);
    else
        Kelvin = 1;     %do not consider Kelvin effect       
    end
    
for i=1:n
    if j==1 %suspended particles
        P(i) = alpha(i)*ROG;  %Production rate of species i (vapor), [=] ug/m3-s
    else
        P(i) = 0;              % Only produce vapors once
    end

    if sum(Cp)<=0 
        MolFrac(i) = 0;
    else
        MolFrac(i)=Cp(i)/(sum(Cp));  %Calculate mole fraction of organic i in particle phase
    end
    
    Phi(i)=UptakeFact*CondSink*(Cv(i)-Kelvin*MolFrac(i)*CStarBasis(i)); % calculte flux assuming no other limitations [=]ug/m3-s    
%    Phi(i) = -Kelvin*MolFrac(i)*CStarBasis(i)*UptakeFact*CondSink; % evap into vaccuum
%   Phi(i) = QPhi(i);
   if Dp*1e9 < 5
       Phi(i) = 0;
   end
   if Cp(i) < 1e-13
       Phi(i) = 0; 
   end
   
  CvAging(i) = 0; % no aging
  dCv(i)=P(i)-Phi(i)+CvAging(i);
  %dCp(i)=Phi(i)+CvAging(i);
  dCp(i) = Phi(i);    
      %Change in gas-phase concentration, or production minus flux to particles, [=] ug/m3-s
      %Change in particle-phase concentration, or flux to particle, [=] ug/m3-s
end
    



