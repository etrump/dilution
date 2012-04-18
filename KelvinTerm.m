function [Ke] = KelvinTerm(j,Dp)

global modelAtm

CurPop = int2str(j);
%eval(['Dp = modelAtm.Pop' CurPop '.Dp;']);
%eval(['shapeFact = modelAtm.Pop' CurPop '.shapeFact;']);
shapeFact = 1;

Ke = exp(4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/...
    (8.314*modelAtm.Temp*shapeFact*modelAtm.SOA.rho*1e-6*Dp));
    

modelAtm.D_K = 4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/...
    (8.314*modelAtm.Temp*shapeFact*modelAtm.SOA.rho*1e-6);