function Find_kwall(Dp)

global modelAtm

Cc = SlipCorr(Dp,modelAtm.Pop1.lambda);
Cc = SlipCorr(Dp,modelAtm.air.lambda);
Vs = modelAtm.Pop1.rho*Dp^2*9.8*Cc/(18*modelAtm.air.viscosity)*1e-9; %1e9 converts ug/m3 to kg/m3 
%settling velocity

x_Debye = pi*Vs/(2*sqrt(modelAtm.EddyDiffCoeff*modelAtm.Pop1.Diff));

D1 = Debye(x_Debye);

modelAtm.kwall = 6*sqrt(modelAtm.EddyDiffCoeff*modelAtm.Pop1.Diff)/...
    (pi*modelAtm.ChamberR)*D1 + Vs/(4*modelAtm.ChamberR/3);

