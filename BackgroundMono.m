function [TotalCS BGCv BGCp TotalSulfMass Mode] = BackgroundCS(Tag, Total, Mode)

%Tag = 1; ;
%Total = 10000;
%Mode = 300;
%logSigma = 0.3;

global modelAtm


rho_SOA = modelAtm.SOA.rho;
rho_Sulf = modelAtm.Sulf.rho;
lambda = modelAtm.Sulf.lambda;
Diff = modelAtm.Sulf.Diff;

%rho = rho_SOA*1e-18; % [=]ug/um3

 % Mode given in nm
Diff = 7e-2; %cm2/s
lambda = 1.8e-8; % [=]m
%Integrate over entire particle distribution
Dp_0 = 10; %[=]nm, start at 0.010 um
Dp_f = 1e5;   %end at 100 um

% Dp = Dp_0:10:Dp_f;
%  for i = 1:length(Dp)
%      logDp(i) = log(Dp(i));
%  end


logSigma = 0.05;  % choose a standard deviation for the distribution (see S&P fig 8.3 pg371 for typical values)
Sigma = exp(logSigma);

% Load particle size range 
logDp_0 = log(Dp_0);
logDp_f = log(Dp_f);
logDp = logDp_0:0.01:logDp_f;
 for i = 1:length(logDp)
     Dp(i) = exp(logDp(i));
 end

 Kn = zeros(length(Dp),1);
 logMode = log(Mode);
 
%  if Tag==1
%      Dpg = Mode
%      DpgS = exp(logMode+2*logSigma^2)
%      DpgV = exp(logMode+3*logSigma^2)
%  elseif Tag==2
%      Dpg = exp(logMode-2*logSigma^2)
%      DpgS = Mode;
%      DpgV = exp(log(Dpg)+3*logSigma^2)
%  elseif (Tag==3)||(Tag==4)
%      Dpg = exp(logMode-3*logSigma^2)
%      DpgS = exp(log(Dpg)+2*logSigma^2)
%      DpgV = Mode
%  end
 
 
  for i = 1:length(Dp)-1
   % if Dp(i)==Mode
    d(i) = Total/((2*pi)^(1/2)*logSigma)*exp(-(logDp(i)-logMode)^2/(2*(logSigma)^2)); 
    %d(i) = Total; %S&P pg 370, eq 8.54
    Kn(i) = 2*lambda*1e9/Dp(i);
    Beta(i) = FuchsC(Kn(i));
    %Beta(i) = 1;
    if Tag==1  %Number distribution given    
        dCS(i) = d(i)*2*pi*Dp(i)*Diff*Beta(i)*1e-7; %[=] 1/s
        n_N(i) = d(i);  %#/cm3 ?
        n_S(i) = d(i)*pi*Dp(i)^2*1e-6; %pi*nm^2*(um/1000nm)^2
        n_V(i) = d(i)*pi/6*Dp(i)^3*1e-9;
    elseif Tag==2    
        dCS(i) = d(i)*2/Dp(i)*Diff*Beta(i)*1e-1; %[=] 1/s
        n_N(i) = d(i)/(pi*Dp(i)^2)*1e6; %#/cm3
        n_S(i) = d(i);  %um2/cm3
        n_V(i) = d(i)*Dp(i)/6*1e-3;
    elseif Tag==3
        dCS(i) = d(i)*12/Dp(i)^2*Diff*Beta(i)*1e2;
        n_N(i) = d(i)*6/pi*1/Dp(i)^3*1e9;
        n_S(i) = d(i)*6/Dp(i)*1e3;
        n_V(i) = d(i);
    elseif Tag==4
        Dpg = 10^(log10(Mode)-3*logSigma^2)
        dCS(i) = d(i)*1e-6/rho*12/Dp(i)^2*Diff*Beta(i)*1e2;
        n_N(i) = d(i)*1e-6/rho*6/pi*1/Dp(i)^3*1e9;
        n_S(i) = d(i)*1e-6/rho*6/Dp(i)*1e3;
        n_V(i) = d(i)*1e-6/rho;
    end
        
    CS(i) = dCS(i)*(logDp(i+1)-logDp(i));
    D(i) = d(i)*(logDp(i+1)-logDp(i));
    V(i) = n_V(i)*(logDp(i+1)-logDp(i));
    N(i) = n_N(i)*(logDp(i+1)-logDp(i));
    S(i) = n_S(i)*(logDp(i+1)-logDp(i));
%     else
%         dCS(i) = 0;
%         n_N(i) = 0;
%         n_S(i) = 0;
%         n_V(i) = 0;
%         CS(i) = 0;
%         D(i) = 0;
%         V(i) = 0;
%         N(i) = 0;
%         S(i) = 0;
%     end
  end  

  dCS(length(Dp)) = 0;
  d(length(Dp)) = 0;
  dN(length(Dp)) = 0;
  dV(length(Dp)) = 0;
  
  %CS(length(Dp))=CS(length(Dp)-1);
  %D(length(Dp))=D(length(Dp)-1);
  %V(length(Dp))=V(length(Dp)-1);
  %S(length(Dp))=S(length(Dp)-1);
  
  n_N(length(Dp))=n_N(length(Dp)-1);
  n_S(length(Dp))=n_S(length(Dp)-1);
  n_V(length(Dp))=n_V(length(Dp)-1);
  
TotalCS = sum(CS)
TotalD = sum(D)  % check if program is working correctly
TotalV = sum(V)
TotalS = sum(S)
TotalN = sum(N)


FracOrg = 0.3;
rho = rho_SOA*1e-18; %ug/um3
TotalOrgMass = FracOrg*TotalV*rho*1e6  % [=] ug/m3
rho = rho_Sulf*1e-18; %ug/um3
TotalSulfMass = (1-FracOrg)*TotalV*rho*1e6; % [=] ug/m3

CStarBasis = modelAtm.CStarBasis;

%BGCp(1) = TotalOrgMass  %BGCp, BGCv [=] ug/m3
BGCp(1) = 0;
BGCp(2) = 0;
%BGCp(3) = 0.25*TotalOrgMass;
%BGCp(3) = 0;
BGCp(3) = TotalOrgMass;
%BGCp(4) = TotalOrgMass-BGCp(3);
BGCp(4) = 0;


for i=1:4
   x(i) = BGCp(i)/sum(BGCp);
   BGCv(i) = x(i)*CStarBasis(i);
  % BGCv(i) = 0;
end
  


figure(1)
semilogx(Dp,n_N)
hold on
xlabel('Diameter (nm)')

if Tag==1
    ylabel('dX/dlogDp')
end

%hold off
%figure(11)
semilogx(Dp,n_S)
semilogx(Dp,n_V)
hold off

%Total2 = sum(TotalCalc)
%figure(2)
%semilogx(Dp,n_S)
%hold on
%semilogx(Dp,n_N)
