function [omega]=WallFluxCorrection(j)

% j is the population of interest
global modelAtm trackOmega trackCov
Pop = modelAtm.Pop;

PopString = int2str(j);



%eval(['lambda = modelAtm.Pop' PopString '.lambda;']);
lambda = modelAtm.SOA.lambda;
alpha = modelAtm.SOA.alpha;
%eval(['alpha = modelAtm.Pop' PopString '.alpha']); % why doesn't this
%work???
%alpha = 1;
eval(['Diff = modelAtm.Pop' PopString '.Diff;']);
eval(['shapeFact = modelAtm.Pop' PopString '.shapeFact;']);
eval(['va = modelAtm.Pop' PopString '.va;']);

eval(['Dp = modelAtm.Pop' PopString '.Dp;']); %m
eval(['NumConc = modelAtm.Pop' PopString '.NumConc;']);

if j>1
    eval(['coverage = modelAtm.Pop' PopString '.coverage;']);
    coverage = 0;
end


Kn =  2*lambda/Dp; %lambda and Dp in [m]
alpha_cr = 8*sqrt(modelAtm.EddyDiffCoeff*Diff)/(pi*va); %all SI units
modelAtm.alpha_cr = alpha_cr;

% 
% total_cov = 0;
% for i=2:Pop
%    CovString = int2str(i);
%    eval(['total_cov = total_cov + modelAtm.Pop' CovString '.coverage;'])
% end

total_cov = modelAtm.total_cov;


%total_covIS = total_cov;

n_a = NumConc/modelAtm.StoV; % #/m2

if (total_cov>modelAtm.Area)
    max = 100000
else
    coverage = n_a*pi*Dp^2*shapeFact;
end

if j<2
    coverage = 0;
end


eval(['modelAtm.Pop' PopString '.n_a = n_a;']);
eval(['modelAtm.Pop' PopString '.coverage = coverage;']);


if shapeFact ==1
    KnFactor = 1;
else
    KnFactor = 2*shapeFact; %for hemisphere
end
Factor = FuchsK(Kn*KnFactor);


uptakecoeffOld = alpha*coverage*FuchsK(Kn*KnFactor);
uptakecoeff = alpha*total_cov*FuchsK(Kn*KnFactor);
%uptakecoeff = uptakecoeffOld;


satcurvOld = (1+uptakecoeffOld/alpha_cr)^(-1);
satcurv = (1+uptakecoeff/alpha_cr)^(-1);
%j;
omega_omegaOld = satcurv/satcurvOld;
omega = FuchsK(KnFactor*Kn)/FuchsK(Kn)*shapeFact*satcurv;
%omega = 1;

%if j==3
%    omega = omega*0.72
% %end
% 
% if omega>0
%     omega;
% else
%     omega=1;
% end


%IF not on wall 
if(eval(['modelAtm.Pop' PopString '.wallYN == 0;']))
    omega = 1;
end


%Index = length(trackOmega(:,1))+1;
%trackOmega(Index,j) = omega;


%Index = length(trackCov(:,1))+1;
%%%%trackCov(Index,1) = t;
%trackCov(Index,2) = j;
%trackCov(Index,3) = total_cov;
%trackCov(Index,4) = coverage;

    
    


