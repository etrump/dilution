function ProcessData()

global modelAtm trackCov MassforPlot MassforPlotAlpha Struct
global Tagg CpT DpT NumConcT CvT DpT2
global trackCondSink

Yagg = Struct.Yagg;
Tagg = Struct.Tagg;

n = modelAtm.NumBins;
TotalPop = modelAtm.Pop;

%trackCondSink = zeros(1,11);

%modelAtm.Dpf = Dpf;
EndIndex(1) = 0;

CvT = Yagg(:,(1:n));
CpT = Yagg(:,n+1:(TotalPop+1)*n);
DpT = Yagg(:,n*(TotalPop+1)+1:n*(TotalPop+1)+TotalPop);
MSulfT = Yagg(:,n*(TotalPop+1)+TotalPop+1:n*(TotalPop+1)+2*TotalPop);
NumConcT = Yagg(:,n*(TotalPop+1)+TotalPop*2+1:n*(TotalPop+1)+3*TotalPop);


%Check mass balance---------------------------------------------------------------------------
TotalSuspMSOA = 0;
TotalPop2MSOA = 0;
TotalVap = 0;
%
%for j=1:n
%    TotalPop2MSOA = TotalPop2MSOA + CpT(length(CpT),n+j);
%    TotalSuspMSOA = TotalSuspMSOA + CpT(length(CpT),j);
%    TotalVap = TotalVap + CvT(length(CpT),j);
%end

TotalMSOA0 = sum(modelAtm.Pop1.Cp0)+sum(modelAtm.Pop2.Cp0);
Tau1 = modelAtm.DecayConstant*3600;

if modelAtm.PulseYN==1
    CheckMassBal = sum(modelAtm.Cv0)+modelAtm.ROG*modelAtm.EmitTime*sum(modelAtm.SOA.alphaProd)+TotalMSOA0-TotalPop2MSOA-TotalSuspMSOA-TotalVap;
else
    CheckMassBal = sum(modelAtm.Cv0)+modelAtm.Injection*(1-exp(-FinalTime/Tau1))*sum(modelAtm.SOA.alphaProd)+TotalMSOA0-TotalPop2MSOA-TotalSuspMSOA-TotalVap;
end

toler = 1e-7;
if abs(CheckMassBal) > 1e-7
    Error = 'Mass balance is not correct!';
else
    Error = 'NO mass balance error!';
end

%%CheckEqm(1) = (CpT(length(CpT),1)+CpT(length(CpT),n+1))/(CvT(length(CpT),1)+CpT(length(CpT),1)+CpT(length(CpT),1*n+1))-...
%    1/(1+modelAtm.CStarBasis(1)/(TotalPop2MSOA+TotalSuspMSOA));
%CheckEqm(2) = NaN;
%CheckEqm(3) = (CpT(length(CpT),3)+CpT(length(CpT),n+3))/(CvT(length(CpT),3)+CpT(length(CpT),3)+CpT(length(CpT),1*n+3))-...
%    1/(1+modelAtm.CStarBasis(3)/(TotalPop2MSOA+TotalSuspMSOA));


CheckEqm(1) = (CpT(length(Tagg),1)+CpT(length(Tagg),n+1))/(CvT(length(Tagg),1)+CpT(length(Tagg),1)+CpT(length(Tagg),1*n+1))-...
    1/(1+modelAtm.CStarBasis(1)/(TotalPop2MSOA+TotalSuspMSOA));
CheckEqm(2) = NaN;
%CheckEqm(3) = (CpT(length(Tagg),3)+CpT(length(Tagg),n+3))/(CvT(length(Tagg),3)+CpT(length(Tagg),3)+CpT(length(Tagg),1*n+3))-...
%    1/(1+modelAtm.CStarBasis(3)/(TotalPop2MSOA+TotalSuspMSOA));

FracMassToSuspended = TotalSuspMSOA/(TotalSuspMSOA+TotalPop2MSOA);

%Calculate mole fractions of organics on particles--------------------------
%for p = 1:TotalPop
%   PopString = int2str(p);
%   XCpT_1 = CpT(:,(p-1)*n+3)./(CpT(:,(p-1)*n+1)+CpT(:,(p-1)*n+3));
%   XCpT_01 = CpT(:,(p-1)*n+1)./(CpT(:,(p-1)*n+1)+CpT(:,(p-1)*n+3));
  % eval(['Struct.P' PopString '.XCpT_1 = XCpT_1;']);
%   eval(['modelAtm.Pop' PopString '.XCpT_1 = XCpT_1;']);
  % eval(['Struct.P' PopString '.XCpT_1 = XCpT_1;']);
%   eval(['modelAtm.Pop' PopString '.XCpT_1 = XCpT_1;']);
%end

%  for i=1:length(Tagg)
%      Cp1Wall(i) = 0;
%      Cp01Wall(i) = 0;
%      Cp_1Wall(i) = 0;
%      Cp1_Wall(i) = 0;
%      for j = 2:TotalPop
%         Cp1Wall(i) = Cp1Wall(i) + CpT(i,4*(j-1)+3);
%         Cp_1Wall(i) = Cp_1Wall(i) + CpT(i,4*(j-1)+2);
%         Cp01Wall(i) = Cp01Wall(i) + CpT(i,4*(j-1)+1);
%         Cp1_Wall(i) = Cp1_Wall(i) + CpT(i,4*(j-1)+4);
%      end
%      MSulfWall(i) = MSulfT(i,2);
%      Cp1Sus(i) = CpT(i,3);
%      Cp1_Sus(i) = CpT(i,4);
%      Cp_1Sus(i) = CpT(i,2);
%      Cp01Sus(i) = CpT(i,1);
%      MSulfSus(i) = MSulfT(i,1);
%  end
% 
% MassforPlot = [Cp01Wall; Cp1Wall; Cp01Sus; Cp1Sus; MSulfWall; MSulfSus];
%MassforPlotAlpha = [Cp01Wall; Cp_1Wall; Cp1Wall; Cp1_Wall; Cp01Sus; Cp_1Sus; Cp1Sus; Cp1_Sus; MSulfWall; MSulfSus];

%MassforPlotAlpha = [0 0 0 0];
%%for j=1:TotalPop
 %   for i = 1:n
 %       MassforPlotAlpha(i) = CpT(
 %   end
%end



% for j = 1:length(Tagg)
% %    TotalOrgMass = sum(CpT(j,1:modelAtm.NumBins));
% %    TotalMass(j) = TotalOrgMass+MSulfT(j,1);
% %    rho = (TotalOrgMass+MSulfT(j,1))/(MSulfT(j,1)/modelAtm.Sulf.rho+TotalOrgMass/modelAtm.SOA.rho);
% %    CondSinkBG(j) = UpdateBackgroundCS(TotalMass(j),DpT(j,1),rho);
% %    BG_MtoCS(j) = TotalMass(j)/CondSinkBG(j);
% %    
% %   
% %    TotalOrgMass = sum(CpT(j,5:2*modelAtm.NumBins));
% %    TotalMass(j) = TotalOrgMass + MSulfT(j,2);
% %    rho = (TotalOrgMass+MSulfT(j,2))/(MSulfT(j,2)/modelAtm.Sulf.rho+TotalOrgMass/modelAtm.SOA.rho);
% %    
% %    ParticleKn = 2*modelAtm.SOA.lambda/DpT(j,2);
% %    BetaP = FuchsC(ParticleKn);
% %    CondSinkNanoP(j) = modelAtm.SOA.alpha*modelAtm.Pop2.NumConc0*2*pi*DpT(j,2)*modelAtm.SOA.Diff*BetaP;
% %    %CondSinkNanoP(j) = UpdateBackgroundCS(TotalMass(j),DpT(j,2),rho);
% %    NanoP_MtoCS(j) = TotalMass(j)/CondSinkNanoP(j);
%    
%    
%    
%    
%    
% end



%%DISPLAY FIGURES----------------------------------------------------


%%CreateFig;
%%CreateFigMassAlpha;
