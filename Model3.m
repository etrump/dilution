function dY=Model3(t,Y)
dY=zeros(2,1);

%load globals -----------------------
global modelAtm trackOmega trackCov trackCondSink trackMode

n = modelAtm.NumBins;
Pop = modelAtm.Pop;
%CStarBasis = modelAtm.CStarBasis;
alpha = modelAtm.SOA.alpha;

if modelAtm.PulseYN==0
    Tau = modelAtm.DecayConstant*3600;
    ROG = modelAtm.Injection/Tau*exp(-t/Tau);
else
    if t>modelAtm.EmitTime
        ROG = 0;
        modelAtm.Sulf.MassConc = eps;
        %xmodelAtm.Sulf.MassConc = 0;
    else
        ROG = modelAtm.ROG;
    end
end


%Load current values of integration variables -----------
for i=1:n
    Cv(i)=Y(i);  

%    for j=1:Pop
%        Cp(j,i)=Y(j*n+i);
%    end
    if Cv(i)<=0 || Cv(i)~=Cv(i) 
        Cv(i) = 0; % Borrowed from Evaporation mfile
    end
    for j=1:Pop
        Cp(j,i)=Y(j*n+i);
        if Cp(j,i)<=0 || Cp(i) ~= Cp(i)
            Cp(j,i) = 0;% Borrowed from Evaporation mfile
        end
    end
end

for j=1:Pop
    Dp(j) = Y((Pop+1)*n+j)*1e-9; %  *1e-9 (nanometers to meters)
    M_Sulf(j) = Y((Pop+1)*n+Pop+j);
    M_SOA(j) = sum(Cp(j,:));
    NumConc(j) = Y((Pop+1)*n+2*Pop+j);
  %  NumConc(j) = modelAtm.Pop1.NumConc0;
 %   Dp(j) =(M_SOA(j)/(NumConc(j)*modelAtm.SOA.rho+eps^2)*6/pi)^(1/3); % calculation diameter based on mass
end

N_max = 0;
j_max = 1;
for j = 1:Pop
    if NumConc(j) > N_max
        N_max = NumConc(j);
        j_max = j;
    end
end

   trackMode(length(trackMode(:,1))+1,1) = t;
   trackMode(length(trackMode(:,1)),2) = Dp(j_max);
   trackMode(length(trackMode(:,1)),3) = j_max;
   trackMode(length(trackMode(:,1)),4) = NumConc(j_max);
   

% NumConc(1) = modelAtm.Pop1.NumConc0;
% NumConc(2) = modelAtm.Pop2.NumConc0;
% NumConc(3) = modelAtm.Pop3.NumConc0;
% NumConc(4) = modelAtm.Pop4.NumConc0;
% NumConc(5) = modelAtm.Pop5.NumConc0;
% NumConc(6) = modelAtm.Pop6.NumConc0;
% NumConc(7) = modelAtm.Pop7.NumConc0;
% NumConc(8) = modelAtm.Pop8.NumConc0;
% NumConc(9) = modelAtm.Pop9.NumConc0;
% NumConc(10) = modelAtm.Pop10.NumConc0;
% if modelAtm.Pop == 20
% NumConc(11) = modelAtm.Pop11.NumConc0;
% NumConc(12) = modelAtm.Pop12.NumConc0;
% NumConc(13) = modelAtm.Pop13.NumConc0;
% NumConc(14) = modelAtm.Pop14.NumConc0;
% NumConc(15) = modelAtm.Pop15.NumConc0;
% NumConc(16) = modelAtm.Pop16.NumConc0;
% NumConc(17) = modelAtm.Pop17.NumConc0;
% NumConc(18) = modelAtm.Pop18.NumConc0;
% NumConc(19) = modelAtm.Pop19.NumConc0;
% NumConc(20) = modelAtm.Pop20.NumConc0;
% end
%
FracSulfSusp = M_Sulf(1)/(M_Sulf(1)+M_SOA(1));

%for j=1:Pop
   % PopString = int2str(j);
    %eval(['modelAtm.Pop' PopString '.Dp = Dp(j);']) % do this need to be
    %done???
   % eval(['NumConc(j) = modelAtm.Pop' PopString '.NumConc0;']); %do this
   % explicitly for each population to save time
%end


%Load associated properties ----------

for j=1:Pop
    TotalMass(j) = M_Sulf(j)+M_SOA(j); %Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
    rho = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.rho+M_SOA(j)/modelAtm.SOA.rho);
    MW = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.MW+M_SOA(j)/modelAtm.SOA.MW);
end

%for i=1:n % is this used???
  
%   XCp(1,i) = Cp(1,i)/sum(Cp(1,:));
%   XCv(i) = Cv(i)/sum(Cv);
%   if sum(Cp(1,:)) == 0
%       XCp(1,i) = 1;
%       XCv(i) = 1;
%   end
   
%end

t

%Calculate condensation sinks------------

for j=1:Pop
    ParticleKn = 2*modelAtm.SOA.lambda/Dp(j);
    BetaP = FuchsC(ParticleKn);
    CondSinkTot(j) = alpha*NumConc(j)*2*pi*Dp(j)*modelAtm.SOA.Diff*BetaP;
   % CondSinkTot(j) = UpdateBackgroundCS(TotalMass(j),Dp(j),rho, j);

   % PopString = int2str(j);
  %  eval(['ParticleNum = modelAtm.Pop' PopString '.NumConc0;']) %No deposition or particle loss
    %ParticleNum = 1; %#/m3
    ParticleNum = NumConc(j);
    CondSinkP(j) = CondSinkTot(j)/ParticleNum;
end

if modelAtm.CoagYN == 1
    [dNCoag, dCpCoag, dSulfCoag] = ParticleCoag(NumConc,Cp,Dp,M_Sulf);
end

%Calculates organic flux to bulk particles---------------------------
for j=1:Pop
  %  trackCov(length(trackCov(:,2)),1) = t;
    [dCv dCp] = ParticleGrowth(Cv,Cp(j,:),CondSinkTot(j),j,modelAtm.KelvinYN,ROG,Dp(j)); 
    
    dCvStack(j,:)=dCv;
    
    for i=1:n
        if NumConc(j)>=1
            if modelAtm.CoagYN==1
                dY(j*n+i) = dCp(i)+dCpCoag(j,i);
            else
                dY(j*n+i) = dCp(i);
            end
                
        else
            dY(j*n+i) = 0;
        end
    end
    J_SOA(j) = sum(dCp);
    J_SOA_Cond(j) = sum(dCp);
end
   
for i=1:n
    dY(i) = sum(dCvStack(:,i));
end

     
%Calculate Sulfate growth on single particle and bulk growth-----------------------------------------     
  
for j=1:Pop
    J_Sulf_Cond(j) = CondSinkTot(j)*modelAtm.Sulf.MassConc;
    J_Sulf(j) = J_Sulf_Cond(j); %[=]ug/m3-s S&P eq12.12 p539
    if modelAtm.CoagYN==1 
        dY((Pop+1)*n+Pop+j) = J_Sulf(j)+dSulfCoag(j);
    else
        dY((Pop+1)*n+Pop+j) = J_Sulf(j);
    end

    rho_T(j) = (J_Sulf_Cond(j)+J_SOA_Cond(j))/(J_Sulf_Cond(j)/modelAtm.Sulf.rho+J_SOA_Cond(j)/modelAtm.SOA.rho);
    rho(j) = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.rho+M_SOA(j)/modelAtm.SOA.rho);
    rho_T(j) = rho(j);%%%%%%%%%%%%%% assumes density of condensing is same as density of particle itself
    
    J_T(j) = CondSinkP(j)/CondSinkTot(j)*(J_Sulf(j)+J_SOA(j));
    J_T_Cond(j) = CondSinkP(j)/CondSinkTot(j)*(J_Sulf_Cond(j)+J_SOA_Cond(j));
    
    dDp0(j) = (2*J_T_Cond(j))/(rho_T(j)*pi*Dp(j)^2); %growth rate of particles in population at start of timestep
    
    
    M_All(j) = M_SOA(j) + M_Sulf(j);
    J_All(j) = J_SOA(j) + J_Sulf(j);
    
    dY((Pop+1)*n+j)=dDp0(j)*1e9; %meters to nanometers
  %  dY((Pop+1)*n+j)=0; % calculation diameter based on mass
    if modelAtm.CoagYN==1 
        dY((Pop+1)*n+2*Pop+j) = dNCoag(j);
    else
        dY((Pop+1)*n+2*Pop+j) = 0;
    end
      
end 



end    
    
