function [dNCoag, dCpCoag, dSulfCoag] = ParticleCoag(NumConc,Cp,Dp,M_Sulf)

global modelAtm

fracorg = 1;
Dp;
nbins = modelAtm.Pop;
temp = modelAtm.Temp;
pres = modelAtm.Press; % check units
dia = Dp;
coag_eff = 1;

bin_mean_mass = modelAtm.SOA.rho*pi/6*Dp.^3;
num_in_bin = NumConc;

nvolbins = length(modelAtm.CStarBasis);
dCpCoag = zeros(nbins,nvolbins);

%% calculate the coagulation kernel
gc = 8.314; % Gas constant in J/mol*K
kb = gc/6.022e23; % boltzman constant
mu=2.5277e-7*temp^0.75302;  % viscosity in kg m-1 s-1
mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gc*temp)));  % mean free path of air S&P eqn 8.6 in m

for k=1:nbins
   ck(k) = (8*kb*temp/(pi*bin_mean_mass(k)*1e-9))^(1/2); % coefficient from s&p table 12.1
   Kn = 2*mfp/dia(k); % Knudson number s&p table 12.1
   Dk(k)=kb*temp/(3.0*pi*mu*dia(k))*((5.0+4.0*Kn+6.0*Kn^2+18.0*Kn^3)/(5.0-Kn+(8.0+pi)*Kn^2));% diffusivity table 12.1 s&p in m2/s
end

% calulate kernel, but don't repeat calculations
for i=1:nbins
   for j=i:nbins
      Kn=4.0*(Dk(i)+Dk(j))/(sqrt(ck(i)^2+ck(j)^2)*(dia(i)+dia(j)));  %Knudson number for 2 particles S&P eqn 12.51
      beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn));          %correction factor for non-continuum S&P eqn 12.50
      %coag_kernel(i,j)=2.0*pi*(dia(i)+dia(j))*(Dk(i)+Dk(j))*beta*1e6;  % coagulation coefficient in cm3/ps  S&P eqn 12.49
      coag_kernel(i,j)=2.0*pi*(dia(i)+dia(j))*(Dk(i)+Dk(j))*beta; % m3/s
      %coag_kernel(i,j)=1e-2;
   end
   
end
coag_kernel;
% fill rest of symetrical matrix
for j=1:nbins
    for i=1:j-1
        coag_kernel(j,i)=coag_kernel(i,j);
    end
end

coag_kernel;
nbins;
%% Apply coag
coag_ch(1:nbins) = 0;
for i=1:nbins
    for j=i:nbins
        coag_rate = coag_eff*coag_kernel(i,j)*num_in_bin(i)*num_in_bin(j); % the number of coagulation events between i and j per second per cm3
        coag_ch(i) = coag_ch(i)-coag_rate;
        coag_ch(j) = coag_ch(j)-coag_rate;
        coag_mass=bin_mean_mass(i)+bin_mean_mass(j);
        bin_mean_mass;
        % find which bin this goes in
        k=j;
        while ((k < nbins) & (bin_mean_mass(k+1) < coag_mass))
            k=k+1;
        end
        if (k < nbins)
            fc=(coag_mass-bin_mean_mass(k))/(bin_mean_mass(k+1)-bin_mean_mass(k));
            coag_ch(k)=coag_ch(k)+coag_rate*(1-fc); % add one particle per cm3 from bin per coagulation event
            coag_ch(k+1)=coag_ch(k+1)+coag_rate*fc;
        elseif (k==nbins)
            %fc=(coag_mass-bin_mean_mass(k))/(bin_mean_mass(k+1)-bin_mean_mass(k));
            %fc=(coag_mass-bin_mean_mass(k-1))/(bin_mean_mass(k)-bin_mean_mass(k-1)); %?????
            coag_ch(k)=coag_ch(k)+coag_rate;
        end
    end
    
    
end

coag_ch;

for j = 1:nbins
    for i = 1:nvolbins
        molfrac(j,i) =  Cp(j,i)/sum(Cp(j,:));
        if sum(Cp(j,:)) == 0
            molfrac(j,i) = 0;
        end  
        dCpCoag(j,i) =  bin_mean_mass(j)*molfrac(j,i)*coag_ch(j)*fracorg; %CHECK THIS           
    end
    dSulfCoag(j) = M_Sulf(j)*coag_ch(j)*bin_mean_mass(j)*(1-fracorg);
    
    dNCoag(j) = coag_ch(j);
end