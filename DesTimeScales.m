
global modelAtm


C_Star = 1e-2; %ug/m3
Uptake = 1;

nm = zeros(1,1);
m = zeros(1,1);
tau = zeros(1,1);

C_Star = [1e-4 1e-3 1e-2 1e-1];
End = length(C_Star);
for j = 1:End
   Tau_cut = 3*3600; %seconds
   Dp_cut(j) = Tau_cut/(modelAtm.SOA.rho)*(3*Uptake*C_Star(j)*modelAtm.SOA.va)*2
end

j = 6;
%loop through diameters

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on');
box('on');
hold('all');

nm = 10.^(0:0.02:3);
%figure(688)
for j = 1:End
for i = 1:length(nm)
   
   %nm(i) = i;
   m(i) = nm(i)*1e-9;
   Ke(i) = exp(4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/...
    (8.314*298*modelAtm.SOA.rho*1e-6*m(i)));
   KeFactor(i) = Ke(i)^(-1);
   %KeFactor(i) = 1;
   tau(i) = KeFactor(i)*modelAtm.SOA.rho*m(i)/2/(3*Uptake*modelAtm.SOA.va*C_Star(j));
 %  tau(i) = 700*m(i)*1e6/C_Star(j); %Approximate as done in basis set gloss
   tau_hr(i) = tau(i)*1/3600;
    
end

ymatrix(j,:) = tau_hr;
%loglog1 = loglog(nm,tau,'Parent',axes1)
%hold on

%set(loglog1(j),'DisplayName',num2str(C_Star(j)));

end

ymatrix(j+1,1:length(tau_hr))=1;
loglog1 = loglog(nm,ymatrix,'Parent',axes1);
ylabel('Desorption Timescale (hr)')
xlabel('Particle size (nm)')

set(loglog1(1),'DisplayName','lC* = -4')
set(loglog1(2),'DisplayName','lC* = -3')
set(loglog1(3),'DisplayName','lC* = -2')
set(loglog1(4),'DisplayName','lC* = -1')
set(loglog1(5),'DisplayName','tau = 1 hr','LineStyle',':','Color',[0 0 0]);
%set(loglog1(6),'DisplayName','lC* = 3')
%set(loglog1(7),'DisplayName','lC* = 4')

legend(axes1,'show');


