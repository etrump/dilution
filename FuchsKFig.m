

global modelAtm

Dp_nm = 1:1e3;
%Dp_nm = 10.^(0:0.02:3)
alpha = 1;

%modelAtm.SOA.lambda2 = 6.6e-8;

for i = 1:length(Dp_nm)
   Dp_m(i) = Dp_nm(i)*1e-9;
   Kn(i) = 2*modelAtm.SOA.lambda/Dp_m(i);
   FuchsK(i) = (Kn(i)+Kn(i)^2)/(Kn(i)^2+Kn(i)+0.283*Kn(i)*alpha+0.75*alpha);
end

figure(888)
plot(Dp_nm,FuchsK)

xlabel('Particle Diameter (nm)');
ylabel('Fuchs corretion (J/J_k)');
title('Fuchs correction');