global modelAtm

modelAtm.Pop88.shapeFact = 1;
nm = zeros(1,1)*NaN;
Ke = zeros(1,1)*NaN;

End = 1e3;

nm = 10.^(0:0.02:3);

for i = 1:length(nm)
    
    m(i) = nm(i)*1e-9;
    modelAtm.Pop88.Dp = m(i);
    
    Ke(i) = KelvinTerm(88);
end

figure(556)
loglog(nm,Ke)
%hold on
%semilogx(nm,1)

xlabel('Particle Diameter (nm)');
ylabel('Kelvin Term');
title('Kelvin');
%axis([0 End 0 6]); 

Ke(7)
Ke(113)