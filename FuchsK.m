function Beta = FuchsK(Kn,alpha)  
%Calculates the Fuchs&Sutugin transition correction correction factor,
%this correction is for the kinetic regime flux - given Kn, alpha


if nargin < 2
    alpha = 1;
end
    
Beta = (Kn+Kn^2)/(Kn^2+Kn+0.283*Kn*alpha+0.75*alpha);


