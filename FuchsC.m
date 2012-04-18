function Beta = FuchsC(Kn,alpha)  
%Calculates the Fuchs&Sutugin transition correction correction factor,
%this correction is for the continuum regime flux - given Kn, alpha

if nargin < 2
    alpha = 1;
end
    
Beta = 0.75*alpha*(1+Kn)/(Kn^2+Kn+0.283*Kn*alpha+0.75*alpha);


