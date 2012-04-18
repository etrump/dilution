function [Caer Cvap] = Partition(C_reacted)
%this function calculates partitioning.  C is the total gas + aerosol conc
%of species i in ug/m^3, Cstar is the saturation concentrations, Coa is a
%guess for the total organic aerosol concentration, maxiter is the maximum
%number of iterations, and tol is the solution tolerance in Coa.

%C_reacted = 1.66e3/0.1;

global modelAtm


Cstar = modelAtm.CStarBasis
C = C_reacted*modelAtm.SOA.alphaProd;
maxiter = 1e5;
tol = 1e-15;
Coa = 1; %initial guess


for n = 1:maxiter
    xi = (1 + Cstar/Coa ) .^ (-1);
    CoaNEW = sum(C .* xi);
    if abs(CoaNEW - Coa) < tol
        break
    else
        Coa = CoaNEW;
    end
end

if n == 1000
    error('Did not Converge')
end

n

Caer = C * CoaNEW ./ Cstar .* (1 + CoaNEW ./ Cstar) .^ (-1)
Cvap = C - Caer

sum(Caer)
