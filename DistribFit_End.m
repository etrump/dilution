global modelAtm

xData = DpT(length(Tagg),:)'
yData = TotalSusp'

N_tot = sum(TotalSusp);

ft = fittype( @(a,b,x) N_tot/((2*pi())^(1/2)*a)*exp(-(log10(x)-b).^2/(2*a.^2)));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.StartPoint = [0.823457828327293 0.694828622975817];
opts.Upper = [Inf Inf];


% Fit model to data.
[fitresult, gof, O] = fit( xData, yData, ft, opts );

coeffval = coeffvalues(fitresult)
logSigma = coeffval(1)


logMode = coeffval(2)