function ans = SlipCorr(Dp,Lambda)

%Make sure to pass Dp, Lambda in equivalent units!
%S&P - eq 9.34 (p407)

ans = 1+2*Lambda/Dp*(1.257+0.4*exp(-1.1*Dp/(2*Lambda)));
%ans = 5 % FOR TESTING