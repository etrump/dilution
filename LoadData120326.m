function [Dp_bin_vect, N_tot_bin] = LoadData120326(TotalPop) %alpha-P 120326 exper
%TotalPop = 100;
global modelAtm
TotalPop = modelAtm.TotalPop;

A1 = [12.2
12.6
13.1
13.6
14.1
14.6
15.1
15.7
16.3
16.8
17.5
18.1
18.8
19.5
20.2
20.9
21.7
22.5
23.3
24.1
25
25.9
26.9
27.9
28.9
30
31.1
32.2
33.4
34.6
35.9
37.2
38.5
40
41.4
42.9
44.5
46.1
47.8
49.6
51.4
53.3
55.2
57.3
59.4
61.5
63.8
66.1
68.5
71
73.7
76.4
79.1
82
85.1
88.2
91.4
94.7
98.2
101.8
105.5
109.4
113.4
117.6
121.9
126.3
131
135.8
140.7
145.9
151.2
156.8
162.5
168.5
174.7
181.1
187.7
194.6
201.7
209.1
216.7
224.7
232.9
241.4
250.3
259.5
269
278.8
289
299.6
310.6
322
333.8
346
358.7
371.8
385.4
399.5
414.2
429.4
445.1
461.4
478.3
495.8
514
532.8
552.3
572.5];
Diam = A1'; %nm;

B1 = [48485.7
31639.2
20043.9
14133.2
10194.8
6214.83
3284.95
1292.58
774.588
440.525
272.931
112.281
137.117
83.6042
79.2988
46.1827
43.1922
26.9899
0
70.9468
22.3869
10.6139
0
19.0474
36.3694
25.8614
37.3618
3.7877
22.6739
14.4626
7.05247
14.5061
11.8724
12.4686
11.9764
18.1615
38.686
58.6312
128.3
288.382
403.236
593.842
1035.68
1327
2153.22
2946.98
4739.52
6620.21
9544.2
13943.5
18756.6
25898.3
34124.6
45581.7
58923.4
74559.3
93626.3
117493
145573
175930
209570
248710
291936
335737
380860
430633
482591
527309
575768
620469
655738
692755
723952
751013
777789
794498
812892
815078
821768
822654
822213
816575
808274
800328
784421
767029
750952
7.25E+05
6.95E+05
6.65E+05
6.21E+05
5.76E+05
524213
470858
413773
352764
295286
241432
190702
145848
108762
78084.6
54961.8
36160.1
23119.1
14191.7
8117.36
4824.12];
dNdlogDp_0 = B1';

rho_SOA = 1.4e12; %ug/m3
for i = 1:length(Diam)-1
   N(i) =  (log10(Diam(i+1))-log10(Diam(i)))*(dNdlogDp_0(i)+dNdlogDp_0(i+1))*1/2;
   Mass_contrib(i) = N(i)*pi()/6*(1/2*1e-9*(Diam(i+1)+Diam(i)))^3*rho_SOA;
end
N(length(Diam)) = N(length(Diam)-1);

N_tot = sum(N)
Mass_tot = sum(Mass_contrib)*1e6 % ug/m3

%% Fit: 'Guassian'.
%[xData, yData] = prepareCurveData( Diam, dNdlogDp_0 );
xData = Diam';
yData = dNdlogDp_0';

% Set up fittype and options.
%ft = fittype( '4.7597e+05/((2*pi())^(1/2)*a)*exp(-(log10(x)-b)^2/(2*a^2))', 'independent', 'x', 'dependent', 'y' );
ft = fittype( @(a,b,x) N_tot/((2*pi())^(1/2)*a)*exp(-(log10(x)-b).^2/(2*a.^2)));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.StartPoint = [0.823457828327293 0.694828622975817];
opts.Upper = [Inf Inf];


% Fit model to data.
[fitresult, gof, O] = fit( xData, yData, ft, opts );

coeffval = coeffvalues(fitresult)

figure(55)
% Plot fit with data.
figure( 'Name', 'Guassian' );
h = plot( fitresult, xData, yData );
%h = plot(xData, yData);
legend( h, 'dNdlogDp_0 vs. Diam', 'Guassian', 'Location', 'NorthEast' );
% Label axes
xlabel( 'Diam' );
ylabel( 'dNdlogDp_0' );
grid on
%hold on
%figure(33)
%plot(Diam,dNdlogDp_0)
%hold off

%% check fit

%logSigma = 0.1597
%logSigma = 0.1847 %small bag, t0
logSigma = coeffval(1)


%logMode = 2.218
%logMode = 2.447 %Small bag, t0
logMode = coeffval(2)

Dp_0m = 50; %[=]nm, start at 0.010 um
Dp_fm = 550;   % match SMPS

logDp_0m = log10(Dp_0m);
logDp_fm = log10(Dp_fm);
logDpm = logDp_0m:0.001:logDp_fm;
 for i = 1:length(logDpm)
     Dpm(i) = 10^(logDpm(i));
 end
 n_N_bin = zeros(5,length(Dpm));
for i = 1:length(Dpm)
   % if Dp(i)==Mode
    d(i) = N_tot/((2*pi)^(1/2)*logSigma)*exp(-(logDpm(i)-logMode)^2/(2*(logSigma)^2)); %dNdlogDp, fitted
    %d(i) = Total; %S&P pg 370, eq 8.54
    %Beta(i) = 1;

    n_N(i) = d(i);  %#/m3 
end
figure(99)
plot( fitresult, xData, yData);
hold on
plot(Dpm,n_N)
hold off

figure(100)
semilogx(Dpm,n_N)

%TotalPop = 20;
xn = 550;
x0 = 50;
delh = (log10(xn)-log10(x0))/TotalPop;
crap = 1:TotalPop;
Bin_split = 10.^(delh*crap + log10(x0));
Dp_bin_vect(1) = (x0+(Bin_split(1)))/2;
for i = 2:TotalPop
    Dp_bin_vect(i) = (Bin_split(i)+Bin_split(i-1))/2;
end

Dp_bin_vect;

N_bin = zeros(10,length(Dpm)-1);

for i = 1:length(Dpm)-1
   N_check(i) =  (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
   M_check(i) =  N_check(i)*pi/6*modelAtm.SOA.rho*(1/2*1e-9*(Dpm(i+1)+Dpm(i)))^3; 
   dump = 0;
   for j=1:TotalPop
       if Dpm(i) < Bin_split(j) && dump~=1
           N_bin(j,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
           dump = 1;
       end
   end
   if dump == 0
       N_bin(1,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
       dump = 1;
   end
           
end
%     if Dpm(i) < Bin_split(1)
%         N_bin(1,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(2)
%         N_bin(2,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(3)
%         N_bin(3,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(4)
%         N_bin(4,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2; 
%     elseif Dpm(i) < Bin_split(5)
%         N_bin(5,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(6)
%         N_bin(6,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(7)
%         N_bin(7,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(8)
%         N_bin(8,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%     elseif Dpm(i) < Bin_split(9)
%         N_bin(9,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(10)
%        N_bin(10,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(11)
%        N_bin(11,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(12)
%        N_bin(12,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2; 
%    elseif Dpm(i) < Bin_split(13)
%        N_bin(13,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(14)
%        N_bin(14,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(15)
%        N_bin(15,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(16)
%        N_bin(16,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(17)
%        N_bin(17,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(18)
%        N_bin(18,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    elseif Dpm(i) < Bin_split(19)
%        N_bin(19,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    else
%        N_bin(20,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
%    end
%end

for j=1:TotalPop
    N_tot_bin(j) = sum(N_bin(j,:));
    M_tot_bin(j) = N_tot_bin(j)*pi/6*(Dp_bin_vect(j)*1e-9)^3*modelAtm.SOA.rho;
end

for j = 1:TotalPop-1
    nN_bin(j) = 1/2*(N_tot_bin(j+1)+N_tot_bin(j))*(log10(Dp_bin_vect(j+1))-log10(Dp_bin_vect(j))); 
end
nN_bin(TotalPop) = 0;

N_tot_bin;
N_tot;
sum(N_check);
M_check_result = sum(M_check) * 1e6
M_check;
M_bin_check = sum(M_tot_bin) * 1e6 % ug /m3

%Binned_tot = sum(N_tot_bin)
Unbinned_tot = sum(N_check)
Binned_tot = sum(N_tot_bin)

%hold on
%%plot(Dpm,n_N)

figure(22)
%plot(Dp_bin_vect,N_tot_bin,Diam,N)

%plot(Dp_bin_vect,nN_bin,Dpm,n_N)
semilogx(Dp_bin_vect,N_tot_bin/max(N_tot_bin),Diam(1:length(N)-1),N(1:length(N)-1)/max(N))

%% Estimate mass
Dp_mode = 10^(logMode)*1e-9; % 2.218 is from lognormal Fit for logDp_mode and 1e-9 converts to meters
Dp_mode_nm = 10^(logMode)

rho_SOA = 1.4e12; %ug/m3
%Mass_est = N_tot*pi()/6*Dp_mode^3*rho_SOA * 1e6 %ug/m3 *very rough
%estimate*


%% hygroscopic growth

%A_factor = 4*18*0.07564/(8.314*300*1e6);
%B_factor = 6*18/(pi()*1e6);
%
%RH = 0.75;
%n_s = Dp_mode^3/B_factor *(A_factor/Dp_mode-log(RH)); % mol
%m_s = n_s*200*1/2; % g
%
%Dp_dry_nm = (m_s*6/pi()*1/1e6)^(1/3)*1e9;




