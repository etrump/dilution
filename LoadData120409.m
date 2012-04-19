function [Dp_bin_vect, N_tot_bin] = LoadData120409(TotalPop) %alpha-P 120326 exper
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

B1 = [562.759
600.121
48.9244
43.3069
78.0673
0
63.9433
58.0967
0
48.2436
66.9853
113.645
122.634
88.5217
180.269
402.308
398.864
553.61
811.825
1481.85
1789.81
2793.77
3705.83
4701.28
6758.33
9203.84
12137.5
15978
20399.6
27134.5
33992.8
43354.8
54087.9
66752.7
82018.4
100145
120315
145523
170849
201141
231679
264870
301508
338640
373092
408941
436676
463386
488269
509548
520646
529223
529039
521209
501172
465797
419714
367445
308587
249371
195712
147543
107959
74606.3
49603.9
31837.3
19358.6
10976.3
6018.35
3178.42
1428.83
600.631
259.286
129.551
65.4817
50.4956
43.2756
26.6624
16.515
26.7543
15.6785
35.7182
18.5801
1.51E+01
1.14E+01
1.87E+01
8.29E+00
1.99E+00
1.52E+01
1.10E+01
1.23E+01
4.18E+00
17.0189
8.6058
8.75739
8.91476
18.1565
17.8632
15.3996
29.4064
8.0322
4.93165
0
10.2045
10.4222
10.5715
5.36263
16.4504];
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

Dp_0m = 30; %[=]nm, start at 0.010 um
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
xn = Dp_fm;
x0 = Dp_0m;
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




