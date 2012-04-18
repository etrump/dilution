
function Ans = Debye(x)

%Tabulated solutions to the Debye function: D1(x) =
%1/x*integral(t*dt/(e^t-1)) .. interation limits: 0,x

X = [0
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
1.1
1.2
1.3
1.4
1.6
1.8
2
2.2
2.4
2.6
2.8
3
3.2
3.4
3.6
3.8
4
4.2
4.4
4.6
4.8
5
5.5
6
6.5
7
7.5
8
8.5
9
9.5
10]';

FX = [1
0.975278
0.951111
0.927498
0.904437
0.881927
0.859964
0.838545
0.817665
0.79732
0.777505
0.758213
0.739438
0.721173
0.703412
0.669366
0.637235
0.606947
0.578427
0.551596
0.526375
0.502682
0.480435
0.459555
0.439962
0.42158
0.404332
0.388148
0.372858
0.358696
0.345301
0.332713
0.320876
0.29424
0.27126
0.251331
0.233948
0.218698
0.205239
0.193294
0.182633
0.173068
0.164443]';

IndexLower = length(find(X<x));
IndexUpper = IndexLower+1;


%Eventually I should interpolate... but for now I will just use the Lower
Index = IndexUpper;

if Index<1
    Index = 1;
    Ans = FX(Index);
elseif Index>length(FX)
    Index = length(FX);
    Ans = FX(Index);
else
    Low = FX(IndexLower);
    High = FX(IndexUpper);
    Ans = Low + (x-X(IndexLower))*(FX(IndexUpper)-FX(IndexLower))/(X(IndexUpper)-X(IndexLower));
end




%Ans = FX(Index);

