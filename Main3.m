function [T, Y] = Main3(tspan, y0)

%*System Parameters*******************************


global modelAtm

entries = length(y0);

fun=@Model3;

options=odeset('RelTol',1e-20,'AbsTol',1e-20,'NonNegative',[1:16]);
options=odeset('RelTol',1e-13,'AbsTol',1e-13,'NonNegative',[1:166]);
options=odeset('RelTol',1e-14,'AbsTol',1e-14,'NonNegative',[1:86]);
options=odeset('RelTol',1e-13,'AbsTol',1e-13,'NonNegative',[1:61]); %1 C-star bin (squalane), 20 size bins
options=odeset('RelTol',1e-13,'AbsTol',1e-13,'NonNegative',[1:entries]); %1 C-star bin (squalane), 100 size bins (401)  ... 200 bins (801)
%options=odeset('RelTol',1e-12,'AbsTol',1e-14);
%options=odeset('RelTol',1e-11,'AbsTol',1e-13);
%options=odeset('RelTol',1e-10,'AbsTol',1e-12,'NonNegative',[1:54]);
%options=odeset('RelTol',1e-8,'AbsTol',1e-11,'NonNegative',[1:187]);
%options=odeset('RelTol',1e-7,'AbsTol',1e-10);
%options=odeset('RelTol',1e-6,'AbsTol',1e-9);
%options=odeset('RelTol',1e-5,'AbsTol',1e-8);
%options=odeset('RelTol',1e-1,'AbsTol',1e-1);
%options=odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',[1:166]);
[T Y]=ode15s(fun,tspan,y0,options);





