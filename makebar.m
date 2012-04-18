function makebar(Y, T, t)

xvector = [0.01 1];
ymatrix1 = [Y(t,5) Y(t,1)
    Y(t,7) Y(t,3)] %Caer Cvap
Top1 = Y(t,9)+Y(t,11);
y11 = Y(t,9)/Top1;
y12 = Y(t,11)/Top1;
Top2 = Y(t,5)+Y(t,7);
y21 = Y(t,5)/Top2;
y22 = Y(t,7)/Top2;
ymatrix2 = [y21 y22
    y11 y12]


%first make bar for C*=-2, background
barcolor = [0 1 0; 1 1 1];
figure1 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'InvertHardcopy','off','Colormap',barcolor,...
    'Color',[1 1 1]);

axes1 = axes('Parent',figure1,'XTickLabel',{'1e-2','1e+0'},...
    'XTick',[1 2],'LineWidth',3,'FontSize',16);
box('on');
hold('all');
bar1 = bar(ymatrix1,'BarLayout','stacked');
set(bar1(1),'DisplayName','C_B_G');
set(bar1(2),'DisplayName','C_v_a_p');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.6754 0.692 0.2616 0.1372]);


title(['time=',num2str(T(t)/3600),' hr'])
ylabel('Concentration (ug/m3)')
xlabel('C* (ug/m3)')

barcolor = [0 0.6 0; 0 1 0];
%ymatrix2 = randn(1,3);
figure2 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'InvertHardcopy','off','Colormap',barcolor,...
    'Color',[1 1 1]);
axes2 = axes('Parent',figure2,'XTickLabel',{'BG','Particle Coating'},...
    'XTick',[1 2],'LineWidth',3,'FontSize',16);
box('on');
hold('all');
bar2 = bar(ymatrix2,'BarLayout','stacked');
set(bar2(1),'DisplayName','C* = 1e-2');
set(bar2(2),'DisplayName','C* = 1e+0');

% Create legend
legend2 = legend(axes2,'show');
set(legend2,'Position',[0.6754 0.692 0.2616 0.1372]);
title(['time=',num2str(T(t)/3600),' hr'])
ylabel('Composition')





