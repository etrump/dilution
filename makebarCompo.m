function makebarCompo(Y)

Total = sum(Y);
y1 = Y(1)/Total;
y2 = Y(2)/Total;

if length(Y)>2
    %xvector = [0.01 0.1 1 10];
    y3 = Y(3)/Total;
    y4 = Y(4)/Total;
    ymatrix1 = [y1 y2 y3 y4; [1 0 0 0]];
else
    %xvector = [0.01 1];
    ymatrix1 = [y1 y2;[1 0]];
end

barcolor = [0 0.4 0; 0 0.6 0; 0 0.8 0; 0 1 0];
figure1 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'InvertHardcopy','off','Colormap',barcolor,...
    'Color',[1 1 1]);

axes1 = axes('Parent',figure1,'YTick',[0 0.5 1],'XTick',zeros(1,0),...
    'LineWidth',3,'FontSize',16);
box('on');
hold('all');
bar1 = bar(ymatrix1,'BarLayout','stacked');
%set(bar1(1),'DisplayName','Composition, COC');


% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.6754 0.692 0.2616 0.1372]);

axis([0.5 1.5 0 1])
ylabel('Composition')
%xlabel('C* (ug/m3)')
 
