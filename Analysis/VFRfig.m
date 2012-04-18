global modelAtm

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

mode_dex = ceil(modelAtm.Pop/2)+1;
mode_dex = length(DpT(1,:))
mode_dex = 12;
mode_dex = 15;
%mode_dex = 1;

Dp_init = DpT(1,mode_dex)

DpT(length(Tagg),mode_dex)

VFR = DpT(:,mode_dex).^3/Dp_init^3;

VFR_inf = VFR(length(Tagg));

%figure(25)
plot1 = plot(Tagg/3600,VFR)

set(plot1(1),'LineWidth',3,'DisplayName','Exper: 120307 \alphaP')
%set(plot1(1),'LineWidth',3,'DisplayName','Exper: 120202 \alphaP')
%set(plot1(1),'LineWidth',3,'DisplayName','Exper: 120326 \alphaP')

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'VFR'},'FontSize',16);

paperSize = [10 8];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast');


axis([0 7 0 1]);
% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/VFR_mode120307.pdf'],'pdf'); 
%saveas(gcf,['./figs/VFR_mode120202.pdf'],'pdf'); 
%saveas(gcf,['./figs/VFR_mode120326.pdf'],'pdf');
hold('off');


figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

%figure(26)
plot2 = plot(Tagg/3600, log(VFR)-log(VFR_inf))

set(plot2(1),'LineWidth',3,'DisplayName','Exper: 120307 \alphaP')
%set(plot2(1),'LineWidth',3,'DisplayName','Exper: 120202 \alphaP')
%set(plot2(1),'LineWidth',3,'DisplayName','Exper: 120326 \alphaP')

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'ln(VFR)-ln(VFR_{final})'},'FontSize',16);

paperSize = [10 8];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast');


%axis([0 7 0 1]);
% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/VFR_timescale120307.pdf'],'pdf');
%saveas(gcf,['./figs/VFR_timescale120202.pdf'],'pdf');
%saveas(gcf,['./figs/VFR_timescale120326.pdf'],'pdf');

figure(77)
plot(Tagg,DpT(:,mode_dex))