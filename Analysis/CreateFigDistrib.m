
function CreateFigDistrib(Time_hr)

global Tagg NumConcT DpT

Indy = length(find(Tagg<Time_hr*3600));

y = [NumConcT(Indy,:)];

x = DpT(Indy,:);

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'XScale','log','LineWidth',3,'FontSize',16);
box(axes1,'on');
hold(axes1,'all');

xlabel({'Dp (nm)'},'FontSize',16);
ylabel({'NumConc (#/m^3)'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)
axis([1 1e4 0 5e8]);
%axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','off');
%box(axes1,'on');
%hold(axes1,'all');

bar(x,y);

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/T_5hr.pdf'],'pdf'); 
