figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg DpT modelAtm DpT2 CpT VFR3

Tagg_hr = Tagg/3600;
DpT_nm = DpT;
PopIndy = ceil(modelAtm.Pop/2)-1;

%for i = 1:length(CpT(:,1))
    
%  DpT2_nm(i,1) = (6/pi*1/modelAtm.SOA.rho*(CpT(i,1)+CpT(i,2)+CpT(i,3)+CpT(i,4))./modelAtm.Pop1.NumConc0).^(1/3)*1e9;
%end

%plot1 = plot(Tagg_hr,DpT_nm(:,14),'LineWidth',3)
plot1 = plot(Tagg_hr,DpT_nm(:,PopIndy-1),Tagg_hr,DpT_nm(:,PopIndy),Tagg_hr,DpT_nm(:,PopIndy+1),'LineWidth',3)

set(plot1(1),'DisplayName','Background','Color',[0 0.6 0]);


% Create xlabel
xlabel({'Time from Dilution (hr)'},'FontSize',16);

% Create ylabel
ylabel({'D_{p,ave} (nm)'},'FontSize',16);

% Create legend
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

%legend1 = legend(axes1,'show');

%axis([0 3.5  0 370]);
%set(legend1,'Position',[0.5034 0.1988 0.37 0.1601]);
%set(legend1,'Location','NorthWest')

% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/Diam.pdf'],'pdf'); 


%V_0 = pi()/6*DpT_nm(1,14)^3;
totalVol = zeros(length(Tagg));

for i = 1:length(Tagg)
for j = 1:modelAtm.Pop
   totalVol(i) = totalVol(i) + TotalSusp(j)*pi()/6*DpT_nm(i,j)^3; 
end
    totalVFR(i) = totalVol(i)/totalVol(1);
end

V_0 = pi()/6*DpT_nm(1,PopIndy)^3;
V_02 = pi()/6*DpT_nm(1,PopIndy-1)^3;
V_03 = pi()/6*DpT_nm(1,PopIndy+1)^3;

for i = 1:length(CpT(:,1))
    
  %VFR4(i,1) = pi()/6*DpT_nm(i,14)^3*1/V_0;
  VFR4(i,1) = pi()/6*DpT_nm(i,PopIndy)^3*1/V_0;
  VFR4(i,2) = pi()/6*DpT_nm(i,PopIndy-1)^3*1/V_02;
  VFR4(i,3) = pi()/6*DpT_nm(i,PopIndy+1)^3*1/V_03;
end

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure2,'LineWidth',3,'FontSize',16);
box('on');
hold('all');
axis([0 7 0 1]);

%plot2 = plot(Tagg_hr,VFR4(:,2),Tagg_hr,VFR4(:,1),Tagg_hr,VFR4(:,3),'LineWidth',3)
plot2 = plot(Tagg_hr,totalVFR,'LineWidth',3)

set(plot2(1),'DisplayName','Background','Color',[0 0.6 0]);

VFR_0p1_3p0 = VFR4;
Tagg_hr_0p1_3p0 = Tagg_hr;

% Create xlabel
xlabel({'Time from Dilution (hr)'},'FontSize',16);

% Create ylabel
ylabel({'VFR'},'FontSize',16);

% Create legend
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/VFR.pdf'],'pdf'); 

VFRdata120409
