function formatgraph(plothan,xlab,ylab,axisbound,tname,fontsize,linewidth)
%This function formats a graph in the desired fashion.  The order of input
%variables are (plothandle,xlabel,ylabel,axisbound,titlename,fontsize,linewidth)
%titlename, fontsize, and linewidth are all optional.  This also makes the
%plot background white.

if nargin < 7; linewidth = 2; end
if nargin < 6; fontsize = 18; end
if nargin < 5; tname = '   '; end
if nargin > 3
    if isempty(axisbound) == 1;
    else
        axis(axisbound); 
    end
end
set(plothan,'LineWidth',linewidth);
set(gcf,'Color',[1 1 1]);
set(gca,'Fontsize',fontsize);
xlabel(xlab);
set(gca,'Fontsize',fontsize);
ylabel(ylab);
set(gca,'Fontsize',fontsize);
title(tname);
set(gca,'Fontsize',fontsize);