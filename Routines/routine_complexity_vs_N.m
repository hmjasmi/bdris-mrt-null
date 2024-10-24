clc
clear
close all
N = 8:8:64;
K = 2;
NG = 2;
G = N./NG;
t = 100;
t1 = 50;
t2 = 20;
t3 = 5;
bigO_Null_FC = t*N.^3;
bigO_Null_NG_2 = t*G*NG.^3;
bigO_Null_SC = t*N;

bigO_MRT_FC = N.^3;
bigO_MRT_NG_2 = G*NG.^3;
bigO_MRT_SC = K*N;

bigO_Joint_FC = t1*t2*N.^3;
bigO_Joint_NG_2 = t1*t2*t3*G*NG.^3;
bigO_Joint_SC = t1*t2*N;

bigO_MaxF_FC = N.^6;
bigO_MaxF_NG_2 = G.^3*NG.^6;
bigO_MaxF_SC = N.^3;
%%
figure
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 3.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 3.5]);
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(gcf,'color','w');
t = tiledlayout(1,1,'TileSpacing','tight');
t.Padding = 'compact';
ax1 = axes(t);
ax1.ColorOrderIndex = 1;
semilogy(ax1,NaN,NaN,'-','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,'-.','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,':','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,'--','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,'sk','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,'^k','LineWidth',1.25); hold on
semilogy(ax1,NaN,NaN,'ok','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
semilogy(ax1,N,bigO_Joint_FC,'-s','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MaxF_FC,'-.s','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_Null_FC,':s','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MRT_FC,'--s','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
semilogy(ax1,N,bigO_Joint_NG_2,'-^','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MaxF_NG_2,'-.^','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_Null_NG_2,':^','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MRT_NG_2,'--^','LineWidth',1.25); hold on
ax1.ColorOrderIndex = 1;
semilogy(ax1,N,bigO_Joint_SC,'-o','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MaxF_SC,'-.o','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_Null_SC,':o','LineWidth',1.25); hold on
semilogy(ax1,N,bigO_MRT_SC,'--o','LineWidth',1.25); hold on

ax1.XGrid = 'on';
ax1.YGrid = 'on';
xlabel(ax1,'Number of REs','interpreter','tex','fontsize',12)
ylabel(ax1,'Computational complexity','interpreter','tex','fontsize',12)
legend(ax1,'Joint','Max-F','Null','MRT','FC','GC (Group Size 2)','SC',...
           'interpreter','tex','fontsize',10,'location','best') 
% ax1.YScale = "log";
ax1.XLim = [min(N) max(N)];
ax1.YLim = [5e0 1e14];
ax1.FontSize = 10;
ax1.LineWidth = 0.75;
ax1.XTick = N;
ax1.YTick = [1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14];

print(gcf,'Figures_pub/Fig_complexity_vs_N.eps','-depsc','-r800')



