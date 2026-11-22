%% Bifurcation diagram in the (r, w_4)-plane (run after the continuation script)
clear
load('Equlibria_with__branching_off.mat')
physonly = true; 
%% Data: (r, w_4) with w_4 = x(end)
xpar_y2_0 = arrayfun(@(x)x.parameter(in.r),eq_y2_0.point);
ypar_y2_0 = arrayfun(@(p)p(end).x(end),eq_y2_0.point);

xpar_1234 = arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234 = arrayfun(@(p)p(end).x(end),eq_1234ref.point);

xpar_2340 = arrayfun(@(x)x.parameter(in.r),eq_2340ref.point);
ypar_2340 = arrayfun(@(p)p(end).x(end),eq_2340ref.point);

xpar_12_340 = arrayfun(@(x)x.parameter(in.r),eq_12_340ref.point);
ypar_12_340 = arrayfun(@(p)p(end).x(end),eq_12_340ref.point);

xpar_13_240 = arrayfun(@(x)x.parameter(in.r),eq_13_240ref.point);
ypar_13_240 = arrayfun(@(p)p(end).x(end),eq_13_240ref.point);

xpar_13_40 = arrayfun(@(x)x.parameter(in.r),eq_13_40br2ref.point);
ypar_13_40 = arrayfun(@(p)p(end).x(end),eq_13_40br2ref.point);

xpar_13_24_br2 = arrayfun(@(x)x.parameter(in.r),eq_13_24_br2ref.point);
ypar_13_24_br2 = arrayfun(@(p)p(end).x(end),eq_13_24_br2ref.point);

%% Remove points with y_i<0 (NaN breaks the lines and hides the markers)
isphys = @(br) arrayfun(@(p)all(p.x(1:4)>-1e-10),br.point) | ~physonly;
ypar_1234(~isphys(eq_1234ref))           = NaN;
ypar_2340(~isphys(eq_2340ref))           = NaN;
ypar_12_340(~isphys(eq_12_340ref))       = NaN;
ypar_13_240(~isphys(eq_13_240ref))       = NaN;
ypar_13_40(~isphys(eq_13_40br2ref))         = NaN;
ypar_13_24_br2(~isphys(eq_13_24_br2ref)) = NaN;

%% Settings
% colours for #unst = 0,1,2,3,>=4
cmap = [0.11 0.27 0.55;    % 0  dark blue  (stable)
        0.25 0.62 0.80;    % 1  light blue
        0.97 0.77 0.22;    % 2  yellow
        0.94 0.50 0.19;    % 3  orange
        0.84 0.19 0.15];   % >=4 red
ms     = 36;     % scatter marker area (pt^2)
stride = 3;      % plot every 3rd point as a marker (increase if crowded)
lw     = 0.6;    % width of the thin black branch lines

%% Plot
figure(1500); clf;
set(gcf,'Units','centimeters','Position',[2 2 17 10],'Color','w', ...
    'PaperUnits','centimeters','PaperSize',[17 10],'PaperPosition',[0 0 17 10]);
ax1 = axes; hold(ax1,'on');

% thin black lines through every branch
plot(xpar_y2_0,      ypar_y2_0,      'k-','LineWidth',lw);
plot(xpar_1234,      ypar_1234,      'k-','LineWidth',lw);
plot(xpar_2340,      ypar_2340,      'k-','LineWidth',lw);
plot(xpar_12_340,    ypar_12_340,    'k-','LineWidth',lw);
plot(xpar_13_240,    ypar_13_240,    'k-','LineWidth',lw);
plot(xpar_13_40,     ypar_13_40,     'k-','LineWidth',lw);
plot(xpar_13_24_br2, ypar_13_24_br2, 'k-','LineWidth',lw);

% markers coloured by #unst
i = 1:stride:numel(xpar_y2_0);
scatter(xpar_y2_0(i), ypar_y2_0(i), ms, cmap(min(nunst_st(i),4)+1,:), ...
    'o','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_1234);
scatter(xpar_1234(i), ypar_1234(i), ms, cmap(min(eq_1234tests(i),4)+1,:), ...
    '>','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_2340);
scatter(xpar_2340(i), ypar_2340(i), ms, cmap(min(eq_2340tests(i),4)+1,:), ...
    'd','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_12_340);
scatter(xpar_12_340(i), ypar_12_340(i), ms, cmap(min(eq_12_340tests(i),4)+1,:), ...
    'v','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_13_240);
scatter(xpar_13_240(i), ypar_13_240(i), ms, cmap(min(eq_13_240tests(i),4)+1,:), ...
    's','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_13_40);
scatter(xpar_13_40(i), ypar_13_40(i), ms, cmap(min(eq_13_40tests(i),4)+1,:), ...
    '^','filled','MarkerEdgeColor','k','LineWidth',0.4);

i = 1:stride:numel(xpar_13_24_br2);
scatter(xpar_13_24_br2(i), ypar_13_24_br2(i), 1.6*ms, cmap(min(eq_13_24_br2tests(i),4)+1,:), ...
    'p','filled','MarkerEdgeColor','k','LineWidth',0.4);

%% Axes
set(ax1,'FontSize',10,'TickLabelInterpreter','latex','Box','on', ...
    'TickDir','in','LineWidth',0.75,'Layer','top', ...
    'XMinorTick','on','YMinorTick','on','XLim',[0 3.5]);
yl = ylim(ax1);  ylim(ax1,[-0.04*yl(2) yl(2)]);
xlabel(ax1,'$r$','Interpreter','latex','FontSize',12);
ylabel(ax1,'$w_4$','Interpreter','latex','FontSize',12);

%% Legend 1: isotropy (marker shapes)
h1 = plot(ax1,NaN,NaN,'ko','MarkerFaceColor','w','MarkerSize',6);
h2 = plot(ax1,NaN,NaN,'kd','MarkerFaceColor','w','MarkerSize',6);
h3 = plot(ax1,NaN,NaN,'kv','MarkerFaceColor','w','MarkerSize',6);
h4 = plot(ax1,NaN,NaN,'ks','MarkerFaceColor','w','MarkerSize',6);
h5 = plot(ax1,NaN,NaN,'k^','MarkerFaceColor','w','MarkerSize',6);
h6 = plot(ax1,NaN,NaN,'kp','MarkerFaceColor','w','MarkerSize',8);
h7 = plot(ax1,NaN,NaN,'k>','MarkerFaceColor','w','MarkerSize',6);
lg1 = legend(ax1,[h1 h2 h3 h4 h5 h6 h7], ...
    {'$(1{=}2{=}3{=}4{=}0)$','$(2{=}3{=}4{=}0)$','$(12)(3{=}4{=}0)$', ...
     '$(13)(2{=}4{=}0)$','$(13)(4{=}0)$','$(13)(24)$','$(1234)$'}, ...
    'Interpreter','latex','FontSize',9,'Location','northwest');
title(lg1,'isotropy','Interpreter','latex');

%% Legend 2: stability (colours) -- needs a second, invisible axes
ax2 = axes('Position',get(ax1,'Position'),'Visible','off'); hold(ax2,'on');
k0 = plot(ax2,NaN,NaN,'s','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k','MarkerSize',7);
k1 = plot(ax2,NaN,NaN,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',7);
k2 = plot(ax2,NaN,NaN,'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerSize',7);
k3 = plot(ax2,NaN,NaN,'s','MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k','MarkerSize',7);
k4 = plot(ax2,NaN,NaN,'s','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerSize',7);
lg2 = legend(ax2,[k0 k1 k2 k3 k4], ...
    {'$\#\mathrm{unst}=0$','$\#\mathrm{unst}=1$','$\#\mathrm{unst}=2$', ...
     '$\#\mathrm{unst}=3$','$\#\mathrm{unst}\geq 4$'}, ...
    'Interpreter','latex','FontSize',9);
title(lg2,'stability','Interpreter','latex');

% place the stability legend just right of the isotropy legend
drawnow;
p1 = lg1.Position;  p2 = lg2.Position;
lg2.Position = [p1(1)+p1(3)+0.01, p1(2)+p1(4)-p2(4), p2(3), p2(4)];

%% Export
 exportgraphics(gcf,'multistrain_r_w4.pdf','ContentType','vector');
% exportgraphics(gcf,'multistrain_r_w4.png','Resolution',600);