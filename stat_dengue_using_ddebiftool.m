clear;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(fullfile('ddebiftool_snapshot'));
format compact
format short g
%%
parnames = {'r',  'e',  'beta',   'mu',   'sigma',   'gamma',  'tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.r,      in.e,     in.beta,        in.mu,      in.sigma,     in.gamma,  in.tau_delay])=... 
    [0.0,         1,        40,           0.02,            10,         0.7         0  ];
[nstrains,nvar]=deal(4,3);
ndim=nstrains*nvar;
x0=zeros(ndim,1); % x is nstrains x nvar
    parbd={'min_bound',[in.gamma, 0;  in.r,0.0],...
    'max_bound',   [in.gamma, 1;     in.r, 3.5],...
    'max_step',    [in.gamma, 0.05;   in.r, 0.05; 0,0.05]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=0.004;%par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% Tracking a branch of equilirium points
% We impose that all components are 0
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
[fun_y2_0,eq_y2_0,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.05,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',alleq0,'outputfuncs',true);
figure(1)
clf; hold on
eq_y2_0=br_contn(fun_y2_0,eq_y2_0,1000);
eq_y2_0=br_rvers(eq_y2_0);
eq_y2_0=br_contn(fun_y2_0,eq_y2_0,300);
[eq_y2_0,nunst_st]=br_stabl(fun_y2_0,eq_y2_0,0,0);
%% plot stability
xpar_y2_0=arrayfun(@(x)x.parameter(in.r),eq_y2_0.point);
ypar_y2_0=arrayfun(@(x)x.x(1),eq_y2_0.point);
biflocation=find(diff(nunst_st));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%% Extract stability information and plot results
figure(1);
clf
hold on 
plot(xpar_y2_0(nunst_st==0), ypar_y2_0(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(nunst_st==1), ypar_y2_0(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_y2_0(nunst_st>=2), ypar_y2_0(nunst_st>=2), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(biflocation), ypar_y2_0(biflocation), 'k.','Linewidth',3,'MarkerSize',30);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
eq_y2_0=br_remove_extracolumns(eq_y2_0);
%% Branching off, Full symmetry (1234)
ex_1234=permfix('stst','p1234',ndim,'x','perm',[1,2,3,4]);
ev_1234=permfix('stst','pv1234',ndim,'v','perm',[1,2,3,4]);
sbargs={'SetupFold.usercond',[alleq0,ev_1234],'SetupFold.initcond',[alleq0,ev_1234]};
[fun_eq1234,eq_1234,suc_1234]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',ex_1234,'outputfuncs',true,'print_residual_info',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
figure(1)
hold on
eq_1234=br_contn(fun_eq1234,eq_1234,1);
eq_1234=br_rvers(eq_1234);
eq_1234=br_contn(fun_eq1234,eq_1234,100);
[eq_1234,nunst_1234]=br_stabl(fun_eq1234,eq_1234);
[eq_1234ref,eq_1234tests,bifs_1234,bifind_1234]=MonitorChange(fun_eq1234,eq_1234,'min_iterations',5);
%% Branch off with  y2=y3=y4=0
cndv_2340=[permfix('stst','pv2340',nstrains,'v','perm',[2,3,4])];
cnd_2340=[permfix('stst','p2340',nstrains,'x','perm',[2,3,4])];
ye234q0=dde_sys_cond_create('name','y24_eq0','reference',false,'fun',@(x)x,'args',{'x',2:4});
sbargs={'SetupFold.usercond',[alleq0,cnd_2340,cndv_2340],'SetupFold.initcond',[alleq0,cndv_2340,cndv_2340]};
[fun_eq2340,eq_2340,suc_2340]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',[cnd_2340],'outputfuncs',true,parbd{:},'print_residual_info',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (12)(34) symmetry\n');
figure(1)
hold on
eq_2340=br_contn(fun_eq2340,eq_2340,100);
eq_2340=br_rvers(eq_2340);
eq_2340=br_contn(fun_eq2340,eq_2340,1);
[eq_2340ref,eq_2340tests,bifs_2340,bifind_2340]=...
    MonitorChange(fun_eq2340,eq_2340,'min_iterations',5);
%%
%% Branching off with y1=y3 and y2=y4=0
y2v0=dde_sys_cond_create('name','y2_v0','reference',false,'fun',@(x)x,'args',{'v',2});
y4eq0=dde_sys_cond_create('name','y4_eq0','reference',false,'fun',@(x)x,'args',{'x',4});
y2eq0=dde_sys_cond_create('name','y4_v0','reference',false,'fun',@(x)x,'args',{'x',2});
cndv_13=permfix('stst','pv13',nstrains,'v','perm',[1,3]);
cnd_13=permfix('stst','p13',nstrains,'x','perm',[1,3]);
sbargs={'SetupFold.usercond',[alleq0,cndv_13,y2v0],'SetupFold.initcond',[alleq0,cndv_13]};
[fun_eq13_240,eq_13_240,suc_13]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',[cnd_13,cndv_13,y4eq0,y2eq0],'initcond',[cnd_13,cndv_13],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)},...
    'print_residual_info',true,parbd{:});
%%
eq_13_240.point(2).x(1:4,:)
figure(1)
hold on
eq_13_240=br_contn(fun_eq13_240,eq_13_240,1);
eq_13_240=br_rvers(eq_13_240);
eq_13_240=br_contn(fun_eq13_240,eq_13_240,100);
[eq_13_240,nunst_13_240]=br_stabl(fun_eq13_240,eq_13_240);
[eq_13_240ref,eq_13_240tests,bifs_13_240,bifind_13_240]=...
    MonitorChange(fun_eq13_240,eq_13_240,'min_iterations',5);
%%  Second Branching Off (13)(24)
cndv_13_24=[permfix('stst','pv13',nstrains,'v','perm',[1,3]),...
           permfix('stst','pv24',nstrains,'v','perm',[2,4])];
cnd_13_24=[permfix('stst','p13',nstrains,'x','perm',[1,3]),...
           permfix('stst','p24',nstrains,'x','perm',[2,4])];
sbargs_13_br2={'SetupFold.usercond',[cndv_13_24],'SetupFold.initcond',[cndv_13_24]};
[fun_eq13_br2,eq_13_24_br2,suc_13_24]=SetupStstFrom_fold(funcs,eq_13_240ref,bifind_13_240(3),...
    'usercond',[cnd_13_24,cndv_13_24],'initcond',[cnd_13_24,cndv_13_24],'outputfuncs',true,...
    'step',0.1,sbargs_13_br2{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)},...
    'print_residual_info',true,parbd{:});
fprintf('Branch off from symmetry breaking, ~ (13)(24) symmetry\n');
eq_13_24_br2.point(2).x(1:4,:)
%%
figure(1)
hold on
eq_13_24_br2=br_contn(fun_eq13_br2,eq_13_24_br2,30);
eq_13_24_br2=br_rvers(eq_13_24_br2);
eq_13_24_br2=br_contn(fun_eq13_br2,eq_13_24_br2,30);
[eq_13_24_br2ref,eq_13_24_br2tests,bifs_13_24,bifind_13_24]=...
    MonitorChange(fun_eq13_br2,eq_13_24_br2,'min_iterations',5);
%%

eq_1234ref=br_remove_extracolumns(eq_1234ref);
eq_1234=br_remove_extracolumns(eq_1234);
eq_2340ref=br_remove_extracolumns(eq_2340ref);
eq_2340=br_remove_extracolumns(eq_2340);
eq_13_240ref=br_remove_extracolumns(eq_13_240ref);
eq_13_240=br_remove_extracolumns(eq_13_240);
eq_13_24_br2ref=br_remove_extracolumns(eq_13_24_br2ref)

%%
xpar_y2_0=arrayfun(@(x)x.parameter(in.r),eq_y2_0.point);
ypar_y2_0=arrayfun(@(p)p(end).x(end),eq_y2_0.point);
biflocation=find(diff(nunst_st));
%%%%%%%%%%%%
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234=arrayfun(@(p)p(end).x(end),eq_1234ref.point);
% xpar_12_34=arrayfun(@(x)x.parameter(in.r),eq_12_34ref.point);
% ypar_12_34=arrayfun(@(p)sum(p.x(end-3:end)),eq_12_34ref.point);
% xpar_12_34=arrayfun(@(x)x.parameter(in.r),eq_12_34ref.point);
% ypar_12_34=arrayfun(@(p)sum(p.x(end-3:end)),eq_12_34ref.point);
%%%%%%%%%%%
xpar_2340=arrayfun(@(x)x.parameter(in.r),eq_2340ref.point);
ypar_2340=arrayfun(@(p)p(end).x(end),eq_2340ref.point);
%%%%%%%%%%
xpar_13_240=arrayfun(@(x)x.parameter(in.r),eq_13_240ref.point);
ypar_13_240=arrayfun(@(p)p(end).x(end),eq_13_240ref.point);
%%%%%%%%%%%

xpar_13_24_br2=arrayfun(@(x)x.parameter(in.r),eq_13_24_br2ref.point);
ypar_13_24_br2=arrayfun(@(p)p(end).x(end),eq_13_24_br2ref.point);
%%
save('Equlibria_with__branching_off.mat')
%%
figure(1500)
clf;
hold on
plot(xpar_y2_0(nunst_st==0), ypar_y2_0(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(nunst_st==1), ypar_y2_0(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_y2_0(nunst_st>=2), ypar_y2_0(nunst_st>=2), 'k:','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(biflocation), ypar_y2_0(biflocation), 'k.','Linewidth',3,'MarkerSize',30);
%%%%%%%%%%%
plot(xpar_1234(eq_1234tests>=1),ypar_1234(eq_1234tests>=1),'-','LineWidth',2)
plot(xpar_2340(eq_2340tests>=1),ypar_2340(eq_2340tests>=1),'LineWidth',2)
plot(xpar_13_240(eq_13_240tests==0),ypar_13_240(eq_13_240tests==0),'-','LineWidth',1)
plot(xpar_13_240(eq_13_240tests>=1),ypar_13_240(eq_13_240tests>=1),'.','LineWidth',1)

plot(xpar_13_24_br2,ypar_13_24_br2,'.','LineWidth',2) 



%%
% y4eq0=dde_sys_cond_create('name','y4_eq0','reference',false,'fun',@(x)x,'args',{'x',4});
% y4veq0=dde_sys_cond_create('name','y4_v0','reference',false,'fun',@(x)x,'args',{'v',4});
% %alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
% cndv_13=[permfix('stst','pv13',nstrains,'v','perm',[1,3])];
% cnd_13=[permfix('stst','p13',nstrains,'x','perm',[1,3])];
% sbargs_13_br2={'SetupFold.usercond',[y4eq0,cnd_13,cndv_13],'SetupFold.initcond',[y4veq0,y4eq0,cndv_13]};
% [fun_eq13_br2,eq_13_br2,suc_13br2]=SetupStstFrom_fold(funcs,eq_13_240ref,bifind_13_240(3),...
%     'usercond',[cnd_13,y4eq0],'initcond',[cnd_13],'outputfuncs',true,...
%     'step',0.1,sbargs_13_br2{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)},...
%     'print_residual_info',true,parbd{:});
% fprintf('Branch off from symmetry breaking, ~ (13)(24) symmetry\n');
% eq_13_br2.point(2).x(1:4,:)
% %%
% figure(1)
% hold on
% eq_13_br2=br_contn(fun_eq13_br2,eq_13_br2,20);

%% % %%%%%%%%%%%%%%%%%%%%
%% plot_multistrain_bifdiag.m
% Publication-quality bifurcation diagram for the D4-symmetric multi-strain
% model. Run AFTER the DDE-Biftool continuation script, so that the
% workspace contains:
%   eq_y2_0, nunst_st, biflocation                      (trivial branch)
%   eq_1234ref,    eq_1234tests,    bifind_1234         (1234)
%   eq_2340ref,    eq_2340tests,    bifind_2340         (234)
%   eq_13_240ref,  eq_13_240tests,  bifind_13_240       (13), y2=y4=0
%   eq_13_24_br2ref, eq_13_24_br2tests, bifind_13_24    (13)(24)
%   funcs, in
%
% Stability convention: tests == 0 -> stable (solid),
%                       tests  > 0 -> unstable (dashed).
% Label convention: (13)(24) means y1 = y3 and y2 = y4, etc.

%% ------------------------------------------------------------------------
%  User settings
%  ------------------------------------------------------------------------
meas    = @(pt) pt.x(1);      % plotted quantity (y_1); change if needed
ylab    = '$y_1$';
xlab    = '$r$';

figW    = 17.0;               % cm  (double column ~17 cm, single ~8.5 cm)
figH    = 9.0;                % cm
fs      = 10;                 % font size (pt) for ticks/labels
fsLeg   = 9;                  % legend font size
lwS     = 1.8;                % line width, stable
lwU     = 1.3;                % line width, unstable
msBif   = 6;                  % marker size, bifurcation points
xl      = [0 3.5];
yl      = [-0.03 0.9];
outname = 'multistrain_bifdiag';   % exported as .pdf (vector) and .png

% Okabe-Ito colour-blind-safe palette (also readable in greyscale)
C.triv  = [0   0   0  ];          % black
C.p1234 = [0   114 178]/255;      % blue
C.p234  = [0   158 115]/255;      % bluish green
C.p13   = [213 94  0  ]/255;      % vermilion
C.p1324 = [204 121 167]/255;      % reddish purple

%% ------------------------------------------------------------------------
%  Collect branches: {branch, stability, bif. indices, colour, legend label}
%  ------------------------------------------------------------------------
ir = in.r;
B = {};
B{end+1} = mkbr(funcs, eq_y2_0,         nunst_st,          biflocation,   C.triv,  ...
    'trivial, $y_j=0$',              ir, meas);
B{end+1} = mkbr(funcs, eq_1234ref,      eq_1234tests,      bifind_1234,   C.p1234, ...
    '$(1234)$',                      ir, meas);
B{end+1} = mkbr(funcs, eq_2340ref,      eq_2340tests,      bifind_2340,   C.p234,  ...
    '$(234)$',                       ir, meas);
B{end+1} = mkbr(funcs, eq_13_240ref,    eq_13_240tests,    bifind_13_240, C.p13,   ...
    '$(13)$, $y_2=y_4=0$',           ir, meas);
B{end+1} = mkbr(funcs, eq_13_24_br2ref, eq_13_24_br2tests, bifind_13_24,  C.p1324, ...
    '$(13)(24)$',                    ir, meas);

%% ------------------------------------------------------------------------
%  Figure and axes
%  ------------------------------------------------------------------------
fig = figure(1500); clf(fig);
set(fig, 'Units','centimeters', 'Position',[2 2 figW figH], 'Color','w', ...
    'PaperUnits','centimeters', 'PaperSize',[figW figH], ...
    'PaperPosition',[0 0 figW figH], 'Renderer','painters');
ax = axes(fig); hold(ax,'on');

% faint reference line y = 0 (drawn first, behind everything)
plot(ax, xl, [0 0], '-', 'Color',[0.85 0.85 0.85], 'LineWidth',0.5, ...
    'HandleVisibility','off');

% coloured branches first, trivial branch last so it sits on top
order = [2:numel(B) 1];
hBr = gobjects(numel(B),1);
for k = order
    hBr(k) = plotbranch(ax, B{k}, lwS, lwU);
end

% bifurcation points on top of all lines
for k = order
    b = B{k};
    idx = b.bif(b.bif>=1 & b.bif<=numel(b.r));
    if isempty(idx), continue; end
    plot(ax, b.r(idx), b.y(idx), 'o', 'MarkerSize',msBif, ...
        'MarkerFaceColor',b.col, 'MarkerEdgeColor','k', 'LineWidth',0.6, ...
        'HandleVisibility','off');
end

% primary (symmetry-breaking) bifurcation from the trivial state
bt = B{1};
if ~isempty(bt.bif)
    plot(ax, bt.r(bt.bif(1)), bt.y(bt.bif(1)), 'o', 'MarkerSize',msBif+2, ...
        'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'LineWidth',0.8, ...
        'HandleVisibility','off');
end

%% ------------------------------------------------------------------------
%  Axes styling
%  ------------------------------------------------------------------------
set(ax, 'FontSize',fs, 'TickLabelInterpreter','latex', ...
    'Box','on', 'TickDir','in', 'TickLength',[0.012 0.012], ...
    'LineWidth',0.75, 'Layer','top', ...
    'XMinorTick','on', 'YMinorTick','on', ...
    'XLim',xl, 'YLim',yl, 'XColor','k', 'YColor','k');
xlabel(ax, xlab, 'Interpreter','latex', 'FontSize',fs+2);
ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',fs+2);

%% ------------------------------------------------------------------------
%  Legend: one colour per symmetry type + line-style key for stability
%  ------------------------------------------------------------------------
hS = plot(ax, NaN, NaN, 'k-',  'LineWidth',lwS);
hU = plot(ax, NaN, NaN, 'k--', 'LineWidth',lwU);
hB = plot(ax, NaN, NaN, 'o', 'MarkerSize',msBif, 'MarkerFaceColor',[.6 .6 .6], ...
    'MarkerEdgeColor','k', 'LineWidth',0.6);
labels = [cellfun(@(b) b.name, B, 'UniformOutput',false), ...
          {'stable','unstable','bifurcation point'}];
lg = legend(ax, [hBr(:); hS; hU; hB], labels, ...
    'Interpreter','latex', 'FontSize',fsLeg, 'Location','northwest', ...
    'NumColumns',2, 'Box','off');
lg.ItemTokenSize = [22 10];

%% ------------------------------------------------------------------------
%  Export (vector PDF for the journal, 600 dpi PNG for slides/preview)
%  ------------------------------------------------------------------------
drawnow;
if exist('exportgraphics','file')                 % R2020a+
    exportgraphics(fig, [outname '.pdf'], 'ContentType','vector', ...
        'BackgroundColor','white');
    exportgraphics(fig, [outname '.png'], 'Resolution',600, ...
        'BackgroundColor','white');
else
    print(fig, outname, '-dpdf', '-painters');
    print(fig, outname, '-dpng', '-r600');
end
fprintf('Saved %s.pdf and %s.png\n', outname, outname);

%% ========================================================================
%  Local functions
%  ========================================================================
function b = mkbr(funcs, br, tests, bifind, col, name, ir, meas)
% Extract (r, y, nunst) from a DDE-Biftool branch.
pts   = br.point;
b.r   = arrayfun(@(p) p.parameter(ir), pts);
b.y   = arrayfun(meas, pts);
b.col = col;
b.name = name;
b.bif = bifind(:)';
b.nunst = getnunst(funcs, br, tests);
end

function nunst = getnunst(funcs, br, tests)
% Number of unstable eigenvalues per point (0 = stable).
n = numel(br.point);
if (isnumeric(tests) || islogical(tests)) && numel(tests) == n
    nunst = double(tests(:)');
    return
end
warning('plot_multistrain_bifdiag:stability', ...
    ['Stability vector length (%d) does not match branch length (%d); ', ...
     'recomputing with br_stabl.'], numel(tests), n);
[~, nunst] = br_stabl(funcs, br, 0, 0);
nunst = double(nunst(:)');
end

function h = plotbranch(ax, b, lwS, lwU)
% Solid where stable, dashed where unstable. Masking with NaN (instead of
% logical indexing) prevents spurious lines joining separate segments;
% the dashed piece is extended by one point at every stability change so
% the two styles meet without a gap.
st = (b.nunst == 0);
ys = b.y;  ys(~st) = NaN;
yu = b.y;  yu(st)  = NaN;
for k = find(diff(st))
    yu(k)   = b.y(k);
    yu(k+1) = b.y(k+1);
end
h = plot(ax, b.r, ys, '-', 'Color',b.col, 'LineWidth',lwS);
plot(ax, b.r, yu, '--', 'Color',b.col, 'LineWidth',lwU, ...
    'HandleVisibility','off');
if all(~st)                      % fully unstable branch: legend shows colour
    h = plot(ax, NaN, NaN, '-', 'Color',b.col, 'LineWidth',lwS);
end
end