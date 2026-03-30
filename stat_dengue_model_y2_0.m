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
    [0.0,         1,        40,           0.02,            10,         0.54         0  ];
[nstrains,nvar]=deal(4,3);
ndim=nstrains*nvar;
x0=zeros(ndim,1); % x is nstrains x nvar
    parbd={'min_bound',[in.gamma, 0;  in.r,0.0],...
    'max_bound',   [in.gamma, 1;     in.r, 3.5],...
    'max_step',    [in.gamma, 0.05;   in.r, 0.05; 0,0.05]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
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
%%
% Extract stability information and plot results
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
%%
eq_y2_0=br_remove_extracolumns(eq_y2_0);
%save('computation_without_delay.mat');
%%
%% Invariant subspaces with respect to the group D_4:
%(1) Fully symetric invariant subspace <(1,1,1,1)> invariant under the
%action of the subgroup <r=(1234)>:={e,r,r^2,r^3}
% odd invariant subspace <(1,-1,1,-1)>  invariant under the
%action of the subgroup :={e,(13)(24),(13),(24)}, expected to detect Hopf
%bifurcation
check_sb=@(cnd)SetupFold(fun_y2_0,eq_y2_0,biflocation(1),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
%% Possible Subgroub of the group D_{4}, here we did not find an invariant subspace...
% with repect to any of the possible subgroups of D4
% {e,(13)(24)}
ev13_24=[permfix('stst','p13',ndim,'v','perm',[1,3]),...
         permfix('stst','p24',ndim,'v','perm',[2,4])];
sv_13_24=check_sb(ev13_24)'
% {e,(14)(23)}
ev14_23=[permfix('stst','p14',ndim,'v','perm',[1,4]),...
         permfix('stst','p23',ndim,'v','perm',[2,3])];
sv_14_23=check_sb(ev14_23)'
% {e,(12)(34)}
ev12_34=[permfix('stst','p12',ndim,'v','perm',[1,2]),...
         permfix('stst','p34',ndim,'v','perm',[3,4])];
sv_12_34=check_sb(ev12_34)'
% {e,(24)}
ev24=permfix('stst','p24',ndim,'v','perm',[2,4]);
sv_24=check_sb(ev24)'
% {e,(13)}
ev13=permfix('stst','p13',ndim,'v','perm',[1,3]);
[sv_13,se_13]=check_sb(ev13)
%
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
sv_1234=check_sb(ev_1234)'
% add additional constraints where y=(z=)0 for some component
y2veq0=dde_sys_cond_create('name','y2_eq0','reference',false,'fun',@(x)x,'args',{'v',2});
sv_13_24_y0=check_sb([ev13_24,y2veq0])'
sv_12_34_y0=check_sb([ev12_34,y2veq0])'
sv_13_y0=check_sb([ev13,y2veq0])'
%% Branch off in three different directions: 1st 1=2, 3=4, y2=0
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
cndv_12_34=[permfix('stst','p12',nstrains,'v','perm',[1,2]),...
           permfix('stst','p34',nstrains,'v','perm',[3,4])];
cnd_12_34=[permfix('stst','p12',nstrains,'x','perm',[1,2]),...
           permfix('stst','p34',nstrains,'x','perm',[3,4])];
y2eq0=dde_sys_cond_create('name','y2_eq0','reference',false,'fun',@(x)x,'args',{'x',2});
sbargs={'SetupFold.usercond',[alleq0,cndv_12_34,y2veq0],'SetupFold.initcond',[alleq0,cndv_12_34,y2veq0]};
[fun_eq12_34,eq_12_34,suc_12_34]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',[cnd_12_34,y2eq0],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (12)(34) symmetry\n');
%%
figure(1)
hold on
eq_12_34=br_contn(fun_eq12_34,eq_12_34,100);
eq_12_34=br_rvers(eq_12_34);
eq_12_34=br_contn(fun_eq12_34,eq_12_34,100);
eq_12_34=br_rvers(eq_12_34);
%[eq_12_34,nunst_12_34]=br_stabl(fun_eq12_34,eq_12_34);
[eq_12_34ref,eq_12_34tests,bifs_12_34,bifind_12_34]=...
    MonitorChange(fun_eq12_34,eq_12_34,'min_iterations',5);
%% Branch off in three different directions: 2nd 1=3, 2=4, y2=0
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
cndv_13_24=[permfix('stst','p13',nstrains,'v','perm',[1,3]),...
           permfix('stst','p24',nstrains,'v','perm',[2,4])];
cnd_13_24=[permfix('stst','p13',nstrains,'x','perm',[1,3]),...
           permfix('stst','p24',nstrains,'x','perm',[2,4])];
y2eq0=dde_sys_cond_create('name','y2_eq0','reference',false,'fun',@(x)x,'args',{'x',2});
sbargs={'SetupFold.usercond',[alleq0,cndv_13_24,y2veq0],'SetupFold.initcond',[alleq0,cndv_13_24,y2veq0]};
[fun_eq13_24,eq_13_24,suc_13_24]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',[cnd_13_24,y2eq0],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (13)(24) symmetry\n');
%%
figure(1)
hold on
eq_13_24=br_contn(fun_eq13_24,eq_13_24,100);
eq_13_24=br_rvers(eq_13_24);
eq_13_24=br_contn(fun_eq13_24,eq_13_24,100);
eq_13_24=br_rvers(eq_13_24);
%[eq_13_24,nunst_13_24]=br_stabl(fun_eq13_24,eq_13_24);
[eq_13_24ref,eq_13_24tests,bifs_13_24,bifind_13_24]=...
    MonitorChange(fun_eq13_24,eq_13_24,'min_iterations',5);
%% Branch off in three different directions: 3rd 1=2=3=4
ex_1234=permfix('stst','p1234',ndim,'x','perm',[1,2,3,4]);
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
sbargs={'SetupFold.usercond',[alleq0,ev_1234],'SetupFold.initcond',[alleq0,ev_1234]};
[fun_eq1234,eq_1234,suc_1234]=SetupStstFrom_fold(funcs,eq_y2_0,biflocation,...
    'usercond',ex_1234,'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (1234) symmetry\n');
%%
figure(1)
hold on
eq_1234=br_contn(fun_eq1234,eq_1234,100);
eq_1234=br_rvers(eq_1234);
eq_1234=br_contn(fun_eq1234,eq_1234,100);
eq_1234=br_rvers(eq_1234);
[eq_1234,nunst_1234]=br_stabl(fun_eq1234,eq_1234);
[eq_1234ref,eq_1234tests,bifs_1234,bifind_1234]=...
    MonitorChange(fun_eq1234,eq_1234,'min_iterations',5);

%%
figure(2);clf
plotopts={'y',@(p)sum(p.x(end-3:end))};
Plot2dBranch(eq_y2_0,plotopts{:});
hold on;
Plot2dBranch(eq_1234ref,plotopts{:},'DisplayName','eq_1234','color',[0,0,0.5])
Plot2dBranch(eq_12_34ref,plotopts{:},'color',[0.5,0,0])
Plot2dBranch(eq_13_24ref,plotopts{:},'DisplayName','eq_13_24ref','color',[0,0.5,0])
xlabel('x');ylabel('sum(w)');
%%
xpar_y2_0=arrayfun(@(x)x.parameter(in.r),eq_y2_0.point);
ypar_y2_0=arrayfun(@(p)sum(p.x(end-3:end)),eq_y2_0.point);
biflocation=find(diff(nunst_st));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%%
% Extract stability information and plot results
clrs=orderedcolors("gem12");
figure(1);
clf
hold on
xpar_12_34=arrayfun(@(x)x.parameter(in.r),eq_12_34ref.point);
ypar_12_34=arrayfun(@(p)sum(p.x(end-3:end)),eq_12_34ref.point);
%figure(33)

plt1_12_34=plot(xpar_12_34(eq_12_34tests==0), ypar_12_34(eq_12_34tests==0),'-','Color','red','Linewidth',3,'MarkerSize',10);
plt2_12_34=plot(xpar_12_34(eq_12_34tests>=1), ypar_12_34(eq_12_34tests>=1), '--','Color',0.5*[1 1 1],'Linewidth',3,'MarkerSize',10);
plot(xpar_12_34(bifind_12_34), ypar_12_34(bifind_12_34), 'k.','Linewidth',3,'MarkerSize',30);
% %%%%%%

zro1=find(eq_1234tests==0,1,'first');
one1=find(eq_1234tests==1,1,'last');
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234=arrayfun(@(p)sum(p.x(end-3:end)),eq_1234ref.point);
plt3_1234=plot(xpar_1234(1:zro1), ypar_1234(1:zro1), '--','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
plt1_1234=plot(xpar_1234(zro1:one1), ypar_1234(zro1:one1), ':','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
plt2_1234=plot(xpar_1234(one1:end), ypar_1234(one1:end), '--','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
plot(xpar_1234(bifind_1234(3)), ypar_1234(bifind_1234(3)), '.','Color',clrs(3,:),'Linewidth',3,'MarkerSize',30);
%%%%%%%%%%%%%%
z1=find(eq_13_24tests==0,1,'first');
z2=find(eq_13_24tests==0,1,'last');
xpar_13_24=arrayfun(@(x)x.parameter(in.r),eq_13_24ref.point);
ypar_13_24=arrayfun(@(p)sum(p.x(end-3:end)),eq_13_24ref.point);
plt1_13_24=plot(xpar_13_24(z1:z2), ypar_13_24(z1:z2), '-','Color',clrs(1,:),'Linewidth',3,'MarkerSize',10);
plt2_13_24=plot(xpar_13_24(1:z1), ypar_13_24(1:z1), '--','Color',clrs(1,:),'Linewidth',2,'MarkerSize',10);
plot(xpar_13_24(z2:end), ypar_13_24(z2:end), '--','Color',clrs(1,:),'Linewidth',2,'MarkerSize',10);
plot(xpar_13_24(bifind_13_24(3)), ypar_13_24(bifind_13_24(3)), '.','Color','r','Linewidth',3,'MarkerSize',30);
hold on 
plot(xpar_y2_0(nunst_st==0), ypar_y2_0(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(nunst_st==1), ypar_y2_0(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_y2_0(nunst_st>=2), ypar_y2_0(nunst_st>=2), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_y2_0(biflocation), ypar_y2_0(biflocation), 'k.','Linewidth',3,'MarkerSize',30);

set_bif={'LineWidth',2,'Box','on','FontSize',14};
ylim([-1,2.5])
set(gca, set_bif{:})
legend('k~(12)(34) unstable','','','k~(1234) unstable','','single-zero (fold)','k~(13)(24) stable','k~(13)(24) unstable',... 
          '','double-zero', 'stable','unstable','bif with 4 zero eiganvalues',...
          'FontSize',17,'FontWeight','normal','location','northwest')
xlabel('$r$','interpreter','latex','FontSize',17,FontWeight='bold')
ylabel('$y_{1}$','interpreter','latex',FontSize=17,FontWeight='bold')
%title('(a)');
%%
eq_1234ref=br_remove_extracolumns(eq_1234ref);
eq_1234=br_remove_extracolumns(eq_1234);
eq_12_34=br_remove_extracolumns(eq_12_34);
eq_12_34ref=br_remove_extracolumns(eq_12_34ref);
eq_13_24=br_remove_extracolumns(eq_13_24);
eq_13_24ref=br_remove_extracolumns(eq_13_24ref);
%% save results

%%
%pt=eq_13_24ref.point(bifind_13_24(3)); % green 
pt2=eq_1234ref.point(bifind_1234(3));   % purple
%%
pt=eq_13_24ref.point(bifind_13_24(3));
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,20000],pt.x+0.01*rand(length(pt.x),1),opts);
%%
figure(88)
clf;
tiledlayout(1,3)
nexttile
hold on
plot(simulode45.x,simulode45.y(1:4,:),'LineWidth',2)
nexttile 
hold on
plot(simulode45.x,simulode45.y(5:8,:),'LineWidth',2)
nexttile
hold on
plot(simulode45.x,simulode45.y(9:12,:),'LineWidth',2)
%%
save('steady_state_y2_0.mat')


%%
pt=eq_13_24ref.point(200);
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'Events', @(t,x) MyEvent(t,x));
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(2)
simulode45=ode45(Fode45,[0,1e5],pt.x+0.05*rand(length(pt.x),1),opts);
%
figure(88)
clf;
tiledlayout(1,3)
nexttile
hold on
plot(simulode45.x,simulode45.y(1:4,:),'LineWidth',2)
nexttile 
hold on
plot(simulode45.x,simulode45.y(5:8,:),'LineWidth',2)
nexttile
hold on
plot(simulode45.x,simulode45.y(9:12,:),'LineWidth',2)
%%
figure(99)
clf; 
ax1 = axes;
hold(ax1,'on')
plot(ax1,simulode45.x,simulode45.y(5,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(7,:),'.','LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(6,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(8,:),'.','LineWidth',2)
xZoom = 10000*[6.5 7.3];
yZoom = [0.452  0.460];

ax2 = axes('Position',[0.22 0.55 0.38 0.28]);   
box(ax2,'on')
hold(ax2,'on')
plot(ax2,simulode45.x,simulode45.y(5:8,:),'LineWidth',2)
grid(ax2,'off')
ax2.XTickLabel = [];
ax2.YTickLabel = [];
xlim(ax2,xZoom)
ylim(ax2,yZoom)
legend(ax1,'$z_{1}$','$z_{2}$','$z_{3}$','$z_{4}$','Interpreter','latex','FontSize',20,'location','best')
ax2.Color='w';
%%
figure(100)
clf; 
ax1 = axes;
hold(ax1,'on')
plot(ax1,simulode45.x,simulode45.y(1,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(3,:),'.','LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(2,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(4,:),'.','LineWidth',2)
ylim(ax1,[0,0.002])
xZoom = 10000*[6.5 7.3];
yZoom = [0.0  0.002];

ax2 = axes('Position',[0.22 0.55 0.38 0.28]);   
box(ax2,'on')
hold(ax2,'on')
plot(ax2,simulode45.x,simulode45.y(1:4,:),'LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(3,:),'.','LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(2,:),'LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(4,:),'.','LineWidth',2)%plot(ax2,simulode45.x,simulode45.y(1:4,:),'LineWidth',2)
grid(ax2,'off')
ax2.XTickLabel = [];
ax2.YTickLabel = [];
xlim(ax2,xZoom)
ylim(ax2,yZoom)
legend(ax1,'$y_{1}$','$y_{2}$','$y_{3}$','$y_{4}$','Interpreter','latex','FontSize',20,'location','best')
ax2.Color='k'
%%
figure(101)
clf; 
ax1 = axes;
hold(ax1,'on')
plot(ax1,simulode45.x,simulode45.y(9,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(11,:),'.','LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(10,:),'LineWidth',2)
plot(ax1,simulode45.x,simulode45.y(12,:),'.','LineWidth',2)
%ylim(ax1,[0,0.002])
xZoom = 10000*[6.5 7.3];
yZoom = [0.6  0.65];

ax2 = axes('Position',[0.22 0.55 0.38 0.28]);   
box(ax2,'on')
hold(ax2,'on')
plot(ax2,simulode45.x,simulode45.y(9:12,:),'LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(11,:),'.','LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(10,:),'LineWidth',2)
%plot(ax2,simulode45.x,simulode45.y(12,:),'.','LineWidth',2)%plot(ax2,simulode45.x,simulode45.y(1:4,:),'LineWidth',2)
grid(ax2,'off')
ax2.XTickLabel = [];
ax2.YTickLabel = [];
xlim(ax2,xZoom)
ylim(ax2,yZoom)
legend(ax1,'$w_{1}$','$w_{2}$','$w_{3}$','$w_{4}$','Interpreter','latex','FontSize',20,'location','best')
ax2.Color='k'
%%
%%
% cndv_13=[permfix('stst','p13',nstrains,'v','perm',[1,3])];
% cnd_13=[permfix('stst','p13',nstrains,'x','perm',[1,3])];
% y1eq0=dde_sys_cond_create('name','y1_eq0','reference',false,'fun',@(x)x,'args',{'x',1});
% sbargs={'SetupFold.usercond',[cndv_13,y1veq0],'SetupFold.initcond',[cndv_13_24,y1veq0]};
% [fun_eq13,eq_13,suc_13]=SetupStstFrom_fold(fun_eq13_24,eq_13_24ref,bifind_13_24(3),...
%     'initcond',[cndv_13,cnd_13],'usercond',[cndv_13,cnd_13],'outputfuncs',true,...
%     'step',0.1,'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)})
% %fprintf('Branch off from symmetry breaking, keeping (13)(24) symmetry\n')
% figure(1)
% hold on
% eq_13=br_contn(fun_eq13,eq_13,100);
% eq_13=br_rvers(eq_13);
% eq_13=br_contn(fun_eq13,eq_13,100);
% %%
% [eq_13ref,eq_13tests,bifs_13,bifind_13]=...
%     MonitorChange(fun_eq13,eq_13,'min_iterations',5);
% %%
% xpar_13=arrayfun(@(x)x.parameter(in.r),eq_13ref.point);
% ypar_13=arrayfun(@(p)sum(p.x(end-3:end)),eq_13ref.point);
% plt1_13=plot(xpar_13(eq_13tests==0), ypar_13(eq_13tests==0), '.','Color',clrs(4,:),'Linewidth',3,'MarkerSize',10);
% plt2_13=plot(xpar_13(eq_13tests>=1), ypar_13(eq_13tests>=1), '.','Color',clrs(5,:),'Linewidth',2,'MarkerSize',10);
%% 