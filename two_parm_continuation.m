clear;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(fullfile('ddebiftool_snapshot'));
format compact
format short g
load('steady_state_y2_0.mat')
%%


%%
check_sb2=@(cnd)SetupFold(fun_eq13_24,eq_13_24ref,bifind_13_24(3),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);

%%
%% Possible Subgroub of the group D_{4}, here we did not find an invariant subspace...
% with repect to any of the possible subgroups of D4
% {e,(13)(24)}
ev13_24=[permfix('stst','p13',ndim,'v','perm',[1,3]),...
         permfix('stst','p24',ndim,'v','perm',[2,4])];
sv_13_24=check_sb2(ev13_24)'
% {e,(14)(23)}
ev14_23=[permfix('stst','p14',ndim,'v','perm',[1,4]),...
         permfix('stst','p23',ndim,'v','perm',[2,3])];
sv_14_23=check_sb2(ev14_23)'
% {e,(12)(34)}
ev12_34=[permfix('stst','p12',ndim,'v','perm',[1,2]),...
         permfix('stst','p34',ndim,'v','perm',[3,4])];
sv_12_34=check_sb2(ev12_34)'
% {e,(24)}
ev24=permfix('stst','p24',ndim,'v','perm',[2,4]);
sv_24=check_sb(ev24)'
% {e,(13)}
ev13=permfix('stst','p13',ndim,'v','perm',[1,3]);
[sv_13,se_13]=check_sb2(ev13)
%
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
sv_1234=check_sb2(ev_1234)'
% add additional constraints where y=(z=)0 for some component
%y2veq0=dde_sys_cond_create('name','y2_eq0','reference',false,'fun',@(v)v,'args',{'x',3});

%%
v2veq0=dde_sys_cond_create('name','y3_eq0','reference',false,'fun',@(v)v,'args',{'v',2});
sv_13_24_y0=check_sb2([ev13_24,v2veq0])';
sv_12_34_y0=check_sb2([ev12_34,y2veq0])';
sv_13_y0=check_sb2([ev13,v2veq0]);
%%
fold_constrains=[ev13,v2veq0];%[ev13_24,y2veq0];
 parbd_fld={'min_bound',[in.gamma, -2;  in.r,-4],...
    'max_bound',   [in.gamma, 4;     in.r, 7],...
    'max_step',    [in.gamma, 0.05;   in.r, 0.05; 0,0.05]};
[funsfold_13,fld_branch_13,suc]=SetupFold(fun_eq13_24,eq_13_24ref,bifind_13_24(3),...
    'initcond',fold_constrains,'usercond',fold_constrains,'outputfuncs',true,...
    'contpar',[in.r,in.gamma],'dir',in.r,'step',0.1,'print_residual_info',1,parbd_fld{:});
%%
figure(9)
clf
hold on
fld_branch_13=br_contn(funsfold_13,fld_branch_13,500);
fld_branch_13=br_rvers(fld_branch_13);
fld_branch_13=br_contn(funsfold_13,fld_branch_13,400);
%%
[fld_branch_13ref,fold_nunst,fold_bifs,fold_bifind]=MonitorChange(funsfold_13,fld_branch_13,...
    'printlevel',1,'print_residual_info',0,'min_iterations',5,'output','branch',...
    'locate_trivial',@(p)[0 0]);
%%
yparf=arrayfun(@(p)sum(p.x(end-3:end)),fld_branch_13ref.point);
clrs=lines();
xpar_fld=arrayfun(@(x)x.parameter(in.r),fld_branch_13ref.point);
ypar_fld=arrayfun(@(x)x.parameter(in.gamma),fld_branch_13ref.point);
figure(11)
clf;
hold on
plot(xpar_fld(fold_nunst==0),ypar_fld(fold_nunst==0),'--','Color','r','LineWidth',2)
plot(xpar_fld(fold_nunst>=1),ypar_fld(fold_nunst>=1),'--','Color','g','LineWidth',2)
%% simulation near the fold-fold bif
pt=eq_13_24ref.point(bifind_13_24(3));
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,1e6],pt.x+0.01*rand(length(pt.x),1),opts);
%%
figure(88)
clf;
tiledlayout(1,3)
nexttile
hold on
plot(simulode45.x,simulode45.y(1,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(3,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(2,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(4,:),'LineWidth',2)
ylabel('$y_{i}$','FontSize',16,'Interpreter','latex')
nexttile 
hold on
plot(simulode45.x,simulode45.y(5,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(7,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(6,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(8,:),'LineWidth',2)

ylabel('$z_{i}$','FontSize',16,'Interpreter','latex')

nexttile
hold on
plot(simulode45.x,simulode45.y(9,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(11,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(10,:),'LineWidth',2)
plot(simulode45.x,simulode45.y(12,:),'LineWidth',2)


ylabel('$w_{i}$','FontSize',16,'Interpreter','latex')
%%
save('two-parameter_computation.mat')