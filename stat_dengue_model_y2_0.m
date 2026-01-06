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
    'max_bound',   [in.gamma, 1;     in.r, 3],...
    'max_step',    [in.gamma, 0.05;   in.r, 0.05; 0,0.05]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% Tracking a branch of equilirium points
% We impose that all components are 0
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
[fun_1234,eq_1234_0,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.05,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',alleq0,'outputfuncs',true);
figure(1)
clf; hold on
eq_1234_0=br_contn(fun_1234,eq_1234_0,1000);
eq_1234_0=br_rvers(eq_1234_0);
eq_1234_0=br_contn(fun_1234,eq_1234_0,300);
[eq_1234_0,nunst_st]=br_stabl(fun_1234,eq_1234_0,0,0);
%% plot stability
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234_0.point);
ypar_1234=arrayfun(@(x)x.x(1),eq_1234_0.point);
biflocation=find(diff(nunst_st));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%%
% Extract stability information and plot results
figure(1);
clf
hold on 
plot(xpar_1234(nunst_st==0), ypar_1234(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(xpar_1234(nunst_st==1), ypar_1234(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_1234(nunst_st>=2), ypar_1234(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
plot(xpar_1234(biflocation), ypar_1234(biflocation), 'r.','Linewidth',2,'MarkerSize',40);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
eq_1234_0=br_remove_extracolumns(eq_1234_0);
%save('computation_without_delay.mat');
%%
%% Invariant subspaces with respect to the group D_4:
%(1) Fully symetric invariant subspace <(1,1,1,1)> invariant under the
%action of the subgroup <r=(1234)>:={e,r,r^2,r^3}
% odd invariant subspace <(1,-1,1,-1)>  invariant under the
%action of the subgroup :={e,(13)(24),(13),(24)}, expected to detect Hopf
%bifurcation
check_sb=@(cnd)SetupFold(fun_1234,eq_1234_0,biflocation(1),...
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
[fun_eq12_34,eq_12_34,suc_12_34]=SetupStstFrom_fold(funcs,eq_1234_0,biflocation,...
    'usercond',[cnd_12_34,y2eq0],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (12)(34) symmetry\n');
%%
figure(1)
hold on
eq_12_34=br_contn(fun_eq12_34,eq_12_34,30);
eq_12_34=br_rvers(eq_12_34);
eq_12_34=br_contn(fun_eq12_34,eq_12_34,30);
eq_12_34=br_rvers(eq_12_34);
[eq_12_34,nunst_12_34]=br_stabl(fun_eq12_34,eq_12_34);
%% Branch off in three different directions: 2nd 1=3, 2=4, y2=0
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
cndv_13_24=[permfix('stst','p13',nstrains,'v','perm',[1,3]),...
           permfix('stst','p24',nstrains,'v','perm',[2,4])];
cnd_13_24=[permfix('stst','p13',nstrains,'x','perm',[1,3]),...
           permfix('stst','p24',nstrains,'x','perm',[2,4])];
y2eq0=dde_sys_cond_create('name','y2_eq0','reference',false,'fun',@(x)x,'args',{'x',2});
sbargs={'SetupFold.usercond',[alleq0,cndv_13_24,y2veq0],'SetupFold.initcond',[alleq0,cndv_13_24,y2veq0]};
[fun_eq13_24,eq_13_24,suc_13_24]=SetupStstFrom_fold(funcs,eq_1234_0,biflocation,...
    'usercond',[cnd_13_24,y2eq0],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (13)(24) symmetry\n');
%%
figure(1)
hold on
eq_13_24=br_contn(fun_eq13_24,eq_13_24,30);
eq_13_24=br_rvers(eq_13_24);
eq_13_24=br_contn(fun_eq13_24,eq_13_24,30);
eq_13_24=br_rvers(eq_13_24);
[eq_13_24,nunst_13_24]=br_stabl(fun_eq13_24,eq_13_24);
%% Branch off in three different directions: 3rd 1=2=3=4
ex_1234=permfix('stst','p1234',ndim,'x','perm',[1,2,3,4]);
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
sbargs={'SetupFold.usercond',[alleq0,ev_1234],'SetupFold.initcond',[alleq0,ev_1234]};
[fun_eq1234,eq_1234,suc_1234]=SetupStstFrom_fold(funcs,eq_1234_0,biflocation,...
    'usercond',ex_1234,'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)});
fprintf('Branch off from symmetry breaking, keeping (1234) symmetry\n');
%%
figure(1)
hold on
eq_1234=br_contn(fun_eq1234,eq_1234,30);
eq_1234=br_rvers(eq_1234);
eq_1234=br_contn(fun_eq1234,eq_1234,30);
eq_1234=br_rvers(eq_1234);
[eq_1234,nunst_1234]=br_stabl(fun_eq1234,eq_1234);
%%
figure(2);clf
plotopts={'y',@(p)sum(p.x(end-3:end))};
Plot2dBranch(eq_1234_0,plotopts{:});
hold on;
Plot2dBranch(eq_1234,plotopts{:})
Plot2dBranch(eq_12_34,plotopts{:})
Plot2dBranch(eq_13_24,plotopts{:})
xlabel('x');ylabel('sum(w)')