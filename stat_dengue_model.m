clear;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% ddebiftool_path(base,'sym');
  base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
addpath([base,filesep,'ddebiftool'],...
    [base,filesep,'ddebiftool_extra_psol'],...
    [base,filesep,'ddebiftool_utilities'],...
    [base,filesep,'ddebiftool_extra_rotsym'],...
    [base,filesep,'ddebiftool_extra_nmfm'],...
    [base,filesep,'ddebiftool_extra_symbolic'],...
    [base,filesep,'ddebiftool_coco']);
format compact
format short 
%%
parnames = {'r',  'e',  'beta',   'mu',   'sigma',   'gamma',  'tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.r,      in.e,     in.beta,        in.mu,      in.sigma,     in.gamma,  in.tau_delay])=... 
    [0.0,         1,        40,           0.02,            10,         0.54         0  ];
ndim=12;
x0=zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.r,0.0],...
    'max_bound',   [in.gamma, 1;     in.r, 3],...
    'max_step',    [in.gamma, 0.05;   in.r, 0.05; 0,0.05]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% Tracking a branch of equilirium points
% We impose symmetry with respect to the permutation (1234) to avoid
% singularity along the branch
ex1234=permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]);
[fun_1234,eq_1234_0,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.05,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
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
%% Possible Subgroub of the group D_{4}, 
% here we did not find an invariant subspace ...
% with repect to any of the possible subgroups of D4 as we can see in the
% following:
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
sv_13=check_sb(ev13)'
%
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
vs_1234=check_sb(ev_1234)'
% 
%% simulation with dde23
% because no an invariant subspace is found with respect to any subgroup of
% $D_4$, we run ode45 and then use the result to compute one-parameter
% branch in $r$
pt=eq_1234_0.point(biflocation-1);
%
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,20000],pt.x+0.01*rand(length(pt.x),1));
%
figure(88)
clf;
hold on
plot(simulode45.x,simulode45.y(9:12,:),'LineWidth',2)
%%
ex2_d45=[permfix('stst','eq2_d45',ndim,'x','perm',[1,3])]%;...
        % permfix('stst','eq2_d45',ndim,'x','perm',[2,4])];
x1=simulode45.y(:,end);
[fun2_d45,eq2_d45_13_24,suc2]=SetupStst(funcs,...
    'parameter',pt.parameter,'x',x1,...
    'contpar',in.r,'step',0.05,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex2_d45,'outputfuncs',true);
%%
figure(1)
hold on
eq2_d45_13_24=br_contn(fun2_d45,eq2_d45_13_24,300);
eq2_d45_13_24=br_rvers(eq2_d45_13_24);
eq2_d45_13_24=br_contn(fun2_d45,eq2_d45_13_24,200);
%
[eq2_d45_13_24,nunst2_st]=br_stabl(fun2_d45,eq2_d45_13_24,0,0);
biflocation2=find(diff(nunst2_st));
xpar_13_24=arrayfun(@(p)p.parameter(in.r),eq2_d45_13_24.point);
ypar_13_24=arrayfun(@(x)x.x(1),eq2_d45_13_24.point);
%%
figure(11);clf;ax11=gca;hold(ax11,'on');grid(ax11,'on');set(ax11,'fontsize',16);
args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax11};
Plot2dBranch(eq2_d45_13_24,'DisplayName','(13)(24)_ic_ode45','color',[0.5,0,0],args{:});
Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.3,0],args{:});
xlabel(ax11,'r');ylabel(ax11,'y_1');
plot(xpar_13_24(biflocation2),ypar_13_24(biflocation2),'.','MarkerSize',20)
plot(xpar_1234(biflocation), ypar_1234(biflocation),'.','MarkerSize',20)
%% Now we check...
% here bifurcation near r=1
check_sb_13_24=@(cnd)SetupFold(fun2_d45,eq2_d45_13_24,biflocation2(2),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true); 
% Possible Subgroub of the group D_{4}, here we did not find an invariant subspace...
sv_13_24=check_sb_13_24(ev13_24)'
% {e,(14)(23)}
sv_14_23=check_sb_13_24(ev14_23)'
% {e,(12)(34)}
sv_12_34=check_sb_13_24(ev12_34)'
% {e,(24)}
sv_24=check_sb_13_24(ev24)'
% {e,(13)}
sv_13=check_sb_13_24(ev13)'
%
sv_1234=check_sb(ev_1234)' % only this works
%
%% here bifurcation near r=2.1
check_sb_13_24=@(cnd)SetupFold(fun2_d45,eq2_d45_13_24,biflocation2(1),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true); 
% Possible Subgroub of the group D_{4}, here we did not find an invariant subspace...
sv2_13_24=check_sb_13_24(ev13_24)'
% {e,(14)(23)}
sv2_14_23=check_sb_13_24(ev14_23)'
% {e,(12)(34)}
sv2_12_34=check_sb_13_24(ev12_34)'
% {e,(24)}
sv2_24=check_sb_13_24(ev24)'  %% this works
% {e,(13)}
sv2_13=check_sb_13_24(ev13)'
%
sv2_1234=check_sb(ev_1234)' %  this works
% %%%%%%%%%%%%%%%%%%%%%%%
eq2_d45_13_24=br_remove_extracolumns(eq2_d45_13_24);
%% Branching off in different directions
%On our branch solutions satisfy (13)(24)
cnd_sb_13_24=[permfix('stst','eq2_d45',ndim,'x','perm',[2,4]),...
              permfix('stst','eq2_d45',ndim,'v','perm',[2,4])];
cnd_24=permfix('stst','p24',ndim,'x','perm',[2,4]);
sbargs={'SetupFold.usercond',cnd_sb_13_24,'SetupFold.initcond',cnd_sb_13_24};
[fun_eq24,eq_24,suc_24]=SetupStstFrom_fold(fun2_d45,eq2_d45_13_24,biflocation2(1),...
    'usercond',cnd_24,'outputfuncs',true,...
    'step',0.1,sbargs{:});
figure(1)
%clf; 
% hold on
eq_24=br_contn(fun_eq24,eq_24,300);
eq_24=br_rvers(eq_24);
eq_24=br_contn(fun_eq24,eq_24,300);
[eq_24,nunst_24]=br_stabl(fun_eq24,eq_24,0,0);
biflocation_24=find(diff(nunst_24))
%%
cnd_sb_13_24=[permfix('stst','eq2_d45',ndim,'x','perm',[1,3]),...
              permfix('stst','eq2_d45',ndim,'x','perm',[2,4]),...
              permfix('stst','eq2_d45',ndim,'v','perm',[1,2,3,4])];
cnd_1234=permfix('stst','p24',ndim,'x','perm',[1,2,3,4]);
sbargs={'SetupFold.usercond',cnd_sb_13_24,'SetupFold.initcond',cnd_sb_13_24};
[fun_eq1234,eq2_1234,suc_1234]=SetupStstFrom_fold(fun2_d45,eq2_d45_13_24,biflocation2(2),...
    'usercond',cnd_1234,'outputfuncs',true,...
    'step',0.1,sbargs{:});
figure(1)
%clf; 
% hold on
eq2_1234=br_contn(fun_eq1234,eq2_1234,400);
eq2_1234=br_rvers(eq2_1234);
eq2_1234=br_contn(fun_eq1234,eq2_1234,100);
[eq2_1234,nunst_1234]=br_stabl(fun_eq1234,eq2_1234,0,0);
biflocation_1234=find(diff(nunst_1234));
%%
figure(3);clf;ax3=gca;hold(ax3,'on');grid(ax3,'on');set(ax3,'fontsize',16);
args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax3};
Plot2dBranch(eq2_d45_13_24,'DisplayName','(13)(24)_ic_ode45','color',[0.2,0,0],args{:});
Plot2dBranch(eq_24,'DisplayName','(24)_ic_ode45','color',[0.,0,0.0],args{:});
Plot2dBranch(eq2_1234,'DisplayName','(1234)_br2','color',[0,0.3,0],args{:});
Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0.3,0.3,0.3],args{:});
plot(xpar_13_24(biflocation2),ypar_13_24(biflocation2),'.','MarkerSize',20)
plot(xpar_1234(biflocation), ypar_1234(biflocation),'.','MarkerSize',20)
xlabel(ax3,'r');ylabel(ax3,'y_1');
%%
eq2_1234=br_remove_extracolumns(eq2_1234);
eq_24=br_remove_extracolumns(eq_24);
eq2_d45_13_24=br_remove_extracolumns(eq2_d45_13_24);
save('demo_multStrain.mat')
% clear;
% % base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% % ddebiftool_path(base,'sym');
%   base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% addpath([base,filesep,'ddebiftool'],...
%     [base,filesep,'ddebiftool_extra_psol'],...
%     [base,filesep,'ddebiftool_utilities'],...
%     [base,filesep,'ddebiftool_extra_rotsym'],...
%     [base,filesep,'ddebiftool_extra_nmfm'],...
%     [base,filesep,'ddebiftool_extra_symbolic'],...
%     [base,filesep,'ddebiftool_coco']);
% format compact
% format short 
% %%
% parnames = {'r',  'e',  'beta',   'mu',   'sigma',   'gamma',  'tau_delay'};
% cind = [parnames; num2cell(1:length(parnames))];
% in = struct(cind{:});
% par([in.r,      in.e,     in.beta,        in.mu,      in.sigma,     in.gamma,  in.tau_delay])=... 
%     [0.0,         1,        40,           0.02,            10,         0.54         0  ];
% ndim=12;
% x0=zeros(ndim,1);
%     parbd={'min_bound',[in.gamma, 0;  in.r,0],...
%     'max_bound',   [in.gamma, 1;     in.r, 3],...
%     'max_step',    [in.gamma, 0.01;   in.r, 0.01; 0,0.05]};
% %%
% % par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
% par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
% funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
% %% Tracking a branch of equilirium points
% % We impose symmetry with respect to the permutation (1234) to avoid
% % singularity along the branch
% ex1234=permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]);
% [fun_1234,eq_1234_0,suc0]=SetupStst(funcs,...
%     'parameter',par,'x',x0,...
%     'contpar',in.r,'step',0.05,parbd{:},...
%     'extracolumns','auto','extra_condition',true,...
%     'usercond',ex1234,'outputfuncs',true);
% figure(1)
% clf; hold on
% eq_1234_0=br_contn(fun_1234,eq_1234_0,1000);
% eq_1234_0=br_rvers(eq_1234_0);
% eq_1234_0=br_contn(fun_1234,eq_1234_0,300);
% %%
% [eq_1234_0,nunst_st]=br_stabl(fun_1234,eq_1234_0,0,0);
% %% plot stability
% r_stst_par=arrayfun(@(x)x.parameter(in.r),eq_1234_0.point);
% x1_axis=arrayfun(@(x)x.x(1),eq_1234_0.point);
% bif_point=find((nunst_st>0),1,'last');
% p1=eq_1234_0.point(bif_point);
% %%
% % Extract stability information and plot results
% figure(1);
% clf
% hold on 
% plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
% plot(r_stst_par(nunst_st==1), x1_axis(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
% plot(r_stst_par(nunst_st>=2), x1_axis(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
% plot(r_stst_par(bif_point), x1_axis(bif_point), 'r.','Linewidth',2,'MarkerSize',40);
% xlabel('r');
% ylabel('y_{1}');
% title('Stability Analysis of Equilibria');
% %%
% eq_1234_0=br_remove_extracolumns(eq_1234_0);
% save('computation_without_delay.mat');
% %%
% %% Invariant subspaces with respect to the group D_4:
% %(1) Fully symetric invariant subspace <(1,1,1,1)> invariant under the
% %action of the subgroup <r=(1234)>:={e,r,r^2,r^3}
% % odd invariant subspace <(1,-1,1,-1)>  invariant under the
% %action of the subgroup :={e,(13)(24),(13),(24)}, expected to detect Hopf
% %bifurcation
% check_sb=@(cnd)SetupFold(fun_1234,eq_1234_0,bif_point,...
%     'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% %% Possible Subgroub of the group D_{4}, here we did not find an invariant subspace...
% % with repect to any of the possible subgroups of D4
% % {e,(13)(24)}
% ev13_24=[permfix('stst','p13',ndim,'v','perm',[1,3]),...
%          permfix('stst','p24',ndim,'v','perm',[2,4])];
% sv_13_24=check_sb(ev13_24)'
% % {e,(14)(23)}
% ev14_23=[permfix('stst','p14',ndim,'v','perm',[1,4]),...
%          permfix('stst','p23',ndim,'v','perm',[2,3])];
% sv_14_23=check_sb(ev14_23)'
% % {e,(12)(34)}
% ev12_34=[permfix('stst','p12',ndim,'v','perm',[1,2]),...
%          permfix('stst','p34',ndim,'v','perm',[3,4])];
% sv_12_34=check_sb(ev12_34)'
% % {e,(24)}
% ev24=permfix('stst','p24',ndim,'v','perm',[2,4]);
% sv_24=check_sb(ev24)'
% % {e,(13)}
% ev13=permfix('stst','p13',ndim,'v','perm',[1,3]);
% sv_13=check_sb(ev13)'
% %
% ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
% vs_1234=check_sb(ev_1234)'
% % 
% %% simulation with dde23
% % because no an invariant subspace is found with respect to any subgroup of
% % $D_4$, we run ode45 and then use the result to compute one-parameter
% % branch in $r$
% pt=eq_1234_0.point(bif_point-1);
% %
% Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
% rng(1)
% simulode4523=ode45(Fode45,[0,20000],pt.x+0.01*rand(length(pt.x),1));
% %%
% figure(88)
% clf;
% hold on
% plot(simulode4523.x,simulode4523.y(9:12,:),'LineWidth',2)
% %%
% ex2_d45=[permfix('stst','eq2_d45',ndim,'x','perm',[1,3]);...
%         permfix('stst','eq2_d45',ndim,'x','perm',[2,4])];
% x1=simulode4523.y(:,end);
% [fun2_d45,eq2_d45_13_24,suc2]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',x1,...
%     'contpar',in.r,'step',0.05,parbd{:},...
%     'extracolumns','auto','extra_condition',true,...
%     'usercond',ex2_d45,'outputfuncs',true);
% eq2_d45_13_24.parameter.max_step(:,2)=[0.02;0.02;0.02];
% eq2_d45_13_24.parameter.max_bound(end)=3;
% %%
% figure(1)
% hold on
% eq2_d45_13_24=br_contn(fun2_d45,eq2_d45_13_24,300);
% eq2_d45_13_24=br_rvers(eq2_d45_13_24);
% eq2_d45_13_24=br_contn(fun2_d45,eq2_d45_13_24,50);
% %%
% [eq2_d45_13_24,nunst2_st]=br_stabl(fun2_d45,eq2_d45_13_24,0,0);
% biflocation2=find(diff(nunst2_st));
% %% %%%%%%%%%%%%%%%%%%%%%%%
% ex3_d45=[permfix('stst','eq3_d45',ndim,'x','perm',[1,2,3,4])];
% x1=simulode4523.y(:,end);
% [fun3_d45,eq3_d45_1234,suc3]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',x1,...
%     'contpar',in.r,'step',0.05,parbd{:},...
%     'extracolumns','auto','extra_condition',true,...
%     'usercond',ex3_d45,'outputfuncs',true);
% eq3_d45_1234.parameter.max_step(:,2)=[0.05;0.05;0.05];
% eq3_d45_1234.parameter.max_bound(end)=3;
% %%
% figure(1)
% hold on
% eq3_d45_1234=br_contn(fun3_d45,eq3_d45_1234,300);
% eq3_d45_1234=br_rvers(eq3_d45_1234);
% eq3_d45_1234=br_contn(fun3_d45,eq3_d45_1234,50);
% %%
% [eq3_d45_1234,nunst3_st]=br_stabl(fun3_d45,eq3_d45_1234,0,0);
% 
% biflocation3=find(diff(nunst3_st));
% %%
% ex4_d45=[permfix('stst','eq3_d45',ndim,'x','perm',[2,4])];
% x1=simulode4523.y(:,end);
% [fun4_d45,eq4_d45_24,suc4]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',x1,...
%     'contpar',in.r,'step',0.05,parbd{:},...
%     'extracolumns','auto','extra_condition',true,...
%     'usercond',ex4_d45,'outputfuncs',true);
% eq4_d45_24.parameter.max_step(:,2)=[0.05;0.05;0.05];
% eq4_d45_24.parameter.max_bound(end)=3;
% %%
% figure(1)
% hold on
% eq4_d45_24=br_contn(fun4_d45,eq4_d45_24,200);
% eq4_d45_24=br_rvers(eq4_d45_24)
% eq4_d45_24=br_contn(fun4_d45,eq4_d45_24,200);
% %%
% [eq4_d45_24,nunst4_st]=br_stabl(fun4_d45,eq4_d45_24,0,0);
% 
% biflocation4=find(diff(nunst4_st));
% %% %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%
% ex5_d45=[permfix('stst','eq3_d45',ndim,'x','perm',[1,3])];
% x1=simulode4523.y(:,end);
% [fun5_d45,eq5_d45_13,suc5]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',x1,...
%     'contpar',in.r,'step',0.05,parbd{:},...
%     'extracolumns','auto','extra_condition',true,...
%     'usercond',ex5_d45,'outputfuncs',true);
% eq5_d45_13.parameter.max_step(:,2)=[0.05;0.05;0.05];
% eq5_d45_13.parameter.max_bound(end)=3;
% figure(1)
% hold on
% eq5_d45_13=br_contn(fun5_d45,eq5_d45_13,200);
% eq5_d45_13=br_rvers(eq5_d45_13);
% eq5_d45_13=br_contn(fun5_d45,eq5_d45_13,200);
% [eq5_d45_13,nunst5_st]=br_stabl(fun5_d45,eq5_d45_13,0,0);
% biflocation5=find(diff(nunst5_st));
% 
% 
% %% %%%%%%%%%%%%%%%%
% figure(3);clf;ax3=gca;hold(ax3,'on');grid(ax3,'on');set(ax3,'fontsize',16);
% args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax3};
% Plot2dBranch(eq2_d45_13_24,'DisplayName','(13)(24)_ic_ode45','color',[0.5,0,0],args{:});
% Plot2dBranch(eq3_d45_1234,'DisplayName','(1234)_ic_ode45','color',[0,0.5,0],args{:});
% Plot2dBranch(eq4_d45_24,'DisplayName','(24)_ic_ode45','color',[0,0,0.5],args{:});
% Plot2dBranch(eq5_d45_13,'DisplayName','(13)_ic_ode45','color',[0.3,0,0],args{:});
% Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.3,0],args{:});
% xlabel(ax3,'r');ylabel(ax3,'y_1');
% %%
% eq2_d45_13_24=br_remove_extracolumns(eq2_d45_13_24);
% eq3_d45_1234=br_remove_extracolumns(eq3_d45_1234);
% eq4_d45_24=br_remove_extracolumns(eq4_d45_24);
% eq5_d45_13=br_remove_extracolumns(eq5_d45_13);
% eq_1234_0=br_remove_extracolumns(eq_1234_0);
% save('computation_without_delay_second_part.mat')
% %%
% 
% figure(4);clf;ax4=gca;hold(ax4,'on');grid(ax4,'on');set(ax4,'fontsize',16);
% args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax4};
% Plot2dBranch(eq3_d45_1234,'DisplayName','(1234)_ic_ode45','color',[0,0.5,0],args{:});
% Plot2dBranch(eq5_d45_13,'DisplayName','(13)_ic_ode45','color',[0.3,0,0],args{:});
% Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.3,0],args{:});
% xlabel(ax4,'r');ylabel(ax4,'y_1');
% %
% figure(5);clf;ax5=gca;hold(ax5,'on');grid(ax5,'on');set(ax5,'fontsize',16);
% args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax5};
% Plot2dBranch(eq2_d45_13_24,'DisplayName','(13)(24)_ic_ode45','color',[0.5,0,0],args{:});
% Plot2dBranch(eq3_d45_1234,'DisplayName','(1234)_ic_ode45','color',[0,0.5,0],args{:});
% Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.3,0],args{:});
% xlabel(ax5,'r');ylabel(ax5,'y_1');
% %
% figure(6);clf;ax6=gca;hold(ax6,'on');grid(ax6,'on');set(ax6,'fontsize',16);
% args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax6};
% Plot2dBranch(eq3_d45_1234,'DisplayName','(1234)_ic_ode45','color',[0,0.5,0],args{:});
% Plot2dBranch(eq4_d45_24,'DisplayName','(24)_ic_ode45','color',[0,0,0.5],args{:});
% Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.3,0],args{:});
% xlabel(ax6,'r');ylabel(ax6,'y_1');
% %%
% check_sb2=@(cnd)SetupFold(fun_d45,eq2_d45_13_24,biflocation(1),...
%     'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% sv_13_24_ode45=check_sb2(ev13_24)'
% sv_13_ode45=check_sb2(ev13)'
% sv_24_od45=check_sb2(ev24)'
% vs_1234_ode45=check_sb2(ev_1234)'
% 
% 
% %%
% check_hopf=@(cnd)SetupHopf(funcs,eq2_d45_13_24,biflocation(1),...
%     'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% 
% ev13_0=permfix('stst','p1234_1_4',ndim,'v',...
%     'perm',[1,3],'rotation',[0,1]);
% sv_13_0=check_hopf(ev13_0)'
% % Permutation (12)_0 only is not enough: we have four zero singular values for time shifts
% %0 or 1/2. Time shift 1/4 would be ok
% ev12_0=permfix('stst','p12',ndim,'v',...
%     'perm',[2,4],'rotation',[0,1]);
% sv_12_0=check_hopf(ev12_0)'
% 
% %%
% % %% Branch off toward different periodic solutions
% % hopfsyms={ev13_0,ev12_0};
% % psol_symmetries.p13_0_1      =[permfix('psol','p13_0_1 ',ndim,'x','perm',[1,3],'shift',[0,1]),...
% %                                permfix('psol','p24_0_1 ',ndim,'x','perm',[2,4],'shift',[0,1])];
% % %
% % psol_symmetries.p24_0_1       =permfix('psol','p24_0_1', ndim,'x','perm',[2,4],  'shift',[0,1]);
% % %
% % 
% % symnames=fieldnames(psol_symmetries);
% % psymrealname.p13_0_1='p$(13)_{0}$';
% % psymrealname.p24_0_1='p$(24)_0$';
% % clear sucpsol
% % psargs={'extracolumns','auto','extra_condition',true,... % important for imposing symmetry
% %     'matrix','sparse','eigmatrix','sparse','degree',4,'intervals',60,... % determines accuracy
% %     'max_step',[0,0.1]}; % stepsize along branch
% % figure(1);
% % hold on
% % for i=1:length(symnames)
% %     sm=symnames{i};
% %     fprintf('Symmetry: %s\n',sm);
% %     hopfargs={'SetupHopf.usercond',hopfsyms{i},'SetupHopf.initcond',hopfsyms{i}};
% %     [funpsol.(sm),psol.(sm),sucpsol.(sm)]=SetupPsol(funcs,eq2_d45_13_24,biflocation(1),...
% %         'outputfuncs',true,hopfargs{:},psargs{:},...
% %         'usercond',psol_symmetries.(sm),'initcond',psol_symmetries.(sm),...
% %         'plot_measure',{@(p)p.parameter(in.r),@(p)max(mean(p.profile,1))});
% %     disp(sucpsol.(sm))
% %     psol.(sm)=br_contn(funpsol.(sm),psol.(sm),10);
% %     [psolref.(sm),psol_tests.(sm),psolbifs.(sm),psolbifind.(sm)]=...
% %         MonitorChange(funpsol.(sm),psol.(sm),'min_iterations',5,'range',2:length(psol.(sm).point));
% % end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Branching off in different directions
% % On our branch solutions satisfy (13)(24)
% % 
% % cnd_sb=[permfix('stst','eq2_d45',ndim,'x','perm',[1,3]),...
% %         permfix('stst','eq2_d45',ndim,'x','perm',[2,4]),...
% %         permfix('stst','pv13',ndim,'v','perm',[1,3])];
% % cnd_13=permfix('stst','px13',ndim,'x','perm',[1,3]);
% % sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
% % [fun_eq13,eq_13,suc9]=SetupStstFrom_fold(funcs,eq2_d45_13_24,biflocation(1),...
% %     'usercond',cnd_13,'outputfuncs',true,...
% %     'step',0.1,sbargs{:})
% % figure(1)
% % %clf; 
% % % hold on
% % eq_13=br_contn(fun_eq13,eq_13,40);
% % eq_13=br_rvers(eq_13)
% % eq_13=br_contn(fun_eq13,eq_13,40);
% % [eq_13,nunst3_st]=br_stabl(fun_eq13,eq_13,0,0);
% % %%
% % cnd_sb_13_24=[permfix('stst','eq2_d45',ndim,'x','perm',[1,3]),...
% %         permfix('stst','eq2_d45',ndim,'x','perm',[2,4]),...
% %         permfix('stst','pv13',ndim,'v','perm',[1,3]),...
% %         permfix('stst','pv24',ndim,'v','perm',[2,4])];
% % cnd_13_24=[permfix('stst','px13',ndim,'x','perm',[1,3]),...
% %     permfix('stst','px24',ndim,'x','perm',[2,4])];
% % sbargs={'SetupFold.usercond',cnd_sb_13_24,'SetupFold.initcond',cnd_sb_13_24};
% % [fun_eq13_24,eq_13_24,suc]=SetupStstFrom_fold(funcs,eq2_d45_13_24,biflocation(1),...
% %     'usercond',cnd_13_24,'outputfuncs',true,...
% %     'step',0.1,sbargs{:})
% % %%
% % figure(1)
% % %clf; 
% % % hold on
% % eq_13_24=br_contn(fun_eq13_24,eq_13_24,40);
% % eq_13_24=br_rvers(eq_13_24)
% % eq_13_24=br_contn(fun_eq13_24,eq_13_24,40);
% % [eq_13_24,nunst13_24]=br_stabl(fun_eq13_24,eq_13_24,0,0);
% %% Continue secondary symmetry breaking in 2 parameters
% % cnd_sb0=[permfix('stst','px13',ndim,'x','perm',[1,2,3,4]);...
% %     permfix('stst','pv13',ndim,'v','perm',[1,3]); ...
% %     permfix('stst','pv13',ndim,'v','perm',[2,4])];
% % [fun_sb2,eq_sb2,suc]=SetupFold(funcs,eq_1234_0,bif_point(1),...
% %     'initcond',cnd_sb0,'usercond',cnd_sb0,'outputfuncs',true,...
% %     'contpar',[in.r,in.gamma],'dir',in.r,'step',0.1,'print_residual_info',1);
% % fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
% % disp(eq_sb2)
% % %%
% % figure(2);
% % clf; hold on
% % eq_sb2=br_contn(fun_sb2,eq_sb2,100);
% 
% %%
% 
% %%
% 
% % 
% % %%
% % figure(3);clf;ax3=gca;hold(ax3,'on');%grid(ax3,'on');
% % args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax3};
% % Plot2dBranch(eq_13_24ref,'DisplayName','(13)','color',[0,0,0.4],args{:});
% % Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.4,0],args{:});
% % xlabel(ax3,'r');ylabel(ax3,'y_1');
% % %%
% % eq_1234_0=br_remove_extracolumns(eq_1234_0);
% % eq_13_24ref=br_remove_extracolumns(eq_13_24ref);
% % eqbranch_13_24=br_remove_extracolumns(eqbranch_13_24);
% % save('stat_dengue_model.mat');
% % %%
% % %% Continue secondary symmetry breaking in 2 parameters
% % cnd_sb2=[permfix('stst','px13',ndim,'x','perm',[1,2,3,4]);...
% %     permfix('stst','pv13',ndim,'v','perm',[1,3]); ...
% %     permfix('stst','pv13',ndim,'v','perm',[2,4])];
% % [fun_sb2,eq_sb2,suc]=SetupFold(funcs,eq_13_24ref,eq_123bifind(3),...
% %     'initcond',cnd_sb2,'usercond',cnd_sb2,'outputfuncs',true,...
% %     'contpar',[in.r,in.gamma],'dir',in.r,'step',0.1,'print_residual_info',1);
% % fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
% % disp(eq_sb2)
% % %%
% % figure(2);
% % %clf; hold on
% % eq_sb2=br_contn(fun_sb2,eq_sb2,100);
% % %%
% % eq_sb2=br_rvers(eq_sb2);
% % eq_sb2=br_contn(fun_sb2,eq_sb2,100);
% % %%
% % [eq_sb2ref,eq_sb2tests,eq_sb2bifs,eq_sb2bifind]=MonitorChange(fun_sb2,eq_sb2,...
% %     'printlevel',1,'print_residual_info',0,'min_iterations',5,'output','branch');
% % %% Continue symmetry breaking in 2 parameters
% % cnd_sb=[permfix('stst','px1234',ndim,'x','perm',[1,3]);...
% %     permfix('stst','pv123',ndim,'v','perm',[1,3])];
% % [fun_sb,eq_sb,suc]=SetupFold(funcs,eq_13_24ref,eq_123bifind(1)+1,...
% %     'initcond',cnd_sb,'usercond',cnd_sb,'outputfuncs',true,...
% %     'contpar',[in.r,in.beta],'dir',in.r,'step',0.1,'print_residual_info',1);
% % fprintf('Initialize two-parameter continuation of  symmetry breaking\n');
% % disp(eq_sb)
% % figure(2);clf;ax2=gca;
% % xlabel(ax2,'a');ylabel(ax2,'c_{ext}');
% % eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
% % eq_sb=br_rvers(eq_sb);
% % eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
% % [eq_sbref,eq_sbtests,eq_sbifs,eq_sbbifind]=MonitorChange(fun_sb,eq_sb,...
% %     'printlevel',1,'print_residual_info',0,'min_iterations',5,'output','branch',...
% %     'locate_trivial',@(p)[0,0,0]);
% 
% % cnd_sb2=[permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]),...
% %         permfix('stst','pv13',ndim,'v','perm',[1,2,3,4])];
% % cnd_1234=permfix('stst','px1234',ndim,'x','perm',[1,2,3,4]);
% % 
% % %%
% % sbargs2={'SetupFold.usercond',cnd_sb2,'SetupFold.initcond',cnd_sb2};
% % [fun_eq1234_bk,eq1234_bk,suc]=SetupStstFrom_fold(funcs,eq_1234_0,bif_point(1),...
% %     'usercond',cnd_1234,'outputfuncs',true,...
% %     'step',0.1,sbargs2{:});
% % fprintf('Branch off from symmetry breaking, keeping (1234) symmetry\n');
% % disp(eq1234_bk);
% % figure(1);
% % eq1234_bk=br_contn(fun_eq1234_bk,eq1234_bk,200);
% % eq1234_bk=br_rvers(eq1234_bk);
% % eq1234_bk=br_contn(fun_eq1234_bk,eq1234_bk,100);
% %%
% % [branch,indbif,indmap,notcorrected]=br_bisection(funcs,branch,inds,'fold')
% 
% %cnd_sb=[permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]),...
% %         permfix('stst','pv13',ndim,'v','perm',[1,2,3,4])];
% % cnd_13=permfix('stst','px13_24',ndim,'x','perm',[1,2,3,4]);
% % sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
% % [fun_eq12_34,eqbranch_13_24,suc]=SetupStstFrom_fold(funcs,eq_1234_0,bif_point(1),...
% %     'usercond',[],'outputfuncs',true,...
% %     'step',0.1,sbargs{:});
% % %fprintf('Branch off from symmetry breaking, keeping (13) symmetry\n');
% % %%
% % eqbranch_13_24.parameter.max_bound(end)=7;
% % disp(eqbranch_13_24);
% % figure(1);
% % eqbranch_13_24=br_contn(fun_eq12_34,eqbranch_13_24,30);
% % eqbranch_13_24=br_rvers(eqbranch_13_24);
% % eqbranch_13_24=br_contn(fun_eq12_34,eqbranch_13_24,200);
% % %% compute stability
% % [eq_13_24ref,eq_123tests,eq_123bifs,eq_123bifind]=...
% %     MonitorChange(fun_eq12_34,eqbranch_13_24,'print_residual_info',0,'min_iterations',5);
% % %% plot stability
% % r_stst_par2=arrayfun(@(x)x.parameter(in.r),eq_13_24ref.point);
% % x1_axis2=arrayfun(@(x)x.x(1),eq_13_24ref.point);
% % 
% % %% Extract stability information and plot results
% % p2=eq_13_24ref.point(eq_123bifind(1));
% % p3=eq_13_24ref.point(eq_123bifind(2));
% % p4=eq_13_24ref.point(eq_123bifind(3));
% % figure(11);
% % clf
% % hold on 
% % plot(r_stst_par2(eq_123tests==0), x1_axis2(eq_123tests==0), 'kx','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par2(eq_123tests==1), x1_axis2(eq_123tests==1), 'r.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par2(eq_123tests==2), x1_axis2(eq_123tests==2), 'g.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par2(eq_123tests>=3), x1_axis2(eq_123tests>=3), 'b.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par(nunst_st==1), x1_axis(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par(nunst_st>=2), x1_axis(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
% % plot(r_stst_par(bif_point), x1_axis(bif_point), 'r.','Linewidth',2,'MarkerSize',40);
% % %plot(r_stst_par2(eq_123bifind(1)),x1_axis2(eq_123bifind(1)),'k.','Linewidth',2,'MarkerSize',20);
% % xlabel('r');
% % ylabel('y_{1}');
% % title('Stability Analysis of Equilibria');
% % %%
% % save('computation_without_delay2.mat');
% % %%
% %%
% % cnd_sb=[permfix('stst','eq2_1234_0',ndim,'x','perm',[1,2,3,4])
% %         permfix('stst','pv13',ndim,'v','perm',[1,3])];
% % cnd_13=[permfix('stst','px13',ndim,'x','perm',[1,3])];
% % sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
% % [fun_eq13,eq_134,suc]=SetupStstFrom_fold(funcs,eq2_d45_13_24,biflocation(1),...
% %     'usercond',cnd_13,'outputfuncs',true,...
% %     'step',0.1,sbargs{:})
% % %%
% % figure(1)
% % %clf; 
% % % hold on
% % eq_134=br_contn(fun_eq13,eq_134,15);
% % eq_134=br_rvers(eq_134)
% % eq_134=br_contn(fun_eq13,eq_134,15);
% % %%
% % %%
% % %%
% % 
