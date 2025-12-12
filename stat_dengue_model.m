clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
%   base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% addpath([base,filesep,'ddebiftool'],...
%     [base,filesep,'ddebiftool_extra_psol'],...
%     [base,filesep,'ddebiftool_utilities'],...
%     [base,filesep,'ddebiftool_extra_rotsym'],...
%     [base,filesep,'ddebiftool_extra_nmfm'],...
%     [base,filesep,'ddebiftool_extra_symbolic'],...
%     [base,filesep,'ddebiftool_coco']);
format compact
format short 
%%
parnames = {'r',  'e',  'beta',   'mu',   'sigma',   'gamma',  'tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.r,      in.e,     in.beta,        in.mu,      in.sigma,     in.gamma,  in.tau_delay])=... 
    [0.1,         1,        40,           0.02,            10,         0.54         0  ];
ndim=12;
x0=zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.r,-0.1],...
    'max_bound',   [in.gamma, 1;     in.r, 1.6],...
    'max_step',    [in.gamma, 0.01;   in.r, 0.01; 0,0.01]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% Tracking a branch of equilirium points
% We impose symmetry with respect to the permutation (1234) to avoid
% singularity along the branch
ex1234=permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4])
[fun_1234,eq_1234_0,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.001,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(1)
clf; hold on
eq_1234_0=br_contn(fun_1234,eq_1234_0,300);
eq_1234_0=br_rvers(eq_1234_0)
eq_1234_0=br_contn(fun_1234,eq_1234_0,100);
%%
[eq_1234_0,nunst_st]=br_stabl(fun_1234,eq_1234_0,0,0);
%% plot stability
r_stst_par=arrayfun(@(x)x.parameter(in.r),eq_1234_0.point);
x1_axis=arrayfun(@(x)x.x(1),eq_1234_0.point);
bif_point=find((nunst_st>0),1,'last');
p1=eq_1234_0.point(bif_point);
%%
% Extract stability information and plot results
figure(1);
clf
hold on 
plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st==1), x1_axis(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st>=2), x1_axis(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(bif_point), x1_axis(bif_point), 'r.','Linewidth',2,'MarkerSize',40);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
eq_1234_0=br_remove_extracolumns(eq_1234_0);
save('computation_without_delay.mat');
%%
cnd_sb=[permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]),...
        permfix('stst','pv13',ndim,'v','perm',[1,2,3,4])];
cnd_13=permfix('stst','px13_24',ndim,'x','perm',[1,2,3,4]);
sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
[fun_eq12_34,eqbranch_13_24,suc]=SetupStstFrom_fold(funcs,eq_1234_0,bif_point(1),...
    'usercond',[],'outputfuncs',true,...
    'step',0.1,sbargs{:});
%fprintf('Branch off from symmetry breaking, keeping (13) symmetry\n');
%%
eqbranch_13_24.parameter.max_bound(end)=7;
disp(eqbranch_13_24);
figure(1);
eqbranch_13_24=br_contn(fun_eq12_34,eqbranch_13_24,30);
eqbranch_13_24=br_rvers(eqbranch_13_24);
eqbranch_13_24=br_contn(fun_eq12_34,eqbranch_13_24,200);
%% compute stability
[eq_13_24ref,eq_123tests,eq_123bifs,eq_123bifind]=...
    MonitorChange(fun_eq12_34,eqbranch_13_24,'print_residual_info',0,'min_iterations',5);
%% plot stability
r_stst_par2=arrayfun(@(x)x.parameter(in.r),eq_13_24ref.point);
x1_axis2=arrayfun(@(x)x.x(1),eq_13_24ref.point);

%% Extract stability information and plot results
p2=eq_13_24ref.point(eq_123bifind(1));
p3=eq_13_24ref.point(eq_123bifind(2));
p4=eq_13_24ref.point(eq_123bifind(3));
figure(11);
clf
hold on 
plot(r_stst_par2(eq_123tests==0), x1_axis2(eq_123tests==0), 'kx','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests==1), x1_axis2(eq_123tests==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests==2), x1_axis2(eq_123tests==2), 'g.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests>=3), x1_axis2(eq_123tests>=3), 'b.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st==1), x1_axis(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st>=2), x1_axis(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(bif_point), x1_axis(bif_point), 'r.','Linewidth',2,'MarkerSize',40);
%plot(r_stst_par2(eq_123bifind(1)),x1_axis2(eq_123bifind(1)),'k.','Linewidth',2,'MarkerSize',20);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
save('computation_without_delay2.mat');
%%


%% Invariant subspaces with respect to the group D_4:
%(1) Fully symetric invariant subspace <(1,1,1,1)> invariant under the
%action of the subgroup <r=(1234)>:={e,r,r^2,r^3}
% odd invariant subspace <(1,-1,1,-1)>  invariant under the
%action of the subgroup :={e,(13)(24),(13),(24)}, expected to detect Hopf
%bifurcation
check_sb=@(cnd)SetupFold(fun_1234,eq_1234_0,bif_point,...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% 1=3, 2=4
ev13_24=[permfix('stst','p13',ndim,'v','perm',[1,3]),...
         permfix('stst','p24',ndim,'v','perm',[2,4])]
sv_13_24=check_sb(ev13_24)'
%1=2,3=4
ev13=permfix('stst','p13',ndim,'v','perm',[1,3]);
sv_12_34=check_sb(ev13)'
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
vs_1234=check_sb(ev_1234)'
%%
ev12_34=[permfix('stst','p12',ndim,'v','perm',[1,2]),...
         permfix('stst','p34',ndim,'v','perm',[3,4])];
vs12_34=check_sb(ev12_34)'
%% Continue secondary symmetry breaking in 2 parameters
% cnd_sb0=[permfix('stst','px13',ndim,'x','perm',[1,2,3,4]);...
%     permfix('stst','pv13',ndim,'v','perm',[1,3]); ...
%     permfix('stst','pv13',ndim,'v','perm',[2,4])];
% [fun_sb2,eq_sb2,suc]=SetupFold(funcs,eq_1234_0,bif_point(1),...
%     'initcond',cnd_sb0,'usercond',cnd_sb0,'outputfuncs',true,...
%     'contpar',[in.r,in.gamma],'dir',in.r,'step',0.1,'print_residual_info',1);
% fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
% disp(eq_sb2)
% %%
% figure(2);
% clf; hold on
% eq_sb2=br_contn(fun_sb2,eq_sb2,100);

%%

%%

% 
% %%
% figure(3);clf;ax3=gca;hold(ax3,'on');%grid(ax3,'on');
% args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax3};
% Plot2dBranch(eq_13_24ref,'DisplayName','(13)','color',[0,0,0.4],args{:});
% Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.4,0],args{:});
% xlabel(ax3,'r');ylabel(ax3,'y_1');
% %%
% eq_1234_0=br_remove_extracolumns(eq_1234_0);
% eq_13_24ref=br_remove_extracolumns(eq_13_24ref);
% eqbranch_13_24=br_remove_extracolumns(eqbranch_13_24);
% save('stat_dengue_model.mat');
% %%
% %% Continue secondary symmetry breaking in 2 parameters
% cnd_sb2=[permfix('stst','px13',ndim,'x','perm',[1,2,3,4]);...
%     permfix('stst','pv13',ndim,'v','perm',[1,3]); ...
%     permfix('stst','pv13',ndim,'v','perm',[2,4])];
% [fun_sb2,eq_sb2,suc]=SetupFold(funcs,eq_13_24ref,eq_123bifind(3),...
%     'initcond',cnd_sb2,'usercond',cnd_sb2,'outputfuncs',true,...
%     'contpar',[in.r,in.gamma],'dir',in.r,'step',0.1,'print_residual_info',1);
% fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
% disp(eq_sb2)
% %%
% figure(2);
% %clf; hold on
% eq_sb2=br_contn(fun_sb2,eq_sb2,100);
% %%
% eq_sb2=br_rvers(eq_sb2);
% eq_sb2=br_contn(fun_sb2,eq_sb2,100);
% %%
% [eq_sb2ref,eq_sb2tests,eq_sb2bifs,eq_sb2bifind]=MonitorChange(fun_sb2,eq_sb2,...
%     'printlevel',1,'print_residual_info',0,'min_iterations',5,'output','branch');
% %% Continue symmetry breaking in 2 parameters
% cnd_sb=[permfix('stst','px1234',ndim,'x','perm',[1,3]);...
%     permfix('stst','pv123',ndim,'v','perm',[1,3])];
% [fun_sb,eq_sb,suc]=SetupFold(funcs,eq_13_24ref,eq_123bifind(1)+1,...
%     'initcond',cnd_sb,'usercond',cnd_sb,'outputfuncs',true,...
%     'contpar',[in.r,in.beta],'dir',in.r,'step',0.1,'print_residual_info',1);
% fprintf('Initialize two-parameter continuation of  symmetry breaking\n');
% disp(eq_sb)
% figure(2);clf;ax2=gca;
% xlabel(ax2,'a');ylabel(ax2,'c_{ext}');
% eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
% eq_sb=br_rvers(eq_sb);
% eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
% [eq_sbref,eq_sbtests,eq_sbifs,eq_sbbifind]=MonitorChange(fun_sb,eq_sb,...
%     'printlevel',1,'print_residual_info',0,'min_iterations',5,'output','branch',...
%     'locate_trivial',@(p)[0,0,0]);





%% simulation with dde23
pt=eq_1234_0.point(bif_point-1);
%
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter)
rng(1)
simulode4523=ode45(Fode45,[0,20000],pt.x+0.01*rand(length(pt.x),1))
%%
figure(88)
clf;
hold on
plot(simulode4523.x,simulode4523.y(9:12,:),'LineWidth',2)























% cnd_sb2=[permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]),...
%         permfix('stst','pv13',ndim,'v','perm',[1,2,3,4])];
% cnd_1234=permfix('stst','px1234',ndim,'x','perm',[1,2,3,4]);
% 
% %%
% sbargs2={'SetupFold.usercond',cnd_sb2,'SetupFold.initcond',cnd_sb2};
% [fun_eq1234_bk,eq1234_bk,suc]=SetupStstFrom_fold(funcs,eq_1234_0,bif_point(1),...
%     'usercond',cnd_1234,'outputfuncs',true,...
%     'step',0.1,sbargs2{:});
% fprintf('Branch off from symmetry breaking, keeping (1234) symmetry\n');
% disp(eq1234_bk);
% figure(1);
% eq1234_bk=br_contn(fun_eq1234_bk,eq1234_bk,200);
% eq1234_bk=br_rvers(eq1234_bk);
% eq1234_bk=br_contn(fun_eq1234_bk,eq1234_bk,100);
%%
% [branch,indbif,indmap,notcorrected]=br_bisection(funcs,branch,inds,'fold')