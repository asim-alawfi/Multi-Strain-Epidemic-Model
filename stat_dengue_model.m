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
    [0.1,         1,        40,           0.02,            10,         0.54         0  ];
ndim=12;
x0=zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.r,-0.5],...
    'max_bound',   [in.gamma, 1;     in.r, 1.6],...
    'max_step',    [in.gamma, 0.01;   in.r, 0.01; 0,0.01]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma))
par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% compute and find trivial equilibri
ex1234=permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4])%,...
     %  permfix('stst','eq13_24_0',ndim,'x','perm',[2,4])];
[fun_1234,eq_1234_0,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.001,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(99)
clf; hold on
eq_1234_0=br_contn(fun_1234,eq_1234_0,300);
eq_1234_0=br_rvers(eq_1234_0)
eq_1234_0=br_contn(fun_1234,eq_1234_0,100);
%%
[eq_1234_0,nunst_st]=br_stabl(fun_1234,eq_1234_0,0,0);
%% plot stability
r_stst_par=arrayfun(@(x)x.parameter(in.r),eq_1234_0.point);
x1_axis=arrayfun(@(x)x.x(1),eq_1234_0.point);
bif_point=find(diff(nunst_st));

% Extract stability information and plot results
figure(1);
clf
hold on 
plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st==1), x1_axis(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st>=2), x1_axis(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%% Invariant subspaces with respect to the group D_4:
%(1) Fully symetric invariant subspace <(1,1,1,1)> invariant under the
%action of the subgroup <r=(1234)>:={e,r,r^2,r^3}
% odd invariant subspace <(1,-1,1,-1)>  invariant under the
%action of the subgroup :={e,(13)(24),(13),(24)}, expected to detect Hopf
%bifurcation
check_sb=@(cnd)SetupFold(fun_1234,eq_1234_0,bif_point,...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% 1=3, 2=4
ev13_24=[permfix('stst','p1234',ndim,'v','perm',[1,3]),...
         permfix('stst','p34',ndim,'v','perm',[2,4])]
sv_13_24=check_sb(ev13_24)'
%1=2,3=4
ev13=permfix('stst','p12',ndim,'v','perm',[1,3]);
sv_12_34=check_sb(ev13)'
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
vs_1234=check_sb(ev_1234)'
%%
cnd_sb=[permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]),...
        permfix('stst','pv13',ndim,'v','perm',[1,3]),...
        permfix('stst','pv24',ndim,'v','perm',[2,4])];
cnd_13=permfix('stst','px13_24',ndim,'x','perm',[1,3]);
sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
[fun_eq12_34,eq_13_24,suc]=SetupStstFrom_fold(funcs,eq_1234_0,bif_point(1),...
    'usercond',cnd_13,'outputfuncs',true,...
    'step',0.1,sbargs{:});
fprintf('Branch off from symmetry breaking, keeping (13)(24) symmetry\n');
disp(eq_13_24);
figure(1);
eq_13_24=br_contn(fun_eq12_34,eq_13_24,200);
eq_13_24=br_rvers(eq_13_24);
eq_13_24=br_contn(fun_eq12_34,eq_13_24,200);
%% compute stability
[eq_13_24ref,eq_123tests,eq_123bifs,eq_123bifind]=...
    MonitorChange(fun_eq12_34,eq_13_24,'print_residual_info',0,'min_iterations',5);
%% plot stability
r_stst_par2=arrayfun(@(x)x.parameter(in.r),eq_13_24ref.point);
x1_axis2=arrayfun(@(x)x.x(1),eq_13_24ref.point);

% Extract stability information and plot results
figure(1);
hold on 
plot(r_stst_par2(eq_123tests==0), x1_axis2(eq_123tests==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests==1), x1_axis2(eq_123tests==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests==2), x1_axis2(eq_123tests==2), 'g.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par2(eq_123tests>=3), x1_axis2(eq_123tests>=3), 'b.','Linewidth',2,'MarkerSize',10);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');

%%
figure(3);clf;ax3=gca;hold(ax3,'on');%grid(ax3,'on');
args={'x',@(p)p.parameter(in.r),'y',@(p)p.x(1),'ax',ax3};
Plot2dBranch(eq_13_24ref,'DisplayName','(13)','color',[0,0,0.5],args{:});
Plot2dBranch(eq_1234_0,'DisplayName','(1234)','color',[0,0.5,0],args{:});
xlabel(ax3,'a');ylabel(ax3,'x');
%%
eq_1234_0=br_remove_extracolumns(eq_1234_0);
eq_13_24ref=br_remove_extracolumns(eq_13_24ref);
eq_13_24=br_remove_extracolumns(eq_13_24);
%
save('stat_dengue_model.mat')
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