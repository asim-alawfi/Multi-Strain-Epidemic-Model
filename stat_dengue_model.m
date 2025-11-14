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
par([in.r,   in.e,     in.beta,        in.mu,            in.sigma,     in.gamma,  in.tau_delay])=... 
    [0.5,         1,        40,         0.02,           10,              0.4           0  ];
ndim=12;
x0=zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.r,-1],...
    'max_bound',   [in.gamma, 1.5;     in.r, 1.5],...
    'max_step',    [in.gamma, 0.01;   in.r, 0.01; 0,0.01]};
%%
% par(in.r)=par(in.beta)/(par(in.mu)+par(in.sigma));
% par(in.e)=par(in.mu)/(par(in.mu)+par(in.sigma));
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% compute and find trivial equilibri
[stat_br,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.01,parbd{:},...
    'print_residual_info',true);
%%
stat_br=br_contn(funcs,stat_br,200);
stat_br = br_rvers(stat_br);
stat_br=br_contn(funcs,stat_br,200);
%%
 [stat_br,nunst_st]=br_stabl(funcs,stat_br,0,0);
%%
r_stst_par=arrayfun(@(x)x.parameter(in.r),stat_br.point);
x1_axis=arrayfun(@(x)x.x(1),stat_br.point);
bif_point=find(diff(nunst_st));

% Extract stability information and plot results
figure(1);
clf
hold on 
plot(r_stst_par(nunst_st==0), x1_axis(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(r_stst_par(nunst_st>=1), x1_axis(nunst_st>=1), 'r.','Linewidth',2,'MarkerSize',10);
xlabel('r');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
ex1234=permfix('stst','eq1234_0',ndim,'x','perm',[1,2,3,4]);
[fun_1234,eq_1234_0,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.r,'step',0.02,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(99)
clf; hold on
eq_1234_0=br_contn(fun_1234,eq_1234_0,100);
%%
eq_1234_0=br_rvers(eq_1234_0)
eq_1234_0=br_contn(fun_1234,eq_1234_0,100);
%%
[eq_1234_0,nunst_st2]=br_stabl(fun_1234,eq_1234_0,0,0);
r_stst_par=arrayfun(@(x)x.parameter(in.r),eq_1234_0.point);
x1_axis=arrayfun(@(x)x.x(1),eq_1234_0.point);
bif_point2=find(diff(nunst_st2));
%%
check_sb=@(cnd)SetupFold(funcs,eq_1234_0,bif_point2,...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% 1=3, 2=4
ev13_24=[permfix('stst','p13',ndim,'v','perm',[1,3])]%,...
    %  permfix('stst','p24',ndim,'v','perm',[2,4])]
sv_13_24=check_sb(ev13_24)'
%1=2,3=4
ev12_34=[permfix('stst','p12',ndim,'v','perm',[1,2]),...
         permfix('stst','p34',ndim,'v','perm',[3,4])];
sv_12_34=check_sb(ev12_34)'
ev_1234=permfix('stst','p1234',ndim,'v','perm',[1,2,3,4]);
vs_1234=check_sb(ev_1234)'
