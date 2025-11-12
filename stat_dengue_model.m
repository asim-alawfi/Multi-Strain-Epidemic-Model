clear;
base=[pwd(),filesep(),'..',filesep(),filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');%[getenv('HOME'),'/sourceforge-ddebiftool/releases/git-2023-06-27'],'sym');

format compact
format short g
%%

parnames = {'r',  'e',  'beta',   'mu',   'sigma',   'gamma',  'tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});

par([in.r,   in.e,     in.beta,      in.mu,            in.sigma,     in.gamma,  in.tau_delay])=... 
    [0.5,         1,        40,         0.02,           10,              0.4           0  ];

x0=zeros(12,1);
%x0=[0;0;0;0;...
    %0;0;0;0;...
    %0;0;0;0];

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
