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



%% Trying to compute the branch with respect to (13) and y4=0
y4eq0=dde_sys_cond_create('name','y4_eq0','reference',false,'fun',@(x)x,'args',{'x',4});
y4veq0=dde_sys_cond_create('name','y4_v0','reference',false,'fun',@(x)x,'args',{'v',4});
%alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
cndv_13=[permfix('stst','pv13',nstrains,'v','perm',[1,3])];
cnd_13=[permfix('stst','p13',nstrains,'x','perm',[1,3])];
sbargs_13_br2={'SetupFold.usercond',[y4eq0,cnd_13,cndv_13],'SetupFold.initcond',[y4veq0,y4eq0,cndv_13]};
[fun_eq13_br2,eq_13_br2,suc_13br2]=SetupStstFrom_fold(funcs,eq_13_240ref,bifind_13_240(3),...
    'usercond',[cnd_13,y4eq0],'initcond',[cnd_13],'outputfuncs',true,...
    'step',0.1,sbargs_13_br2{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)},...
    'print_residual_info',true,parbd{:});
fprintf('Branch off from symmetry breaking, ~ (13)(24) symmetry\n');
eq_13_br2.point(2).x(1:4,:)
%%
figure(1)
hold on
eq_13_br2=br_contn(fun_eq13_br2,eq_13_br2,20);

