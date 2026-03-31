% Extract stability information and plot results
clear;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(fullfile('ddebiftool_snapshot'));
format compact
format short g
load('steady_state_y2_0.mat')
%%
clrs=orderedcolors("gem12");
figure(10);
clf
hold on
xpar_12_34=arrayfun(@(x)x.parameter(in.r),eq_12_34ref.point);
ypar_12_34=arrayfun(@(p)sum(p.x(end-3:end)),eq_12_34ref.point);
%figure(33)

%plt1_12_34=plot(xpar_12_34(eq_12_34tests==0), ypar_12_34(eq_12_34tests==0),'-','Color','red','Linewidth',3,'MarkerSize',10);
%plt2_12_34=plot(xpar_12_34(eq_12_34tests>=1), ypar_12_34(eq_12_34tests>=1), '--','Color',0.5*[1 1 1],'Linewidth',3,'MarkerSize',10);
%plot(xpar_12_34(bifind_12_34), ypar_12_34(bifind_12_34), 'k.','Linewidth',3,'MarkerSize',30);
% %%%%%%

zro1=find(eq_1234tests==0,1,'first');
one1=find(eq_1234tests==1,1,'last');
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234=arrayfun(@(p)sum(p.x(end-3:end)),eq_1234ref.point);
%plt3_1234=plot(xpar_1234(1:zro1), ypar_1234(1:zro1), '--','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
%plt1_1234=plot(xpar_1234(zro1:one1), ypar_1234(zro1:one1), ':','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
%plt2_1234=plot(xpar_1234(one1:end), ypar_1234(one1:end), '--','Color',clrs(3,:),'Linewidth',3,'MarkerSize',10);
%plot(xpar_1234(bifind_1234(3)), ypar_1234(bifind_1234(3)), '.','Color',clrs(3,:),'Linewidth',3,'MarkerSize',30);
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
legend('k~(13)(24) stable','k~(13)(24) unstable',... 
          '','double-zero', 'stable','unstable','bif with 4 zero eiganvalues',...
          'FontSize',17,'FontWeight','normal','location','northwest')
xlabel('$r$','interpreter','latex','FontSize',17,FontWeight='bold')
ylabel('$y_{1}$','interpreter','latex',FontSize=17,FontWeight='bold')
%title('(a)');
%%
pt=eq_13_24ref.point(bifind_13_24(3));
pt.parameter(in.r)=pt.parameter(in.r)+0.2 
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,1e5],pt.x+0.001*rand(length(pt.x),1),opts);
%%
figure(88)
clf;
tiledlayout(1,4)
nexttile
hold on
plot(simulode45.x,simulode45.y(5,:),'LineWidth',2)
nexttile 
hold on
plot(simulode45.x,simulode45.y(6,:),'LineWidth',2)
nexttile
hold on
plot(simulode45.x,simulode45.y(7,:),'LineWidth',2)
nexttile
hold on
plot(simulode45.x,simulode45.y(8,:),'LineWidth',2)
%%
y3veq0=dde_sys_cond_create('name','y3_eq0','reference',false,'fun',@(x)x,'args',{'v',3});
y1veq0=dde_sys_cond_create('name','y1_eq0','reference',false,'fun',@(x)x,'args',{'v',1});
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
cndv_13=[permfix('stst','p13',nstrains,'v','perm',[1,3])]%,...
           %permfix('stst','p24',nstrains,'v','perm',[2,4])];
cnd_13=[permfix('stst','p13',nstrains,'x','perm',[1,3])]%,...
          %permfix('stst','p24',nstrains,'x','perm',[2,4])];
y3eq0=dde_sys_cond_create('name','y3_eq0','reference',false,'fun',@(x)x,'args',{'x',3});
y1eq0=dde_sys_cond_create('name','y1_eq0','reference',false,'fun',@(x)x,'args',{'x',1});

sbargs={'SetupFold.usercond',[cndv_13,y3veq0],'SetupFold.initcond',[cndv_13,y3veq0]};
[fun_eq24,eq_24,suc_24]=SetupStstFrom_fold(funcs,eq_13_24ref,bifind_13_24(3),...
    'usercond',[cnd_13,y3eq0],'outputfuncs',true,...
    'step',0.1,sbargs{:},'plot_measure',{@(p)p(end).parameter(in.r),@(p)p(end).x(end)},'print_residual',true)
fprintf('Branch off from symmetry breaking, keeping (13) symmetry\n');
%%
check_sb=@(cnd)SetupFold(funcs,eq_13_24ref,bifind_13_24(3),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
sv_13_y0=check_sb([cndv_13,y3veq0])'
%%
% %
% figure(1)
% hold on
% eq_13_24=br_contn(fun_eq13_24,eq_13,100);
%%
