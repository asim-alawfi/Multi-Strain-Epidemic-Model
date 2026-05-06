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
pt.parameter(in.r)=pt.parameter(in.r)+0.2;
opts = odeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,x)MyEvent(x));
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,1e5],pt.x+0.001*rand(length(pt.x),1),opts);
%%
figure(88)
clf;
%tiledlayout(1,4)
%nexttile
hold on
plot(simulode45.x,simulode45.y(5:8,:),'LineWidth',2)


%%
save('run_ode45.mat')
%%
Rsym=zeros(12);
Rsym(1,3)=1;
Rsym(3,1)=1;
disp(Rsym)
e_v=zeros(1,ndim);
e_v(1)=1;
e_v(3)=1;
% xdim=length(x0);
% % symmetric conditios
psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,4],...
    'condprojint',linspace(0.1,0.25,2)'*[1,1],'condprojmat',e_v);  %'condprojmat',[1,0,0,0,0,0]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
% %% picking up a periodic solution from dde23 and continuing computation in PR (here df=0.73).
% [funcs_s2,po2_symmetry,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df,in.PR,par1,...
%     'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
%     'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
%     'adapt_mesh_after_correct',1,...
%     'usercond',{psolsym},'outputfuncs',true);
%%
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
% cndv_13_24=[permfix('stst','p13',nstrains,'v','perm',[1,3]),...
%            permfix('stst','p24',nstrains,'v','perm',[2,4])];
cnd_13_24=[permfix('psol','p13',nstrains,'profile','perm',[1,3]),...
           permfix('psol','p24',nstrains,'profile','perm',[2,4])];
[funcs_psol,psol_12_34,suc]=branch_from_ode45sol(fun_eq13_24,simulode45,in.r,par,...
    'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',6,'intervals',80,'newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',{psolsym},'outputfuncs',true);
%%
figure(17)
clf
hold on
psol_12_34=br_contn(funcs_psol,psol_12_34,5);