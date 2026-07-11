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
pt.parameter(in.r)=pt.parameter(in.r)+0.3;
%pt.parameter(in.gamma)=0.4;
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@(t,x)MyEvent(x));
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,1e5],pt.x+1e-5*rand(length(pt.x),1),opts);
%%
figure(88)
clf;
%tiledlayout(1,4)
%nexttile
hold on
plot(simulode45.x,simulode45.y(5,:),'ko','LineWidth',1)
plot(simulode45.x,simulode45.y(7,:),'bx','LineWidth',1)

plot(simulode45.x,simulode45.y(6,:),'r.','MarkerSize',10,'LineWidth',1)
plot(simulode45.x,simulode45.y(8,:),'s','LineWidth',1)
%
figure(90)
clf;
%tiledlayout(1,4)
%nexttile
hold on
plot(simulode45.x,simulode45.y(1,:),'ko','LineWidth',2)
plot(simulode45.x,simulode45.y(3,:),'bx','LineWidth',2)

plot(simulode45.x,simulode45.y(2,:),'r.','MarkerSize',10,'LineWidth',2)
plot(simulode45.x,simulode45.y(4,:),'s','LineWidth',2)

%%
save('run_ode45.mat')
%%
Rsym=zeros(12);
R1=Rsym(1:2,:);
R1(1,7)=1;
R1(2,5)=1;
R2=Rsym(1:2,:);
R2(1,8)=1;
R2(2,6)=1;
%disp(R1)
disp(R2)
%%
% %% picking up a periodic solution from dde23 and continuing computation in PR (here df=0.73).
% [funcs_s2,po2_symmetry,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df,in.PR,par1,...
%     'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
%     'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
%     'adapt_mesh_after_correct',1,...
%     'usercond',{psolsym},'outputfuncs',true);

ev1=zeros(2,ndim);
ev1(1,5)=1;
ev1(2,7)=1;
ev2=zeros(2,ndim);
ev2(1,6)=1;
ev2(2,8)=1;
%%
psolsym=@(p)dde_psol_lincond(p,'profile','trafo',R1,'shift',[1,2],...
    'condprojint',linspace(0.1,0.25,2)'*[1,1],'condprojmat',ev1);  %'condprojmat',[1,0,0,0,0,0]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);

%%


cnd_13=permfix('psol','p13',nstrains,'x','perm',[1,3],'shift',[0,1])
cnd_24=permfix('psol','p24',nstrains,'x','perm',[2,4],'shift',[0,1]);
impos_cond={'usercond',[cnd_13,cnd_24],'initcond',[cnd_13,cnd_24],...
           'extracolumns','auto','extra_condition',true,};
%%
mycnd_13=Mypermfix('psol','p13',nstrains,'x','trafomat',R1,'stateproj',ev1,'shift',[0,1]);
mycnd_24=Mypermfix('psol','p13',nstrains,'x','trafomat',R2,'stateproj',ev1,'shift',[0,1]);
impos_cond2={'usercond',[mycnd_13,mycnd_24],'extracolumns','auto','extra_condition',true,};
%%
[funcs_psol,psol_12_34,suc]=branch_from_ode45sol(funcs,simulode45,in.r,pt.parameter,...
    'phase_condition',1,'stability.eigmatrix','sparse',parbd{:},'step',1e-5,...
    'newton_max_iterations',8, 'adapt_mesh_after_correct',1, ...
    'matrix','sparse','eigmatrix','sparse','degree',6,'intervals',150,'max_step',[0,1e-5],...
    'print_residual_info',1,'outputfuncs',true)%,...
  %  'plot_measure',{@(p)p.parameter(in.r),@(p)max(mean(p.profile,1))})%,impos_cond2{:});
%%
figure(17)
clf
hold on
psol_12_34=br_contn(funcs_psol,psol_12_34,3);
%%


%%
%po=br_stabl(funcs_psol,psol_12_34);
%
for i=1:length(psol_12_34.point)
po=psol_12_34.point(i);
figure(11)
clf;
hold on
plot(po.mesh*po.period,po.profile(5:8,:),'LineWidth',2)
end
%%
[psol_12_34_wbifs,nunst,bifs,bifind]=MonitorChange(funcs_psol,psol_12_34,...
    'range',2:length(psol_12_34.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
%%
[psol_12_34,unst]=br_stabl(funcs_psol,psol_12_34,0,0)
%% [nunst2_sym,dom,triv2_defect_sym,po2_symmetry.point]=GetStability(po2_symmetry,'funcs',funcs_s2,...
%     'exclude_trivial',true);%,'recompute',true);
chang2_stb_sym=find(diff(nunst));
rp2_x=arrayfun(@(x)x.parameter(in.r),psol_12_34_wbifs.point);
max2_y=arrayfun(@(x)max(x.profile(1,:)),psol_12_34_wbifs.point);
min2_y2=arrayfun(@(x)min(x.profile(1,:)),psol_12_34_wbifs.point);
% plotting stability
figure(2);clf;hold on;grid on
plot(rp2_x(nunst==0),max2_y(nunst==0),'bo',rp2_x(nunst>=1),max2_y(nunst>=1),'kx')%,'LineWidth',1)
plot(rp2_x(nunst==0),min2_y2(nunst==0),'bo',rp2_x(nunst>=1),min2_y2(nunst>=1),'kx')%,
%%
ptf=eq_13_24ref.point(bifind_13_24(3)+1)