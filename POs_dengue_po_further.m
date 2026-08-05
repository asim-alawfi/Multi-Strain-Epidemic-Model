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

%% %%%%

zro1=find(eq_1234tests==0,1,'first');
one1=find(eq_1234tests==1,1,'last');
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234=arrayfun(@(p)sum(p.x(end-3:end)),eq_1234ref.point);
%% %%%%%%%%%%%%
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
load('run30_ode45.mat')
%
npt=length(pt);
figure(1);clf;ax=gca;
s=linspace(1,3,9);
for i=1:length(s)    
    pt_in=eq_13_24ref.point(bifind_13_24(3));
    pt{npt+i}=pt_in;
    pt{npt+i}.parameter(in.r)=pt_in.parameter(in.r)+s(i); %0.5
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@(t,x)MyEvent(x));
    Fode45=@(t,X)funcs.wrap_rhs(X,pt{npt+i}.parameter);
    rng(npt+i)
    simulode45{npt+i}=ode45(Fode45,[0,1e5],simulode45{npt+i-1}.ye(:,end),opts);
    sol=simulode45{npt+i};
    t=linspace(sol.xe(end-1),sol.xe(end),10001);
    y=deval(sol,t);
    plot(ax,t-t(1),y(1:4,:));
    title(ax,sprintf('i=%d of %d, p=%g',i,length(s),pt{i}.parameter(in.r)));
    drawnow
end
%%

%%
cnd_13_24=[permfix('psol','p13',nstrains,'profile','perm',[1,3]),...
           permfix('psol','p24',nstrains,'profile','perm',[2,4]),...
           permfix('psol','p24',nstrains,'profile','perm',[1,2],'shift',[1,2]),...
           permfix('psol','p24',nstrains,'profile','perm',[3,4],'shift',[1,2])];
sim_index=2;%length(simulode45);
sol=simulode45{sim_index};
par0=pt{sim_index}.parameter;
Wm=sparse(dim,dim);
Wm(5,5)=1;
pnorm5=@(p)sqrt(dde_coll_dot(p,p,'free_par_ind',[],'period',false,'inner_matrix',Wm));
dim=size(sol.y,1);
parbd(5:end)=[];
[funcs_po,po_branch,sucin]=correct_ini_ode45sol(funcs,sol,[],par0,...
    'corpar',[1,2],'phase_condition',1,'matrix','sparse',parbd{:},'step',1e-2,...
    'degree',5,'intervals',100,'print_residual_info',true,...
    'use_tangent',true,'minimal_angle',0,'steplength_minimum',1e-5,...
    'adapt_mesh_after_correct',1,'outputfuncs',true,'extracolumns','auto',...
    'usercond',cnd_13_24,...
    'max_number_of_eigenvalues',dim,...
    'plot_measure',{@(p)p(end).period,@(p)pnorm5(p(end))},...
    'stops',{@(p)p(end).period>sol.xe(end)-sol.xe(end-1)+20});
%%
figure(1);clf;ax=gca;
xlabel(ax,'period');
ylabel(ax,'$\|x_5(\cdot)\|_2$','Interpreter','latex');
[po_branch,suc]=br_contn(funcs_po,po_branch,100,'ax',ax);
po_branch=br_rvers(po_branch);
[po_branch,suc,fail,rjct]=br_contn(funcs_po,po_branch,100,'ax',ax);
%%
[po_branch,nunst]=br_stabl(funcs_po,po_branch,'recompute',true);
%%
ini_period=sol.xe(end)-sol.xe(end-1);
[~,ix]=min(abs(arrayfun(@(p)p.period,po_branch.point)-ini_period));
sim2_i=ix;%length(simulode45);
psol=po_branch.point(sim2_i);
Fode45=@(t,X)funcs.wrap_rhs(X,psol.parameter);
rng(1)
[hdev,ntrial,dim]=deal(1e-5,10,size(psol.profile,1));
figure(1);clf;ax=gca;hold(ax,'on');
scal=diag([1/psol.parameter(in.e),1,1]);
clear y dsol
for i=2:-1:1
    dev=(rand(dim,1)-0.5)*hdev;
    dsol(i)=ode45(Fode45,[0,5e4],psol.profile(:,end)+dev,opts);
    %if length(dsol(i).xe)<length(sol.xe)
    %    continue
    %end
    plot(dsol(i).x,scal*dsol(i).y([1,5,9],:))
    ylim([0,1])
    drawnow
end
plot((psol.mesh-1)*psol.period,scal*psol.profile([1,5,9],:))
%%
balance=@(z,w)(psol.parameter(in.r)*(1-w-(1-psol.parameter(in.gamma))*z)-1);
figure(2);clf;
plot((psol.mesh-1)*psol.period,[balance(psol.profile(5,:),psol.profile(9,:));scal*psol.profile([1,5,9],:)]);
hold on;
plot(dsol(i).x,[balance(dsol(i).y(5,:),dsol(i).y(9,:));scal*dsol(i).y([1,5,9],:)]);
grid on;
%%
save(sprintf('po_branch_%02d.mat',sim_index));
%%
