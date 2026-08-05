% Extract stability information and plot results
clear;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(fullfile('..','ddebiftool_snapshot'));
format compact
format short g
stst=load(fullfile('..','stst','stst_sym_e_0_019.mat'));
[nstrains,dim]=deal(stst.nstrains,stst.dim);
%% Define parameters
parnames = {'r','e','g'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
parbd={'min_bound',[in.g,  0;  in.r,0.0],...
    'max_bound',   [in.g,  1;  in.r,3.5],...
    'max_step',    [   0,100;  in.r,0.1]};
funcs=set_symfuncs(@dengue_fever_sym,'sys_tau',@()[]);
%% find Hopf bifurcation for (13) equilibrium
eq13=df_brnch(funcs,in.r,'stst');
sel=~isnan(stst.x13(1,:,2));
eq13.point=cellfun(@(x,p)dde_stst_create('x',x,'parameter',reshape(p,1,[])),...
    num2cell(stst.x13(:,sel,2),1),num2cell(stst.par13(:,sel),1));
[eq13,nunsteq13]=br_stabl(funcs,eq13,'max_number_of_eigenvalues',nstrains*dim);
%%
indhopf=find(abs(diff(nunsteq13))==2,1,'first');
[psol13,suc]=SetupPsol(funcs,eq13,indhopf,'degree',5,'intervals',40,'print_residual_info',1,...
    'plot_measure',{@(p)p(end).parameter(in.r),@(p)max(p.profile(1,:))},...
    'stops',{@(p)p(end).period*p(end).parameter(in.e)>100});
%%
figure(2);clf;ax2=gca;
xlabel(ax2,'r');
ylabel(ax2,'max(y_1)')
psol13=br_contn(funcs,psol13,100,'ax',ax2);
[psol13,pnunst13]=br_stabl(funcs,psol13,'exclude_trivial',true,'max_number_of_eigenvalues',nstrains*dim,'recompute',true);
%% plot example orbit
figure(3);clf;
ptplot=psol13.point(end-5);plot(ptplot.mesh*ptplot.period,ptplot.profile([1,2,5,6,9,10],:))
%% run simulation at parameter value where no stable equilibrium exists
pt=eq13.point(indhopf+2);
opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@(t,x)event_cross_poincare(x,1,0.2,1));
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1);
h=1e-2;
dev=2*(rand(size(pt.x))-0.5);
sol=ode45(Fode45,[0,1e5],pt.x+h*dev,opts);
t=linspace(sol.xe(end-1),sol.xe(end),10001);
y=deval(sol,t);
figure(3);clf;pax=gca;
plot(pax,t-t(1),y(1:4,:));
title(pax,sprintf('r=%g',pt.parameter(in.r)));
drawnow
%%
%cnd_13_24=[permfix('psol','p13',nstrains,'profile','perm',[1,3]),...
%    permfix('psol','p24',nstrains,'profile','perm',[2,4]),...
%    permfix('psol','p24',nstrains,'profile','perm',[1,2],'shift',[1,2]),...
%    permfix('psol','p24',nstrains,'profile','perm',[3,4],'shift',[1,2])];
par0=pt.parameter;
alldim=size(sol.y,1);
[po13_24,suc]=correct_ini_ode45sol(funcs,sol,in.r,par0,...
    'corpar',in.r,'matrix','sparse',parbd{:},'step',1e-2,...
    'degree',5,'intervals',100,'print_residual_info',true,...
    'use_tangent',true,'minimal_angle',0,'steplength_minimum',1e-5,...
    ...'extracolumns','auto','usercond',cnd_13_24,...
    'max_number_of_eigenvalues',alldim,...
    'plot_measure',{@(p)p.parameter(in.r),@(p)max(p.profile(1,:))},...
    'stops',{@(p)p(end).period*p(end).parameter(in.e)>100});
%%
figure(2);clf;pax=gca;
xlabel(pax,'$r$','Interpreter','latex');
ylabel(pax,'$\max y_1$','Interpreter','latex');
[po13_24,suc]=br_contn(funcs,po13_24,100,'ax',pax);
po13_24=br_rvers(po13_24);
po13_24=br_contn(funcs,po13_24,100,'ax',pax);
%%
[po13_24,pnunst13_24]=br_stabl(funcs,po13_24,'recompute',true);
%%
save('po13_24_results.mat')
%% plot example orbit
figure(3);clf;tl=tiledlayout(1,2,'TileSpacing','tight');
pax(1)=nexttile(tl,1);pax(2)=nexttile(tl,2);
rpo13_24=arrayfun(@(p)p.parameter(in.r),po13_24.point);
e_unit=@(k)double(1:alldim==k);
p_amp=@(p)max(e_unit(1)*dde_coll_eva(p,dde_coll_roots(p,e_unit(1),'diff',1)));
amp_po13_24=arrayfun(p_amp,po13_24.point);
plot(pax(2),rpo13_24,amp_po13_24,'.-');
hold(pax(2),'on');
pl=plot(pax(2),rpo13_24(end),amp_po13_24(end),'ko','MarkerFaceColor','k','MarkerSize',8);
for i=length(po13_24.point):-1:1
    ptplot=po13_24.point(i);
    hold(pax(1),'off')
    plot(pax(1),ptplot.mesh*ptplot.period,ptplot.profile([1,5,9],:),'.-');
    hold(pax(1),'on');
    plot(pax(1),ptplot.mesh*ptplot.period,ptplot.profile([2,6,10],:),'--','LineWidth',2);
    axis(pax(1),'tight');
    title(pax(1),sprintf('r=%g',ptplot.parameter(in.r)))
    pl.XData=ptplot.parameter(in.r);
    pl.YData=p_amp(ptplot);
    pause(0.1)
    drawnow
end
%%
figure(1);
ax=findall(get(gcf,'Children'),'Type','axes');
plot(ax(2),rpo13_24,amp_po13_24,'ko-','MarkerFaceColor','k','MarkerSize',5,...
    'DisplayName','po(13)(24)_{1/2}');
set(ax(2),'FontSize',12);
set(ax(1),'FontSize',12);
figpos=get(gcf,'Position');
set(gcf,'Position',[figpos(1:2),1064,692]);
%%
exportgraphics(figure(1),'bifdiagram_r.pdf','ContentType','vector')
exportgraphics(figure(1),'bifdiagram_r.html')