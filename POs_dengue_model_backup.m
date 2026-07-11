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

 %%%%%%

zro1=find(eq_1234tests==0,1,'first');
one1=find(eq_1234tests==1,1,'last');
xpar_1234=arrayfun(@(x)x.parameter(in.r),eq_1234ref.point);
ypar_1234=arrayfun(@(p)sum(p.x(end-3:end)),eq_1234ref.point);
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
s=linspace(0,1,12);
simulode45 = cell(length(s),1);
pt=cell(length(s),1);
for i=1:length(s)    
pt_in=eq_13_24ref.point(bifind_13_24(3));
pt{i}=pt_in;
pt{i}.parameter(in.r)=pt_in.parameter(in.r)+s(i); %0.5
opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@(t,x)MyEvent(x));
Fode45=@(t,X)funcs.wrap_rhs(X,pt{i}.parameter);
rng(i)
simulode45{i}=ode45(Fode45,[0,1e5],pt{i}.x+1e-4*rand(length(pt{i}.x),1),opts);
end
%%
save('run30_ode45.mat')
%% Plotting the simulation from ode45 
% here we plot z fully immune population ( that is y(5:8,:))
% 12 simulations have been computed and plotted into 2 Figures 
numPerFigure = 6;
numFigures   = 2;
for figIdx = 1:numFigures
    figure(30 + figIdx)
    clf
    tiledlayout(3,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    for panelIdx = 1:numPerFigure
        % Simulation number: 1,...,30
        simIdx = (figIdx-1)*numPerFigure + panelIdx;
            sol = simulode45{simIdx};
            nexttile
            hold on
            plot(sol.x, sol.y(5,:), 'ko', ...
                'LineWidth',2)
            plot(sol.x, sol.y(7,:), 'bx', ...
                'LineWidth',2)
            plot(sol.x, sol.y(6,:), 'r.', ...
                'MarkerSize',10, ...
                'LineWidth',2)
            plot(sol.x, sol.y(8,:), 's', ...
                'LineWidth',2)
            title(sprintf('Simulation %d, s = %.3g', ...
                simIdx, s(simIdx)))
            xlabel('Time')
            ylabel('Solution')
            grid on
            box on
        if panelIdx == 1
            legend('y_5','y_7','y_6','y_8', ...
                'Location','best')
        end
    end
    sgtitle(sprintf('Simulations %d–%d', ...
        (figIdx-1)*numPerFigure + 1, ...
        figIdx*numPerFigure))
end

%%

%%
alleq0=dde_sys_cond_create('name','all_eq0','reference',false,'fun',@(x)x,'args',{'x',1:ndim});
% cndv_13_24=[permfix('stst','p13',nstrains,'v','perm',[1,3]),..clear.
%            permfix('stst','p24',nstrains,'v','perm',[2,4])];
cnd_13_24=[permfix('psol','p13',nstrains,'profile','perm',[1,3]),...
           permfix('psol','p24',nstrains,'profile','perm',[2,4])];
sim_index=6;
[funcs_po,po_branch,sucin]=branch_from_ode45sol(funcs,simulode45{sim_index},in.r,pt{sim_index}.parameter,...
    'phase_condition',1,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},'step',1e-7,...
    'degree',6,'intervals',150,'newton_max_iterations',8,'print_residual_info',true,...
    'adapt_mesh_after_correct',1,'outputfuncs',true);
figure(22)
clf 
hold on
[po_branch,suc,fail,rjct]=br_contn(funcs_po,po_branch,2);
%%
po=br_stabl(funcs_po,po_branch)





%%
% figure(88)
% clf;
% %tiledlayout(1,4)
% %nexttile
% hold on
% plot(simulode45.x,simulode45.y(5,:),'ko','LineWidth',2)
% plot(simulode45.x,simulode45.y(7,:),'bx','LineWidth',2)
% 
% plot(simulode45.x,simulode45.y(6,:),'r.','MarkerSize',10,'LineWidth',2)
% plot(simulode45.x,simulode45.y(8,:),'s','LineWidth',2)
% ylabel('Z','LineWidth',15)
% %%
% figure(90)
% clf;
% %tiledlayout(1,4)
% %nexttile
% hold on
% plot(simulode45.x,simulode45.y(1,:),'ko','LineWidth',2)
% plot(simulode45.x,simulode45.y(3,:),'bx','LineWidth',2)
% 
% plot(simulode45.x,simulode45.y(2,:),'r.','MarkerSize',10,'LineWidth',2)
% plot(simulode45.x,simulode45.y(4,:),'s','LineWidth',2)
% ylabel('y','LineWidth',15)































































% Rsym=zeros(12);
% Rsym(1,3)=1;
% Rsym(3,1)=1;
% disp(Rsym)
% e_v=zeros(1,ndim);
% e_v(1)=1;
% e_v(3)=1;
% % xdim=length(x0);
% % % symmetric conditios
% psolsym=@(p,pref)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,4],...
%     'condprojint',linspace(0.1,0.25,2)'*[1,1],'condprojmat',e_v);  %'condprojmat',[1,0,0,0,0,0]);
% addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
%     'UniformOutput',false),args(2:2:end)),1,[]);
% % %% picking up a periodic solution from dde23 and continuing computation in PR (here df=0.73).
% % [funcs_s2,po2_symmetry,suc2]=branch_from_sol_audi_sym(funcs_audi,sol23_df,in.PR,par1,...
% %     'extra_condition',true,'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
% %     'degree',6,'intervals',80,'extracolumns','auto','newton_max_iterations',8,...
% %     'adapt_mesh_after_correct',1,...
% %     'usercond',{psolsym},'outputfuncs',true);