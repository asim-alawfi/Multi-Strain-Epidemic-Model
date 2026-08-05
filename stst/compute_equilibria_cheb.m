%% equilibrium values (symbolic formulas + polynomial root finding)
clear
startup_coco(fullfile('..','coco_2025January28'));
format compact
format short g
parnames = {'r','e','g'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
syms(parnames{:});
par = reshape(sym(parnames),[],1);
%% Symbolic setup
[dim,nstrains]=deal(3,4);
y = sym('y', [nstrains,1]);  % infectious strains y1..y4
z = sym('z', [nstrains,1]);  % immune responses z1..z4
w = sym('w', [nstrains,1]);  % cross-immunity w1..w4
%% matrix A
A = [1 1 0 1;
    1 1 1 0;
    0 1 1 1;
    1 0 1 1];
Lambda=r*y;
AL=A*Lambda;
dydt=Lambda.*(1-w-(1-g).*z)-y; %% need to check the term gamma*w or just gamma
dzdt=e*(Lambda.*(1-z)-z);
dwdt=e*((1-w).*AL-w);
%% impose assumptions on all quantities
assume([r>1,e>0,g>0,g<1]);
for i=1:nstrains
    assume([y(i)>=0,z(i)>=0,w(i)>=0,z(i)<=1,w(i)<=1]);
end
%% for checking stability numerically
F=sco_gen(@dengue_coco);
FJ=F('x');
[g_n,e_n]=deal(0.54,0.001996); % set numerical parameter values
par0([in.r,in.e,in.g])=...
    [    1,  e_n, g_n];
par0=par0(:);
%% for general equilibrium: w and z are function of y
ws=solve(dwdt,w);
zs=solve(dzdt,z);
y_eqn=expand(subs(subs(dydt,zs),ws));
%% mesh of r values at which we compute equilibria
r_rgfine=linspace(1.001,3.5,3500);
%% general equilibria, no entry=0
y_eqn_g=simplify(y_eqn./y);
y2s=solve(y_eqn_g(1),y(2));
y_eqn2_g=numden(subs(y_eqn_g(2),y(2),y2s));
y3s=solve(y_eqn2_g,y(3));
y_eqn34_g=numden(subs(y_eqn_g(3:4),y(2:3),[y2s;y3s]));
clear y_eqn_34_fun
for i=2:-1:1
    y_eqn_34_fun{i}=matlabFunction(y_eqn34_g(i),'Vars',[y(1),y(4),r,g]);
end
y2fun=matlabFunction(y2s,'Vars',[y(1),y(4),r,g]);
y3fun=matlabFunction(y3s,'Vars',[y(1),y(4),r,g]);
y_eqn_34_funall=matlabFunction(y_eqn34_g,'Vars',[y(1),y(4),r,g]);
y_eqn_funall=matlabFunction(y_eqn_g,'Vars',{y(:),[r;g]});
%%
addpath(fullfile('..','chebfun'));
%%
r_rg=linspace(1.01,3.5,350);
dom=[0,1,0,1];
info_g=struct('y_eqn_fun',{y_eqn_34_fun},'dom',dom,'y2fun',y2fun,'y3fun',y3fun);
clear y_g
parfor i=1:length(r_rgfine)
    y_g{i}=get_eq_cheb2(r_rg(i),g_n,info_g,true);
    fprintf('i=%d of %d\n',i,length(r_rgfine));
end
%%
for i=1:length(r_rg)
    y_g{i}=y_g{i}(:,~any(isnan(y_g{i})));
end
%%
save('stst_sym_cheb.mat')