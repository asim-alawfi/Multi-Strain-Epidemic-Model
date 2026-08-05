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
npar=length(parnames);
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
%% trivial equilibria
r00=linspace(0.5,3.5,4001);
y00=zeros(nstrains,length(r00));
[x00,par00,stab00]=x_eq_from_y(r00,y00,A,par0,in,FJ);
%% y1=y3,y2=y4=0 equilibria and their bifurcation with non-zero y1
y_eqn_13_240=numden(simplify(subs(simplify(y_eqn(1:2)./y(1:2)),y(2:4),[0;y(1);0])));
ys_13_240=solve(y_eqn_13_240(1),y(1)); % find non-zero solution y
yfun_13_240=matlabFunction(ys_13_240,'Vars',[r,g]);
y_13_240=[1;0;1;0]*yfun_13_240(r_rgfine,g_n);
[x13_240,par13_240,stab13_240]=x_eq_from_y(r_rgfine,y_13_240,A,par0,in,FJ);
fprintf('Equilibrium with y1=y3,y2=y4=0: y1=%s\n',char(ys_13_240));
Js_13_24=simplify(subs(jacobian(y_eqn,y),y,[ys_13_240;0;ys_13_240;0]));
dJ=factor(det(Js_13_24));
% find first positive point with det(J)=0
csol=array_solve(dJ,r);
ind=find(cellfun(@(c)isAlways(c>1),csol),1,'first');
%% bifurcation of y1=y3,y2=y4=0 equilibrium at r>0
rs=csol{ind}; % bifurcation value of r
y1bif=simplify(subs(ys_13_240,r,rs)); % bifurcation value of y1 and y3 (y2=y4=0)
fprintf('Bifurcation with r=%s, y1=%s\n',char(rs),char(y1bif));
%% nullspace of Jacobian in bifurcation (only y component)
Jbif=subs(Js_13_24,r,rs);
nspace=null(Jbif);
fprintf('Nullspace at r=%s: {%s, %s}\n',char(rs),char(nspace(:,1)),char(nspace(:,2)));
%% Look for equilibria with y4=0, y1=y3
% defining equations
y_eqn_13_40=numden(simplify(subs(y_eqn(1:2)./y(1:2),y(3:4),[y(1);0])));
info13_40=reduce_grobner(y_eqn_13_40,y(1:2),[r,g]);
[info13_40.target,info13_40.factor,info13_40.source]=deal(...
           [3,4],           [1,0],           [1,1]);
y13_40=arrayfun(@(r){get_eq_grobner(r,g_n,info13_40,true)},r_rgfine);
y13_40=permute(cat(3,y13_40{:}),[1,3,2]);
[x13_40,par13_40,stab13_40]=x_eq_from_y(r_rgfine,y13_40,A,par0,in,FJ);
plt13_40args=[repmat({r_rgfine},size(y13_40,3),1),squeeze(num2cell(y13_40(1,:,:),[1,2]))].';
figure(1);clf;ax=gca;hold(ax,'on');
plot(ax,[1,r_rgfine],yfun_13_240([1,r_rgfine],g_n));
plot(ax,plt13_40args{:})
%% check for more general equilibria of the form y1==y3
y_eqn_13=subs(y_eqn([1,2,4])./y([1,2,4]),y(3),y(1));
info13=reduce_grobner(y_eqn_13,y([1,2,4]),[r,g]);
[info13.target,info13.factor,info13.source]=deal(...
        [4,3],          [1,1],        [3,1]);
y13=arrayfun(@(r){get_eq_grobner(r,g_n,info13,true)},r_rgfine);
y13=permute(cat(3,y13{:}),[1,3,2]);
[x13,par13,stab13]=x_eq_from_y(r_rgfine,y13,A,par0,in,FJ);
plt13args=[repmat({r_rgfine},size(y13,3),1),squeeze(num2cell(y13(1,:,:),[1,2]))].';
plot(ax,plt13args{:})
nr=length(r_rgfine);
rs1=@(z)reshape(z,nr,[]);
rs2=@(z)reshape(z,size(z,1),nr,[]);
% separate (1234) out from (13)
sel1234=squeeze(all(diff(x13(1:4,:,:),[],1)<1e-12));
[x1234,par1234,stab1234]=deal(rs2(x13(:,sel1234)),par13,rs1(stab13(sel1234)));
[x13,stab13]=deal(rs2(x13(:,~sel1234)),rs1(stab13(~sel1234)));
%% check for general equilibria with y4=0
y_eqn_40g=subs(y_eqn([1,2,3])./y([1,2,3]),y(4),0);
info40g=reduce_grobner(y_eqn_40g,y([1,2,3]),[r,g]);
[info40g.target,info40g.factor,info40g.source]=deal(...
         4,              0,             1);
y40g=arrayfun(@(r){get_eq_grobner(r,g_n,info40g,true)},r_rgfine);
y40g=permute(cat(3,y40g{:}),[1,3,2]);
[x40g,par40g,stab40g]=x_eq_from_y(r_rgfine,y40g,A,par0,in,FJ);
plt40gargs=[repmat({r_rgfine},size(y40g,3),1),squeeze(num2cell(y40g(1,:,:),[1,2]))].';
plot(ax,plt40gargs{:})
%% check for general equilibria with y3=y4
y_eqn_34=subs(y_eqn([1,2,3])./y([1,2,3]),y(4),y(3));
info34=reduce_grobner(y_eqn_34,y([1,2,3]),[r,g]);
[info34.target,info34.factor,info34.source]=deal(...
      4,              1,             3);
y34=arrayfun(@(r){get_eq_grobner(r,g_n,info34,true)},r_rgfine);
y34=permute(cat(3,y34{:}),[1,3,2]);
[x34,par34,stab34]=x_eq_from_y(r_rgfine,y34,A,par0,in,FJ);
% turns out all these equilibria have y1=y2=y3=y4

%%
save('stst_sym_grobner.mat')
%%
function dpoly=poly_deg(expr,y)
dpoly = -ones(size(expr));
for i=1:numel(expr)
    try
        dpoly(i)=polynomialDegree(expr(i),y);
    end
end
end
%%
%%
function sol=array_solve(expr,y)
sol = num2cell(NaN(size(expr)));
for i=1:numel(expr)
    try
        sol{i}=solve(expr(i),y);
        if isempty(sol{i})
            sol{i}=sym(NaN);
        end
    end
end
end