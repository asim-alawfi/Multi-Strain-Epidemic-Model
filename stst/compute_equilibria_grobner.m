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
%% for general equilibrium: w and z are functions of y -> substitute
ws=solve(dwdt,w);
zs=solve(dzdt,z);
y_eqn=expand(subs(subs(dydt,zs),ws));
%% mesh of r values at which we compute equilibria
r_rgfine=linspace(1.001,3.5,3500);
%% trivial equilibria
r00=linspace(0.5,3.5,4001);
y00=zeros(nstrains,length(r00));
[x00,par00,stab00]=x_eq_from_y(r00,y00,A,par0,in,FJ);
plt00args=[repmat({r00},size(y00,3),1),squeeze(num2cell(y00(1,:,:),[1,2]))].';
figure(1);clf;ax=gca;hold(ax,'on');
plot(ax,plt00args{:});drawnow
%% fully symmetric equilibria, y1=y2=y3=y4>0
y_eqn_1234=subs(simplify(y_eqn(1)./y(1)),y(2:end),y(1)*ones(nstrains-1,1));
y_eqn_1234=numden(simplify(y_eqn_1234));
info1234=reduce_grobner(y_eqn_1234,y(1),[r,g]);
[info1234.target,info1234.factor,info1234.source]=deal(...
         [2,3,4],        [1,1,1],         [1,1,1]);
y1234=arrayfun(@(r){get_eq_grobner(r,g_n,info1234,true)},r_rgfine);
y1234=permute(cat(3,y1234{:}),[1,3,2]);
[x1234,par1234,stab1234]=x_eq_from_y(r_rgfine,y1234,A,par0,in,FJ);
[x1234,par1234,stab1234]=br_remove_sym(x1234,par1234,stab1234,...
    @(x)abs(x(1:4,:,:))<1e-12);
plt1234args=[repmat({r_rgfine},size(y1234,3),1),squeeze(num2cell(y1234(1,:,:),[1,2]))].';
plot(ax,plt1234args{:});drawnow
%% y1=y3,y2=y4=0 equilibria and their bifurcation with non-zero y1
y_eqn_13_240=numden(simplify(subs(simplify(y_eqn(1)./y(1)),y(2:4),[0;y(1);0])));
y_eqn_13_240=numden(simplify(y_eqn_13_240));
info13_240=reduce_grobner(y_eqn_13_240,y(1),[r,g]);
[info13_240.target,info13_240.factor,info13_240.source]=deal(...
           [2,3,4],          [0,1,0],          [1,1,1]);
y13_240=arrayfun(@(r){get_eq_grobner(r,g_n,info13_240,true)},r_rgfine);
y13_240=permute(cat(3,y13_240{:}),[1,3,2]);
[x13_240,par13_240,stab13_240]=x_eq_from_y(r_rgfine,y13_240,A,par0,in,FJ);
plt13_240args=[repmat({r_rgfine},size(y13_240,3),1),squeeze(num2cell(y13_240(1,:,:),[1,2]))].';
plot(ax,plt13_240args{:});drawnow
%% bifurcation of y1=y3,y2=y4=0 equilibrium at r>0
ys_13_240=solve(y_eqn_13_240(1),y(1)); % find non-zero solution y
fprintf('Equilibrium with y1=y3,y2=y4=0: y1=%s\n',char(ys_13_240));
Js_13_24=simplify(subs(jacobian(y_eqn,y),y,[ys_13_240;0;ys_13_240;0]));
dJ=factor(det(Js_13_24));
% find first positive point with det(J)=0
csol=array_solve(dJ,r);
ind=find(cellfun(@(c)isAlways(c>1),csol),1,'first');
rs=csol{ind}; % bifurcation value of r
y1bif=simplify(subs(ys_13_240,r,rs)); % bifurcation value of y1 and y3 (y2=y4=0)
fprintf('Bifurcation with r=%s, y1=%s\n',char(rs),char(y1bif));
%% nullspace of Jacobian in bifurcation (only y component)
Jbif=subs(Js_13_24,r,rs);
fprintf('Nullspace at r=%s:\n',char(rs));
nspace=null(Jbif)
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
plot(ax,plt13_40args{:});drawnow
%% check for more general equilibria of the form y1==y3
y_eqn_13=subs(y_eqn([1,2,4])./y([1,2,4]),y(3),y(1));
info13=reduce_grobner(y_eqn_13,y([1,2,4]),[r,g]);
[info13.target,info13.factor,info13.source]=deal(...
        [4,3],          [1,1],        [3,1]);
y13=arrayfun(@(r){get_eq_grobner(r,g_n,info13,true)},r_rgfine);
y13=permute(cat(3,y13{:}),[1,3,2]);
[x13,par13,stab13]=x_eq_from_y(r_rgfine,y13,A,par0,in,FJ);
% separate (1234) and (13)4=0 and (13)2=0 out from (13)
[x13,par13,stab13]=br_remove_sym(x13,par13,stab13,...
    @(x)abs(diff(x(1:4,:,:)))<1e-12|abs(x(4,:,:))<1e-12|abs(x(2,:,:))<1e-12);
plt13args=[repmat({r_rgfine},size(y13,3),1),squeeze(num2cell(y13(1,:,:),[1,2]))].';
plot(ax,plt13args{:});drawnow
%% check for general equilibria with y4=0 (no extra solutions exist)
y_eqn_40g=subs(y_eqn([1,2,3])./y([1,2,3]),y(4),0);
info40g=reduce_grobner(y_eqn_40g,y([1,2,3]),[r,g]);
[info40g.target,info40g.factor,info40g.source]=deal(...
         4,              0,             1);
y40g=arrayfun(@(r){get_eq_grobner(r,g_n,info40g,true)},r_rgfine);
y40g=permute(cat(3,y40g{:}),[1,3,2]);
[x40g,par40g,stab40g]=x_eq_from_y(r_rgfine,y40g,A,par0,in,FJ);
% separate (13)4=0 and (13)2=0 out from 4=0
[x40g,par40g,stab40g]=br_remove_sym(x40g,par40g,stab40g,...
    @(x)abs(diff(x([1,3],:,:),[],1))<1e-12|abs(x(2,:,:))<1e-12);
assert(isempty(x40g));
%% check for general equilibria with y3=y4
y_eqn_34=subs(y_eqn([1,2,3])./y([1,2,3]),y(4),y(3));
info34=reduce_grobner(y_eqn_34,y([1,2,3]),[r,g]);
[info34.target,info34.factor,info34.source]=deal(...
      4,              1,             3);
y34=arrayfun(@(r){get_eq_grobner(r,g_n,info34,true)},r_rgfine);
y34=permute(cat(3,y34{:}),[1,3,2]);
[x34,par34,stab34]=x_eq_from_y(r_rgfine,y34,A,par0,in,FJ);
% separate (1234) (34)
[x34,par34,stab34]=br_remove_sym(x34,par34,stab34,...
    @(x)abs(diff(x(1:4,:,:),[],1))<1e-12);
assert(isempty(x34));
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