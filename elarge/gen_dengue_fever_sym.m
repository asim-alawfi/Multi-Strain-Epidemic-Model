clear;
ddebiftool_path(fullfile('..','ddebiftool_snapshot'),'sym');
%%
format compact
format short g
%%
%% Define parameters
parnames = {'r','e','g'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
syms(parnames{:});
par = sym(parnames);
%% Symbolic setup
[dim,nstrains,ntau]=deal(3,4,0);
y = sym('y', [nstrains, ntau+1]);  % infectious strains y1..y4
z = sym('z', [nstrains, ntau+1]);  % immune responses z1..z4
w = sym('w', [nstrains, ntau+1]);  % cross-immunity w1..w4
%% matrix A
A = [1 1 0 1;
     1 1 1 0;
     0 1 1 1;
     1 0 1 1];
Lambda=r*y;
AL=A*Lambda;
dydt=Lambda.*(1-g*w-(1-g).*z)-y; %% need to check the term gamma*w or just gamma
dzdt=e*(Lambda.*(1-z)-z);
dwdt=e*((1-w).*AL-w);
%%
[fstr,erives]=dde_sym2funcs([dydt;dzdt;dwdt],[y;z;w],par,'filename','dengue_fever_sym',...
    'maxorder',2,'directional_derivative',true);
%%


