clear;
base=[pwd(),filesep(),'..',filesep(),filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');%[getenv('HOME'),'/sourceforge-ddebiftool/releases/git-2023-06-27'],'sym');
%%
format compact
format short g
%%
%% Define parameters
parnames = {'r','e','beta','mu','sigma','gamma','tau_delay'};
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
dydt=Lambda.*(1-w-(1-gamma).*z)-y; %% need to check the term gamma*w or just gamma
dzdt=Lambda.*(1-z)-e*z;
dwdt=(1-w).*AL-e*w;
%%
[fstr,erives]=dde_sym2funcs([dydt;dzdt;dwdt],[y;z;w],par,'filename','dengue_fever_sym',...
    'maxorder',2,'directional_derivative',true);
%%


