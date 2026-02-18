clear
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(fullfile('ddebiftool_snapshot'));
format compact
format short g
load('steady_state_y2_0.mat')

%% Parameters (from Figure 4 in the paper)
beta  = 40;
sigma = 10;
mu    = 0.02;
gamma = 0.8; %0.58 0.6 0.8 %  figure4
%%
par(in.beta)=40;
par(in.sigma)=10;
par(in.mu)=0.2;
par(in.gamma)=0.53;
par(in.r)=beta/(mu + sigma);
par(in.e) = mu/(mu + sigma);
r=par(in.r);
e=par(in.e);
%% Fully symmetric equilibrium E = (Y,Y,Y,Y, Z,Z,Z,Z, W,W,W,W)
% From eq. (7): Y
term   = r*(3 - 2*gamma) - 4;
Y = (e/(6*r)) * ( term + sqrt(term^2 + 12*(r - 1)) );

% From eq. (6): Z and W
Z = (r*Y) / (r*Y + e);
W = (3*r*Y) / (3*r*Y + e);

E = [Y;Y;Y;Y;  Z;Z;Z;Z;  W;W;W;W];   % 12x1

%% 
eps0 = 1e-8;  
vOdd = [ 1; -1; 1; -1;  0;0;0;0;  0;0;0;0 ];  % perturb y's only

x0 = E + eps0*vOdd;

disp("Equilibrium E:");
disp(E.');
disp("Initial guess x0 (perturbed E):");
disp(x0.');
%%
Fode45=@(t,X)funcs.wrap_rhs(X,par)
%% Example call (your Fode45 should accept (t,x))
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'Events', @(t,x) MyEvent(t,x));
sol = ode45(Fode45, [0 1e5], x0, opts);
%%

figure(88)
clf;
% tiledlayout(1,3)
% nexttile
hold on
plot(sol.x,sol.y(5:8,:),'LineWidth',2)


%%
[funcs_s2,po2_symmetry,suc2]=branch_from_ode45sol(fun_eq13_24,sol,in.r,par,...
    'phase_condition',0,'matrix','sparse','stability.eigmatrix','sparse',parbd{:},...
    'degree',6,'intervals',80,'newton_max_iterations',8,...
    'adapt_mesh_after_correct',1,...
    'usercond',[],'outputfuncs',true);

