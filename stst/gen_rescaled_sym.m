clear;
clear
startup_coco(fullfile('..','coco_2025January28'))
parnames = {'r','e','g'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
syms(parnames{:});
par = sym(parnames);
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
dydt=Lambda.*(1-g*w-(1-g).*z)-y; %% need to check the term gamma*w or just gamma
dzdt=e*(Lambda.*(1-z)-z);
dwdt=e*((1-w).*AL-w);
%%
F=sco_sym2funcs([dydt;dzdt;dwdt],{[y;z;w],par(:)},{'x','p'},'filename','dengue_coco');
%%


