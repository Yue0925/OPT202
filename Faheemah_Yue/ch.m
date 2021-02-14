% Prog principal
global A B L R S FOUT

%% les test de TP4
%test4a; FOUT = fopen('test4a.res','w');	
%test4b; FOUT = fopen('test4b.res','w');	
test4c; FOUT = fopen('test4c.res','w');	


lme = [];
lmi = [];


%% calculer lambda lme initial par moindre carrées
[e, ce, ci, g, ae, ai, hl, indic] = chs(4, xy, lme, lmi);
lme = -ae'\g; 
%% calculer lambda lmi initial
lmi = repelem(0,length(xy)/2);

options.tol = [1.e-6 1.e-6 1.e-6]; %précision
options.maxit = 100; %nb max itérations

%options.rl = 0; % 0: avec recherche lineair; 1: le pas unitaire
%options.verb = 1; %1: les sorties de chaque itération; 2: détailées de la boucle de RL sur la recherche de alpha

[x, lme, lmi, info] = sqp(@chs, xy, lme, lmi, options);
info
fclose(FOUT);
