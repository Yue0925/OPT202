% Prog principal
global A B L R S

%% les test de TP4
%test4a;
%test4b;
test4c;


lme = [];
lmi = [];


%% calculer lambda lme initial par moindre carrées
[e, ce, ci, g, ae, ai, hl, indic] = chs(4, xy, lme, lmi);
lme = -ae'\g; 
%% 
lmi = repelem(0,length(xy)/2);

%%
%indic = 5;
%[e, ce, ci, g, ae, ai, hl, indic] = chs(indic, xy, lme, lmi)


%---------------------------------------------
%small = 1.e-5;
%big = 1.e+5;

%[L, d, flag] = cholmod(hl, small, big)
%L*diag(d)*L'

options.tol = [1.e-6 1.e-6 1.e-6]; %précision
options.maxit = 100; %nb max itérations

%options.rl = 0; % 0: avec recherche lineair; 1: le pas unitaire
%options.verb = 1; %1: les sorties de chaque itération; 2: détailées de la boucle de RL sur la recherche de alpha

[x, lme, lmi, info] = sqp(@chs, xy, lme, lmi, options);
info

