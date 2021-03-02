% Prog principal
global A B L R S FOUT FIG_XMIN FIG_XMAX FIG_YMIN FIG_YMAX

%% les test de TP4
%test4a; FOUT = fopen('test4a.res','w');	
%test4b; FOUT = fopen('test4b.res','w');	
%test4c; FOUT = fopen('test4c.res','w');	

%% les tests de TP5
%test5a; FOUT = fopen('test5a.res','w');	
test5b; FOUT = fopen('test5b.res','w');	


lme = [];
lmi = [];

%% calculer lambda lme initial par moindre carrées
[e, ce, ci, g, ae, ai, hl, indic] = chs(4, xy, lme, lmi);
lme = -ae'\g; 
%% calculer lambda lmi initial
lmi = -ai'\g; 

options.tol = [1.e-6 1.e-6 1.e-6]; %précision
options.maxit = 100; %nb max itérations
options.deriv = 2; % 1: algo Quasi-Newton; 2: méthode de Newton


[x, lme, lmi, info] = sqp(@chs, xy, lme, lmi, options);
info
fclose(FOUT);
