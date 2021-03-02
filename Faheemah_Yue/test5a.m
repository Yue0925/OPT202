
%% cas-test 5a
A = 0;
B = 0;
L = [0.5 0.3 0.4 1.2 0.3 0.3]'; % longeur des barres fixées

xy = [0.2   0.5   0.8   1.0   1.2   ... 
      -0.4  -0.6  -0.4  -0.2   0.]'; 
      
R = [-1];
S = [-0.1];

FIG_XMIN=-1; FIG_XMAX=2; FIG_YMIN=-1.5; FIG_YMAX=2;

% option Quasi-Newton : 22 ite
% option Méthode Newton : 51 ite