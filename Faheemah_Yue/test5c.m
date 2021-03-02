
%% cas-test 5b
A = 0;
B = 0;
L = [0.1 0.2 0.3 0.4 0.5 0.4 0.3 0.1]'; % longeur des barres fixées

xy = [-0.1   0   -0.2  -0.1    0.2    0.1   0.2   ... 
      0.1   0.3   0.5   0.9    1.4     1   0.7]'; 
      
R = [-1.0; -0.2; -1.0];
S = [-7.0; 0.0; 7.0];

FIG_XMIN=-1; FIG_XMAX=1; FIG_YMIN=-1; FIG_YMAX=2;

% option Méthode Newton : 
% option Quasi Newton : 