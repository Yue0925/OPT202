function [e, ce, ci, g, ae, ai, hl, indic] = chs(indic, xy, lme, lmi)
  
  global R S
  global EXIT_SUCESS = 0; % sortie normale 
  global EXIT_FAILURE = 1; % paramètre(s) d’entrée non correct(s)
  
  %% initialisation des output de fonction
  e = [];
  ce = [];
  ci = [];
  g = [];
  ae = [];
  ai = [];
  hl = [];
  
  %% pré-calculs et vérifier les données en entrées
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  %% vérification la bonne entrées
  if mod(length(xy),2)==1
    error("ERROR: La longeur de xy d'entrée est impaire.\n");
    indic = EXIT_FAILURE; % length(xy) impaire
    return
  endif
  
  %% vérification paramètres des planchers
  if length(R) ~= length(S)
    error("ERROR: La longeur de R et S d'entrée est différente.\n");
    indic = EXIT_FAILURE;
    return
  endif
  
  %% pilote indic
  switch (indic)
    
    case 1 % tracer la chaîne
      plot(x_complet, y_complet, '-b') % tracé de la chaîne
      indic = EXIT_SUCESS;
      return
      
    case 2 % calcul de e, ce, ci
      [e] = cal_e(xy);
      [ce] = cal_ce(xy);
      [ci] = cal_ci(xy);
      indic = EXIT_SUCESS;
      return
      
    case 4 % calcul de e, ce, ci, g, ae, ai
      [e] = cal_e(xy);
      [ce] = cal_ce(xy);
      [ci] = cal_ci(xy);
      [g] = cal_g(xy);
      %verify_gradient(xy);
      [ae] = cal_ae(xy);
      [ai] = cal_ai(xy);
      %verify_Jacobian(xy);
      indic = EXIT_SUCESS;
      return

    case 5 % calcul de hl hessien de lagrangien
      [hl] = calcul_hl(xy,lme);  % calcul de hl
      %verify_Hess_Lag(xy, lm);
      indic = EXIT_SUCESS;
      return
      
    otherwise
      error("ERROR: invalid indic number %d.\n ", indic);
      indic = EXIT_FAILURE;
      return
   endswitch
endfunction



%% pre-calcul / traitements les donées en entrées
function [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy)
  global A B R
  
  nn = length(xy)/2; % nombre noeuds(sauf extrémités)
  nb = nn+1; % nombre barres
  p = length(R);
  
  x = xy(1:nn); 
  y = xy(nn+1 : end);
  
  x_complet = [0;x;A]; % les abscisses en ajoutant les deux extrémités
  y_complet = [0;y;B]; % les ordonnées en ajoutant les deux extrémités
endfunction



%% fonction qui calcule e énergie potentielle
function [e] = cal_e(xy)
  global L
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  e = sum((y_complet(2:nb+1) + y_complet(1:nn+1))./2 .* L(1:nb)); % énergie
endfunction



%% fonction qui calcule ce les contraintes d'égalités
function [ce] = cal_ce(xy)
  global L
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  ce = (x_complet(2:nb+1) - x_complet(1:nb)).^2 + (y_complet(2:nb+1) - y_complet(1:nb)).^2 - L(1:nb).^2; % contrainte
endfunction



%% fonction qui calcule ci les contraintes d'inégalités
function [ci] = cal_ci(xy)
  global L R S
  ci = [ ];
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  for i=1:nn
    ci = [ci, [R + x(i).*S - y(i)]'];
  endfor
  ci = ci';
endfunction



%% fonction qui calcule g le gradient
function [g] = cal_g(xy)
  global L
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  g = [zeros(nn,1); (L(1:nn)+L(2:nb))/2]; % gradient de xy
endfunction



% fonction qui calcule a Jacobians des contraintes d'égalitées
function [ae] = cal_ae(xy)
  global A B
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  ax = spdiags([2*(x_complet(2:nb) - x_complet(3:nb+1)) 2*(x_complet(2:nb) - x_complet(1:nn))], -1:0, nb, nn);
  ay = spdiags([2*(y_complet(2:nb) - y_complet(3:nb+1)) 2*(y_complet(2:nb) - y_complet(1:nn))], -1:0, nb, nn);
  ae = [ax ay];
endfunction



% fonction qui calcule a Jacobians des contraintes d'inégalitées
function [ai] = cal_ai(xy)
  global A B S
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);

  ai = sparse(p*nn, 2*nn);
  if p==0
    return
  endif
  
  for i=1:nn
    ai((i-1)*p+1:p*i, i) = S;
    ai((i-1)*p+1:p*i, i+nn) = repelem(-1,p);
  endfor 
endfunction



% calcul de la matrice hessien de lagrangien
function [hl] = calcul_hl(xy,lme)
  [nn, nb, p, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  v1 = 2.*(lme(1:nn) + lme(2:nb));
  v2 = -2.*lme(2:nn);
  
  hl = spdiags([[v2' 0 v2' 0]' [v1' v1']' [0 v2' 0 v2']'], -1:1, 2*nn, 2*nn); 
  return
endfunction



% calcul dérivées de i-ième composant
function [df, t_i] = cal_df(phi, i, xy)
  t_i = sqrt(eps) * max(1, abs(xy(i)));
  e_i = double(1:length(xy) == i); % vector ligne
  df = (phi(xy + (t_i*e_i)') - phi(xy - (t_i*e_i)') ) / (2 * t_i); % d'abord transformer en vec colonne
  return
endfunction



% Vérification gradient affichage un composant par un
function verify_gradient(xy)
  printf("Vérification vecteur gradient \n");
  printf("i \t pas \t \t f'(i) \t \t DF \t \t erreur \t \t < 10e-8? \n");
  g = cal_g(xy); % gradient
  for i=1:length(xy)
    [df, t_i] = cal_df(@cal_e, i, xy);
    printf("%d \t %e \t", i, t_i);
    printf("%e \t %e \t ", g(i), df);
    erreur = abs(df - g(i)); % erreur absolue
    if abs(g(i))!=0
      erreur = erreur/abs(g(i)); % erreur relaive
      printf("%e rel \t", erreur);
    else
      printf("%e abs \t", erreur);
    endif
    if erreur < 10e-8
      printf("TRUE \n");
    else
      printf("FALSE \n");
    endif
  endfor
endfunction



% calcul Jacobians
function [jacobian] = cal_Jacob(phi, xy)
  lin = length(phi(xy));
  col = length(xy);
  jacobian = sparse(lin, col);

  for i=1:col
    [jacobian(1:end, i), t_i] = cal_df(phi, i, xy);
  endfor
  return
endfunction



function verify_Jacobian(xy)
  printf("Vérification Jacobians \n");
  [jacobian] = cal_Jacob(@cal_c, xy);
  [a] = cal_a(xy);
  erreurs = abs(jacobian - a) ./ abs(a);
  %erreurs
  printf("erreurs relatives < 10e-8? (1:TRUE; 0:FALSE) \n");
  erreurs < 10e-8
endfunction



global I=0;

% calcul hessien de lagrangien de i-ème composant
% on vérifie pour chaque composant de Jacob son n composants (en supposant dérivées premiers correctes)
function [hess_lag] = cal_Hess_Lag(xy, lm)
  global I
  nn = length(xy);
  hess_lag = sparse(nn,nn);
  
  for j=1:length(lm)
    I = j;
    tmp = sparse(nn,nn);
    for i=1:nn
      [HL_compo, t_i] = cal_df(@tmp, i, xy);
      tmp(1:end,i) = HL_compo;
    endfor
    hess_lag = hess_lag + tmp *lm(j);
  endfor
  %hess_lag
  return
endfunction



%i-ème ligne du Jacobians
function [comp_i] = tmp(xy)
  global I
  [jacobian] = cal_Jacob(@cal_c, xy);
  comp_i = jacobian(I,1:end)';
endfunction



function verify_Hess_Lag(xy, lm)
  [hess_lag] = cal_Hess_Lag(xy, lm);
  [hl] = calcul_hl(xy,lm);
  erreur = abs(hess_lag - hl)/abs(hl);
  %hess_lag
  erreur
  erreur < 10e-8
endfunction


