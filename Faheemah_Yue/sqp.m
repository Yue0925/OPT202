function [x, lme, lmi, info] = sqp(simul, x, lme, lmi, options)
  global R
  
  %% initializer les valeurs sorties
  info.status = 0; %sortie normale
  info.niter = 1; %nb itérations
  small = 1.e-5;
  big = 1.e+5;

  p = length(R);
  n = length(x);
  me = length(lme);
  mi = length(lmi);
  
  %% vérification vars entrées
  if n==0
    error("ERROR: valeur initial x est vide. \n");
    info.status = 1; %inconsistance des arguments d’entrée
    return
  endif

  if options.maxit<1
    error("ERROR: le maximum d'itération est inférieur à 1 ici. \n");
    info.status = 1; %inconsistance des arguments d’entrée
    return
  endif
  
  %% les données initiales
  printf("\n Afficher les données initiales\n");
  x
  lme
  lmi
  figure(1); clf(1) ;
  simul(1, x, lme, lmi); % tracer la chaîne
  
  [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
  [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
  grad_lag = g + ae'*lme + ai'*kron(lmi, ones(p,1));%% A REFLECHIR!!!!!!!!!! p*nn nn*p
  [L, d, flag] = cholmod(hl, small, big);  %%%%%%%%%%%%%%
  M = L*diag(d)*L';
  
  x_h=[x]; % historique des x
  lme_h=[lme]; % historique des lme
  lmi_h=[lmi]; % historique des lme
  
  %% algo quasi-Newthon
  while 1
    
    %% sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(ce, Inf) <= options.tol(2) && norm(min(min(lmi, -ci)), Inf) <= options.tol(3)
      info.status = 0; % terminaison normale
      printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n "); 
      printf("On trouve un point stationnaire z_*, voici son hessien lagrangien "); 
      break
    endif
    
    %% calculer la solution primale-duale par OQS
    [d,lm,infoqp] = qpalm(g, M, repelem(-inf, n), repelem(inf, n), ai, repelem(-inf, mi*p), -ci, ae, -ce, [ ], [repelem(0, n)' kron(lmi, ones(p,1))' lme']', struct()); %%%%%%%%%%%%%%%%%%
    x = x + d;
    lme = lm(n+1+p*mi:n+p*mi+me);
    lmi = lm(n+1:n+p*mi);
    if p>=2
      lmi = lmi(1:size(ones(p,1),1):end,1:size(ones(p,1),2):end)./ones(p,1)(1);
    endif
    
    %% calculer la vouvelle itération
    [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
    [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
    grad_lag = g + ae'*lme + ai'*kron(lmi, ones(p,1)); %% A REFLECHIR!!!!!!!!!! p*nn nn*p
    [L, d, flag] = cholmod(hl, small, big); %%%%%%%%%%%%%%
    M = L*diag(d)*L';   
    
    info.niter = info.niter + 1;
    
    x_h=[x_h, x]; % historique des x
    lme_h=[lme_h, lme]; % historique des lme
    lmi_h=[lmi_h, lmi]; % historique des lme
    
    %% sortir si max itération atteint
    if info.niter >= options.maxit
      info.status = 2;
      break
    endif
    
  endwhile

  figure(2); clf(2) ;
  simul(1, x, lme, lmi); % tracer la chaîne

  printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n "); 
  printf("Affichage les résultats historiques:\n "); 
  x_h
  lme_h 
  lmi_h 

  return
endfunction


