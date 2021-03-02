function [x, lme, lmi, info] = sqp(simul, x, lme, lmi, options)
  global R FOUT
  
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
  
  %% la chaîne initiale
  figure(1); clf(1) ;
  simul(1, x, lme, lmi); % tracer la chaîne
  
  [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
  [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
  grad_lag = g + ae'*lme + ai'*lmi;
  if options.deriv==1 % option Quasi-Newton
    M= speye(n,n); % M_1 matrice d'identité
  else
    [L, d, flag] = cholmod(hl, small, big);
    M = L*diag(d)*L';   
  endif
  
  x_h=[x]; % historique des x
  lme_h=[lme]; % historique des lme
  lmi_h=[lmi]; % historique des lme
  
  fprintf(FOUT, "------------------------------------------- \n"); 
  fprintf(FOUT, "iter\t\t |gl|\t\t |ce|\t (ci,lmi)\t |x|\t\t |lm|\t\t Powell\t\t cond(M) \n"); 
  
  %% algo
  while 1
    %% sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(ce, Inf) <= options.tol(2) && norm(min(min(lmi, -ci)), Inf) <= options.tol(3) %lmi
      info.status = 0; % terminaison normale
      printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n On trouve une solution minimum globale\n "); 
      break
    endif
    
    if p==0
      fprintf(FOUT, '%d\t %e\t %e\t %e\t %e\t %e\t ', info.niter, norm(grad_lag, Inf), norm(ce, Inf), 0, norm(x, Inf), norm([lme' lmi']', Inf)); 
    else
      fprintf(FOUT, '%d\t %e\t %e\t %e\t %e\t %e\t ', info.niter, norm(grad_lag, Inf), norm(ce, Inf), norm(min(min(lmi, -ci)), Inf), norm(x, Inf), norm([lme' lmi']', Inf)); 
    endif
    
    %% calculer la solution primale-duale par OQS
    [d,lm,infoqp] = qpalm(g, M, repelem(-inf, n), repelem(inf, n), ai, repelem(-inf, mi), -ci, ae, -ce, [ ], [repelem(0, n)' lmi' lme']', struct()); 
    lmePQ = lm(n+1+mi:n+mi+me);
    lmiPQ = lm(n+1:n+mi);    
    
    %% calculer la vouvelle itération
    if options.deriv==2 % option Méthode Newton
      x = x + d;
      lme = lmePQ;
      lmi = lmiPQ;
      [L, d, flag] = cholmod(hl, small, big); 
      M = L*diag(d)*L';
    else % 1: Quasi-Newton avec recherche linéaire
      x_old = x;
      x = x + d;
      lme = lme + (lmePQ - lme);
      lmi = lmi + (lmiPQ - lmi);
      delta = d;
      % calcul BFGS
      gamma_gl = calcul_gamma_gl(x, lme, lmi, simul, x_old);
      % correction de Powell
      if gamma_gl'*delta < 0.2*delta'*M*delta
        theta = 0.8*( (delta'*M*delta)/(delta'*M*delta - gamma_gl'*delta));
      else
        theta = 1;
      endif
      gamma = (1-theta) * M * delta + theta * gamma_gl;
      % mis à jour M
      if info.niter == 1
        eta = norm(gamma, 2)^2 / (gamma'*delta);
        M = M*eta; % motif M_1
      endif
      M = M - ( (M*delta*delta'*M)/(delta'*M*delta) ) + ( (gamma*gamma')/(gamma'*delta) );
    endif
    
    if options.deriv== 1 % Quasi-Newton
      fprintf(FOUT, '%e\t %e\t \n', theta, cond(M));
    else
      fprintf(FOUT, '\t\t %e\t \n', cond(M));
    endif
    
    [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
    [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
    grad_lag = g + ae'*lme + ai'*lmi;
   
    x_h=[x_h, x]; % historique des x
    lme_h=[lme_h, lme]; % historique des lme
    lmi_h=[lmi_h, lmi]; % historique des lme
    
    info.niter = info.niter + 1;
    
    %% sortir si max itération atteint
    if info.niter >= options.maxit
      info.status = 2;
      break
    endif
  endwhile

  %% la chîne finale
  figure(2); clf(2) ;
  simul(1, x, lme, lmi); % tracer la chaîne
  plancher;

  fprintf(FOUT, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"); 
  fprintf(FOUT, "Affichage les résultats historiques:\n"); 
  
  for i = 1:info.niter
    fprintf(FOUT, 'x%d: [', i);
    fprintf(FOUT, '%g, ', x_h(1:end)((i-1)*n+1:n*i-1));
    fprintf(FOUT, '%g]\n', x_h(1:end)(n*i));
  endfor
  fprintf(FOUT, "\n"); 
  for i = 1:info.niter
    fprintf(FOUT, 'lme%d: [', i);
    fprintf(FOUT, '%g, ', lme_h(1:end)((i-1)*me+1:me*i-1));
    fprintf(FOUT, '%g]\n', lme_h(1:end)(me*i));
  endfor
  if p==0
    return
  endif
  fprintf(FOUT, "\n"); 
  for i = 1:info.niter
    fprintf(FOUT, 'lmi%d: [', i);
    fprintf(FOUT, '%g, ', lmi_h(1:end)((i-1)*mi+1:mi*i-1));
    fprintf(FOUT, '%g]\n', lmi_h(1:end)(mi*i));
  endfor

endfunction


function gamma_gl = calcul_gamma_gl(x, lme, lmi, simul, x_old)
  [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
  gl1 = g + ae'*lme + ai'*lmi; 
  [e, ce, ci, g, ae, ai] = simul(4, x_old, lme, lmi); 
  gl2 = g + ae'*lme + ai'*lmi;
  gamma_gl = gl1 - gl2;
  return
endfunction



