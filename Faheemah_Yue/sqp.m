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
  
  format_sortie = ""; % affichage sur l'écran
  
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
  
  %% la chaîne initiale
  figure(1); clf(1) ;
  simul(1, x, lme, lmi); % tracer la chaîne
  
  [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
  [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
  grad_lag = g + ae'*lme + ai'*kron(lmi, ones(p,1));
  [L, d, flag] = cholmod(hl, small, big);
  M = L*diag(d)*L';
  
  x_h=[x]; % historique des x
  lme_h=[lme]; % historique des lme
  lmi_h=[lmi]; % historique des lme
  
  fprintf(FOUT, "------------------------------------------- \n "); 
  fprintf(FOUT, "iter\t\t |gl|\t\t |ce|\t\t (ci,lmi) \n "); 
  
  %% algo quasi-Newthon
  while 1
    if p==0
      fprintf(FOUT, '%d\t %e\t %e\t %e\t \n', info.niter, norm(grad_lag, Inf), norm(ce, Inf), 0); 
    else
      fprintf(FOUT, '%d\t %e\t %e\t %e\t \n', info.niter, norm(grad_lag, Inf), norm(ce, Inf), norm(min(min(kron(lmi, ones(p,1)), -ci)), Inf)); 
    endif
    
    %fprintf(FOUT, '%d\t %e\t %e\t %e\t \n', info.niter, norm(grad_lag, Inf), norm(ce, Inf), lmi'*ci); 
    
    %% sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(ce, Inf) <= options.tol(2) && norm(min(min(lmi, -ci)), Inf) <= options.tol(3)
      info.status = 0; % terminaison normale
      printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n On trouve une solution minimum globale\n "); 
      break
    endif
    
    %% calculer la solution primale-duale par OQS
    [d,lm,infoqp] = qpalm(g, M, repelem(-inf, n), repelem(inf, n), ai, repelem(-inf, mi*p), -ci, ae, -ce, [ ], [repelem(0, n)' kron(lmi, ones(p,1))' lme']', struct()); 
    x = x + d;
    lme = lm(n+1+p*mi:n+p*mi+me);
    lmi = lm(n+1:n+p*mi);
    if p>=2
      lmi = lmi(1:size(ones(p,1),1):end,1:size(ones(p,1),2):end)./ones(p,1)(1);
    endif
    
    %% calculer la vouvelle itération
    [e, ce, ci, g, ae, ai, hl, indic] = simul(5, x, lme, lmi); % calcul hessien de lagrangien
    [e, ce, ci, g, ae, ai] = simul(4, x, lme, lmi); % calcul ce, ci, g, ae, ai
    grad_lag = g + ae'*lme + ai'*kron(lmi, ones(p,1));
    [L, d, flag] = cholmod(hl, small, big); 
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

  %% la chîne finale
  figure(2); clf(2) ;
  simul(1, x, lme, lmi); % tracer la chaîne

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


