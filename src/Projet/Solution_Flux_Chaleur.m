function Solution_flux_chaleur(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot)
    % --- Définition des constantes ---
    % Ces valeurs sont passées directement en paramètres de la fonction
    
    % Appel de la fonction pour obtenir les profils de température numériques
    [T_numerique, ~] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

    % Discrétisation (espacement entre les nœuds)
    dx = Longueur / (Ntot - 1);

    % Application de la différence centrale d'ordre 2 pour le gradient de température à x=0
    % Différence centrale d'ordre 2 : (T_2 - T_(-1)) / (2 * dx)
    T0 = T_numerique(1);  % Température en x=0 (T_0)
    T1 = T_numerique(2);  % Température en x=dx (T_1)
    T2 = T_numerique(3);  % Température en x=2*dx (T_2)
    
    gradT_x0_num = (-3*T0 + 4*T1 - T2) / (2*dx);  % Approximation de la différence centrale

    % Calcul du flux de chaleur à x=0 en utilisant la loi de Fourier
    Ac = (pi() * D^2) / 4;   % Aire de la section transversale de l'ailette
    q_num = -k_cuivre * Ac * gradT_x0_num;  % Flux de chaleur à x=0
    
    % Calcul symbolique du flux de chaleur analytique
    syms x m L
    T_expr = T_inf + (Tm - T_inf) * cosh(m * (L - x)) / cosh(m * L);
    dTdx = diff(T_expr, x);

    % Calcul de m et évaluation de la dérivée symbolique en x=0
    P = pi * D;
    m_val = (h * P) / (k_cuivre * Ac);
    gradT_x0_ana = double(subs(dTdx, [x, m, L], [0, m_val, Longueur]));
    q_ana = -k_cuivre * Ac * gradT_x0_ana;

    % Affichage des résultats
    fprintf('Flux de chaleur numérique à x = 0 : %.3f W\n', q_num);
    fprintf('Flux de chaleur analytique à x = 0 : %.3f W\n', q_ana);
end

