function [T_numerique, T_analytique, q_numerique, q_analytique] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot)
    %{
        Cette fonction résout l'équation de diffusion de la chaleur dans une ailette cylindrique
        selon un modèle stationnaire. Elle retourne la solution numérique et la solution analytique
        ainsi que les flux de chaleur numérique et analytique.

        Paramètres :
        D         - Diamètre du pilier (m)
        Longueur  - Longueur du pilier (m)
        k_cuivre  - Conductivité thermique du cuivre (W/m.K)
        h         - Coefficient de transfert thermique (W/m².K)
        T_inf     - Température ambiante (°C)
        Tm        - Température à l’entrée du pilier (°C)
        Ntot      - Nombre total de nœuds pour la discrétisation spatiale
    %}

    % Calcul des paramètres géométriques
    P = pi() * D;                      % Périmètre (m)
    Ac = (pi() * D^2) / 4;             % Aire de section transversale (m²)
    m = (h * P) / (k_cuivre * Ac);     % Coefficient

    % Discrétisation spatiale
    dx = Longueur / (Ntot - 1);

    % Initialisation des matrices A et B
    A = zeros(Ntot-2, Ntot-2); 
    B = zeros(Ntot-2, 1);       

    % Remplissage de la matrice A
    A(1, 1) = (-m^2 * dx^2) - 2; 
    A(1, 2) = 1;

    for i = 2:Ntot-3
        A(i, i-1) = 1;
        A(i, i) = (-m^2 * dx^2) - 2;
        A(i, i+1) = 1;
    end

    A(Ntot-2, Ntot-3) = 2/3;
    A(Ntot-2, Ntot-2) = (4/3) - m^2 * dx^2 - 2;

    % Remplissage de la matrice B
    B(1) = (-m^2 * dx^2 * T_inf) - Tm;
    for i = 2:Ntot-2
        B(i) = -m^2 * dx^2 * T_inf;
    end

    % Résolution du système
    T = A \ B; 

    % Ajout des conditions limites
    T4 = (1/3)*(4*T(end) - T(end-1)); % Gear retardé
    T_numerique = [Tm; T; T4];

    % Coordonnées des nœuds
    x = linspace(0, Longueur, Ntot);

    % Solution analytique
    T_analytique = T_inf + (Tm - T_inf) * cosh(m * (Longueur - x)) / cosh(m * Longueur);

    % --- Flux de chaleur ---
    % Calcul du flux de chaleur numérique à x=0 avec une différence avant de premier ordre
    T0 = T_numerique(1);  % Température en x=0 (T_0)
    T1 = T_numerique(2);  % Température en x=dx (T_1)
    
    % Approximation de la différence avant (forward difference) de premier ordre
    gradT_x0_num = (T1 - T0) / dx;  % First-order forward difference

    % Calcul du flux de chaleur à x=0 en utilisant la loi de Fourier
    Ac = (pi() * D^2) / 4;   % Aire de la section transversale de l'ailette
    q_numerique = -k_cuivre * Ac * gradT_x0_num;  % Flux de chaleur à x=0
    
    % Calcul du flux de chaleur analytique
    syms x m L
    T_expr = T_inf + (Tm - T_inf) * cosh(m * (L - x)) / cosh(m * L);
    dTdx = diff(T_expr, x);

    % Calcul de m et évaluation de la dérivée symbolique en x=0
    m_val = (h * P) / (k_cuivre * Ac);
    gradT_x0_ana = double(subs(dTdx, [x, m, L], [0, m_val, Longueur]));
    q_analytique = -k_cuivre * Ac * gradT_x0_ana;

end
