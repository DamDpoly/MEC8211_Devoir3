

function [T_numerique, T_analytique] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot)
   
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

    % AFFICHAGE
    figure;
    plot(x, T_numerique, 'ro-', 'LineWidth', 1.5)
    hold on
    plot(x, T_analytique, 'b--', 'LineWidth', 1.5)
    xlabel('Position x (m)')
    ylabel('Température (°C)')
    legend('Numérique', 'Analytique')
    title('Comparaison des solutions numérique et analytique')
    grid on
end
