clear;

% --- Définition des constantes ---
D = 0.05;           % Diamètre du pilier (m)
Longueur = 1;       % Longueur du pilier (m)
k_cuivre = 385;     % Conductivité thermique (W/m.K)
h = 20;             % Coefficient d'échange thermique (W/m².K)
T_inf = 25;         % Température ambiante (°C)
Tm = 100;           % Température à la base (°C)

%
Ntot_values = [5, 10, 15, 20];
dx_values = Longueur ./ (Ntot_values - 1);

% Préparation de la figure pour les sous-graphes
figure;

for i = 1:length(Ntot_values)
    Ntot = Ntot_values(i);
    dx = dx_values(i);

    % Appel de la fonction pour obtenir les profils de température et les flux de chaleur
    [T_numerique, T_analytique, q_num, q_ana] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

    % Affichage des profils de température
    subplot(2, 2, i);
    
    % Affichage des profils de température numérique vs analytique
    plot(linspace(0, Longueur, Ntot), T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(linspace(0, Longueur, Ntot), T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    
    % Affichage du flux de chaleur correspondant à x=0
    % Calcul du flux de chaleur numérique et analytique à x=0
    Ac = (pi * D^2) / 4; % Aire de la section transversale
    dx = Longueur / (Ntot - 1);
    T0 = T_numerique(1);
    T1 = T_numerique(2);
    T2 = T_numerique(3);
    gradT_x0_num = (-3*T0 + 4*T1 - T2) / (2*dx);  % Méthode des différences centrales pour le gradient à x=0
    q_num_x0 = -k_cuivre * Ac * gradT_x0_num;
    
    % Calcul du flux de chaleur analytique
    syms x m L
    T_expr = T_inf + (Tm - T_inf) * cosh(m * (L - x)) / cosh(m * L);
    dTdx = diff(T_expr, x);
    m_val = (h * pi * D) / (k_cuivre * Ac);
    gradT_x0_ana = double(subs(dTdx, [x, m, L], [0, m_val, Longueur]));
    q_ana_x0 = -k_cuivre * Ac * gradT_x0_ana;
    
    % Affichage du flux numérique et analytique à x=0 dans le titre
    title(sprintf('Ntot = %d\nq_{num} = %.3f W, q_{ana} = %.3f W', Ntot, q_num_x0, q_ana_x0));
    
    % Ajouter les labels et la légende
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end

% --- Affichage général ---
sgtitle('Profils de température et flux de chaleur pour différentes valeurs de Ntot');

% --- Calcul des erreurs pour les profils de température et les flux de chaleur ---
Ntot_values = 5:1:500;

dx_values = Longueur ./ (Ntot_values - 1);

[L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot_values);

% --- Compilation des erreurs en colonnes (exportable vers Excel) ---
Erreurs_T = [dx_values', L1_error_T, L2_error_T, Linf_error_T];
Erreurs_q = [dx_values', Erreur_q];

% --- Affichage graphique pour la température et le flux de chaleur ---
% Graphique des erreurs pour la température
figure;
loglog(dx_values, L1_error_T, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Temp)');
hold on;
loglog(dx_values, L2_error_T, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2 (Temp)');
loglog(dx_values, Linf_error_T, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf (Temp)');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreurs (L1, L2, Linf) pour la température');
legend show;
grid on;

% Graphique de l'erreur pour le flux de chaleur
figure;
loglog(dx_values, Erreur_q, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Flux)');

xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreur pour le flux de chaleur');
legend show;
grid on;
