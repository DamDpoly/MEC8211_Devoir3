clear;

% --- Définition des constantes ---
D = 0.5;             % Diamètre du pilier (m)
L = 1;               % Longueur du pilier (m)
k = 385;             % Conductivité thermique (W/m.K)
h = 500;             % Coefficient d'échange thermique (W/m².K)
T_inf = 25;          % Température ambiante (°C)
Tm = 100;            % Température à la base (°C)


Bi = (h * D) / k;    % Nombre de Biot


% Valeurs de Ntot pour les quatres profils de concentration
Ntot_values = [5, 10, 15, 1000];
dx_values = L ./ (Ntot_values - 1);

figure;

for i = 1:length(Ntot_values)
    Ntot = Ntot_values(i);
    dx = dx_values(i);

    % Appel de la fonction pour obtenir les profils de température et les flux de chaleur
    [T_numerique, T_analytique, q_numerique, q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);

    % Affichage des profils de température
    subplot(2, 2, i);
    
    % Affichage des profils de température numérique vs analytique
    plot(linspace(0, L, Ntot), T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(linspace(0, L, Ntot), T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    
    % Affichage du flux numérique et analytique à x=0 dans le titre
    title(sprintf('Ntot = %d, Bi = %.4f\nq_{num} = %.3f W, q_{ana} = %.3f W', Ntot, Bi, q_numerique, q_analytique));
    
    % Ajouter les labels et la légende
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end

% --- Affichage général ---
sgtitle('Profils de température et flux de chaleur pour différentes valeurs de Ntot');

% --- Calcul des erreurs pour les profils de température et les flux de chaleur ---
Ntot_values = 5:1:200;

dx_values = L ./ (Ntot_values - 1);

[L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values);

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

