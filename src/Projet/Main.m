clear;

% --- Définition des constantes ---
D = 0.050;             % Diamètre du pilier (m)
L = 0.400;             % Longueur du pilier (m)
k = 150;               % Conductivité thermique (W/m.K)
h = 15;                % Coefficient d'échange thermique (W/m².K)
T_inf = 300;           % Température ambiante (°C)
Tm = 400;              % Température à la base (°C)

Bi = (h * D) / k;      % Nombre de Biot

% --- Partie 1 : Visualisation pour quelques Ntot ---
Ntot_values_plot = [5, 10, 15, 500];
dx_values_plot = L ./ (Ntot_values_plot - 1);

figure;

for i = 1:length(Ntot_values_plot)
    Ntot = Ntot_values_plot(i);
    dx = dx_values_plot(i);

    % Calcul des profils de température et flux
    [T_numerique, T_analytique, q_numerique, q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);

    subplot(2, 2, i);
    plot(linspace(0, L, Ntot), T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(linspace(0, L, Ntot), T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    title(sprintf('Ntot = %d, Bi = %.4f\nq_{num} = %.3f W, q_{ana} = %.3f W', Ntot, Bi, q_numerique, q_analytique));
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end

sgtitle('Profils de température et flux de chaleur pour différentes valeurs de Ntot');

% --- Partie 2 : Calculs d'erreur pour Ntot = 5 à 200 ---
Ntot_values = 5:1:200;
dx_values = L ./ (Ntot_values - 1);

% Calcul des normes d'erreur
[L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values);

% Calcul des flux q_num et q_ana pour chaque Ntot
q_numerique_all = zeros(length(Ntot_values), 1);
q_analytique_all = zeros(length(Ntot_values), 1);

for i = 1:length(Ntot_values)
    Ntot = Ntot_values(i);
    [~, ~, q_num, q_ana] = Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot);
    q_numerique_all(i) = q_num;
    q_analytique_all(i) = q_ana;
end

% --- Compilation des erreurs en matrices complètes ---
Erreurs_T = [Ntot_values', dx_values', L1_error_T, L2_error_T, Linf_error_T];
Erreurs_q = [Ntot_values', dx_values', Erreur_q, q_numerique_all, q_analytique_all];

% --- Affichage graphique des erreurs de température ---
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

% --- Affichage graphique des erreurs de flux ---
figure;
loglog(dx_values, Erreur_q, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Flux)');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreur pour le flux de chaleur');
legend show;
grid on;

% --- (Facultatif) Affichage des tableaux dans la console ---
disp('--- Erreurs Température (Ntot, dx, L1, L2, Linf) ---');
disp(array2table(Erreurs_T, 'VariableNames', {'Ntot','dx','L1','L2','Linf'}));

disp('--- Erreurs Flux (Ntot, dx, Erreur_q, q_num, q_ana) ---');
disp(array2table(Erreurs_q, 'VariableNames', {'Ntot','dx','Erreur_q','q_num','q_ana'}));

