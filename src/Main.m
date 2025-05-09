%Ce code est sous-divisé en trois parties. La première permet la
%visualisation des résultats FDS vs le résultat analytique, la deuxième
%mesure la convergence des deux schémas et de la SRQ pour la FDS, et la
%troisième calcule la convergence du modèle FEM.

clear;
clc

% --- Définition des constantes ---
D = 0.050;             % Diamètre du pilier (m)
L = 0.400;             % Longueur du pilier (m)
k = 150;               % Conductivité thermique (W/m.K)
h = 15;                % Coefficient d'échange thermique (W/m².K)
T_inf = 300;           % Température ambiante (°C)
Tm = 400;              % Température à la base (°C)

Bi = (h * D) / k;      % Nombre de Biot

%% ---------------------- PARTIE 1 : VISUALISATION ---------------------- %%
Ntot_values_plot = [5, 10, 15, 500];
dx_values_plot = L ./ (Ntot_values_plot - 1);

figure;
for i = 1:length(Ntot_values_plot)
    Ntot = Ntot_values_plot(i);
    dx = dx_values_plot(i);
    
    [T_numerique, T_analytique, Q_numerique, Q_analytique] = ...
        Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);

    subplot(2, 2, i);
    plot(linspace(0, L, Ntot), T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(linspace(0, L, Ntot), T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    title(sprintf('Ntot = %d, Bi = %.4f\nQ_{num} = %.3f W, Q_{ana} = %.3f W', ...
        Ntot, Bi, Q_numerique, Q_analytique));
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end
sgtitle('Profils de température et débit de chaleur pour différentes valeurs de Ntot');

figure;
for i = 1:length(Ntot_values_plot)
    Ntot = Ntot_values_plot(i);
    dx = dx_values_plot(i);
    
    [T_numerique, T_analytique, Q_numerique, Q_analytique] = ...
        Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot);

    subplot(2, 2, i);
    x = linspace(0, L, Ntot);
    plot(x, T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(x, T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    title(sprintf('Ntot = %d, Bi = %.4f\nQ_{num} = %.3f W, Q_{ana} = %.3f W', ...
        Ntot, Bi, Q_numerique, Q_analytique));
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end
sgtitle('Profils de température et débit de chaleur (schéma 1) pour différentes valeurs de Ntot');

%% ------------------ PARTIE 2 : Convergence FDS (Ntot = 5:200) ------------------ %%
Ntot_values = 5:1:200;
dx_values = L ./ (Ntot_values - 1);

[L1_error_T, L2_error_T, Linf_error_T, Erreur_Q, ~, ~, ~, ~] = ...
    Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values, 1, 1, true, false);

Q_numerique_all = zeros(length(Ntot_values), 1);
Q_analytique_all = zeros(length(Ntot_values), 1);

for i = 1:length(Ntot_values)
    Ntot = Ntot_values(i);
    [~, ~, Q_num, Q_ana] = Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot);
    Q_numerique_all(i) = Q_num;
    Q_analytique_all(i) = Q_ana;
end

Erreurs_T = [Ntot_values', dx_values', L1_error_T, L2_error_T, Linf_error_T];
Erreurs_Q = [Ntot_values', dx_values', Erreur_Q, Q_numerique_all, Q_analytique_all];

% Graphique erreurs température (FDS)
figure;
loglog(dx_values, L1_error_T, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Temp)');
hold on;
loglog(dx_values, L2_error_T, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2 (Temp)');
loglog(dx_values, Linf_error_T, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf (Temp)');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreurs (L1, L2, Linf) pour la température (FDS)');
legend show;
grid on;

% Graphique erreurs débit de chaleur (FDS)
figure;
loglog(dx_values, Erreur_Q, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Débit)');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreur pour le débit de chaleur (FDS)');
legend show;
grid on;

disp('--- Erreurs Température FDS ---');
disp(array2table(Erreurs_T, 'VariableNames', {'Ntot','dx','L1','L2','Linf'}));

disp('--- Erreurs Débit FDS ---');
disp(array2table(Erreurs_Q, 'VariableNames', {'Ntot','dx','Erreur_Q','Q_num','Q_ana'}));

%% ------------------- PARTIE 3 : Convergence FEM (H = 0.00040:0.00200) ------------------- %%
Hmin = 0.00080;
Hmax = 0.00200;
dH   = 0.00020;
H_values = Hmin:dH:Hmax;
num_z = 50;
Ntot = 50;

[T_FEM, Q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, Hmin, num_z, true);
disp('Q_FEM:');
disp(Q_FEM);

[~, ~, ~, ~, L1_error_T_FEM, L2_error_T_FEM, Linf_error_T_FEM, Erreur_Q_FEM] = ...
    Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot, H_values, num_z, false, true);

[~, ~, ~, Q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);
disp('Q_analytique:');
disp(Q_analytique);

Erreurs_T_FEM = [H_values', L1_error_T_FEM, L2_error_T_FEM, Linf_error_T_FEM];
Erreurs_Q_FEM = [H_values', Erreur_Q_FEM];

% Graphique erreurs température (FEM)
figure;
loglog(H_values, L1_error_T_FEM, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Temp)');
hold on;
loglog(H_values, L2_error_T_FEM, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2 (Temp)');
loglog(H_values, Linf_error_T_FEM, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf (Temp)');
xlabel('Maillage H');
ylabel('Erreur');
title('Erreurs (L1, L2, Linf) pour la température (FEM)');
legend show;
grid on;

% Graphique erreurs débit de chaleur (FEM)
figure;
loglog(H_values, Erreur_Q_FEM, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Débit)');
xlabel('Maillage H');
ylabel('Erreur');
title('Erreur pour le débit de chaleur (FEM)');
legend show;
grid on;
