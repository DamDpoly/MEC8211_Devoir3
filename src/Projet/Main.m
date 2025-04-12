clear;

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
    
    % Adjusted the call to Solution_numerique_ailette to match expected arguments
    [T_numerique, T_analytique, q_numerique, q_analytique] = ...
        Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);  % No 'false' argument

    subplot(2, 2, i);
    plot(linspace(0, L, Ntot), T_numerique, '-o', 'DisplayName', 'Température Numérique', 'LineWidth', 2);
    hold on;
    plot(linspace(0, L, Ntot), T_analytique, '-s', 'DisplayName', 'Température Analytique', 'LineWidth', 2);
    title(sprintf('Ntot = %d, Bi = %.4f\nq_{num} = %.3f W, q_{ana} = %.3f W', ...
        Ntot, Bi, q_numerique, q_analytique));
    xlabel('Longueur (m)');
    ylabel('Température (°C)');
    legend('show');
    grid on;
end
sgtitle('Profils de température et flux de chaleur pour différentes valeurs de Ntot');

%% ------------------ PARTIE 2 : Convergence FDS (Ntot = 5:200) ------------------ %%
Ntot_values = 5:1:200;
dx_values = L ./ (Ntot_values - 1);

% Call the error calculation function for FDS
[L1_error_T, L2_error_T, Linf_error_T, Erreur_q, ~, ~, ~, ~] = ...
    Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values, 1, 1, true, false);

q_numerique_all = zeros(length(Ntot_values), 1);
q_analytique_all = zeros(length(Ntot_values), 1);

for i = 1:length(Ntot_values)
    Ntot = Ntot_values(i);
    [~, ~, q_num, q_ana] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);
    q_numerique_all(i) = q_num;
    q_analytique_all(i) = q_ana;
end

Erreurs_T = [Ntot_values', dx_values', L1_error_T, L2_error_T, Linf_error_T];
Erreurs_q = [Ntot_values', dx_values', Erreur_q, q_numerique_all, q_analytique_all];

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

% Graphique erreurs flux (FDS)
figure;
loglog(dx_values, Erreur_q, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Flux)');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreur pour le flux de chaleur (FDS)');
legend show;
grid on;

disp('--- Erreurs Température FDS ---');
disp(array2table(Erreurs_T, 'VariableNames', {'Ntot','dx','L1','L2','Linf'}));

disp('--- Erreurs Flux FDS ---');
disp(array2table(Erreurs_q, 'VariableNames', {'Ntot','dx','Erreur_q','q_num','q_ana'}));

%% ------------------- PARTIE 3 : Convergence FEM (H = 0.00040:0.00200) ------------------- %%
Hmin = 0.00040;
Hmax = 0.00200;
dH   = 0.00020;
H_values = Hmin:dH:Hmax;
num_z = 50;
Ntot = 50;

[T_FEM, q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, Hmin, num_z, true);
disp('q_FEM:');
disp(q_FEM);

% Call the error calculation function for FEM
[~, ~, ~, ~, L1_error_T_FEM, L2_error_T_FEM, Linf_error_T_FEM, Erreur_q_FEM] = ...
    Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot, H_values, num_z, false, true);

[~, ~, ~, q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);
disp('q_analytique:');
disp(q_analytique);

Erreurs_T_FEM = [H_values', L1_error_T_FEM, L2_error_T_FEM, Linf_error_T_FEM];
Erreurs_q_FEM = [H_values', Erreur_q_FEM];

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

% Graphique erreurs flux (FEM)
figure;
loglog(H_values, Erreur_q_FEM, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Flux)');
xlabel('Maillage H');
ylabel('Erreur');
title('Erreur pour le flux de chaleur (FEM)');
legend show;
grid on;
