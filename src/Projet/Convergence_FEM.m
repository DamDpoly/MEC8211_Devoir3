% Données d'entrée
D = 0.05;             % Diamètre du pilier (m)
L = 0.4;              % Longueur du pilier (m)
k = 150;              % Conductivité thermique (W/m.K)
h = 15;               % Coefficient d'échange thermique (W/m².K)
T_inf = 300;           % Température ambiante (°C)
Tm = 400;             % Température à la base (°C)
Hmin = 0.00020;
Hmax = 0.00100;
dH   = 0.00020;          % Taille du maillage
H_values = Hmin:dH:Hmax;          % Taille du maillage
num_z = 50;
Ntot = 50;

% Appeler la fonction FEM avec plot_flag défini sur vrai
[T_FEM, q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, Hmin, num_z, true);

disp(q_FEM)  % Afficher le flux de chaleur calculé

[L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur_FEM(D, L, k, h, T_inf, Tm, Ntot, H_values, num_z);
[~, ~, ~, q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);
disp(q_analytique)

Erreurs_T = [H_values', L1_error_T, L2_error_T, Linf_error_T];
Erreurs_q = [H_values', Erreur_q];

% --- Affichage graphique pour la température et le flux de chaleur ---
% Graphique des erreurs pour la température
figure;
loglog(H_values, L1_error_T, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Temp)');
hold on;
loglog(H_values, L2_error_T, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2 (Temp)');
loglog(H_values, Linf_error_T, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf (Temp)');
xlabel('Maillage H');
ylabel('Erreur');
title('Erreurs (L1, L2, Linf) pour la température');
legend show;
grid on;

% Graphique de l'erreur pour le flux de chaleur
figure;
loglog(H_values, Erreur_q, 'o--', 'LineWidth', 2, 'DisplayName', 'Norme L1 (Flux)');

xlabel('Maillage H');
ylabel('Erreur');
title('Erreur pour le flux de chaleur');
legend show;
grid on;