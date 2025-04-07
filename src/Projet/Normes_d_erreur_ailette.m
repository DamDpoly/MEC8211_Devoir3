clear;

% --- Définition des constantes ---
D = 0.05;           % Diamètre du pilier (m)
Longueur = 1;       % Longueur du pilier (m)
k_cuivre = 385;     % Conductivité thermique (W/m.K)
h = 20;             % Coefficient d'échange thermique (W/m².K)
T_inf = 25;         % Température ambiante (°C)
Tm = 100;           % Température à la base (°C)

% --- Erreur en fonction de la discrétisation spatiale ---
Ntot_values = 5:1:100;
espace_values = Longueur ./ (Ntot_values - 1);  % Just used for plotting and integration

L1_error = zeros(length(Ntot_values), 1);
L2_error = zeros(length(Ntot_values), 1);
Linf_error = zeros(length(Ntot_values), 1);

for idx = 1:length(Ntot_values)
    Ntot = Ntot_values(idx);
    pas = espace_values(idx);

    % Appel de la fonction pour obtenir les profils de température et les flux de chaleur
    [T_numerique, T_analytique, q_num, q_ana] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

    % Calcul de l'erreur entre la solution numérique et analytique
    erreur = abs(T_numerique - T_analytique');

    % Calcul des différentes normes d'erreur
    L1_error(idx) = sum(erreur) * pas;
    L2_error(idx) = sqrt(sum(erreur.^2) * pas);
    Linf_error(idx) = max(erreur);
    
    % Affichage des flux de chaleur (optionnel)
    fprintf('Pour Ntot = %d: Flux numérique = %.3f W, Flux analytique = %.3f W\n', Ntot, q_num, q_ana);
end

% --- Compilation des erreurs en colonnes (exportable vers Excel) ---
Erreurs = [L1_error, L2_error, Linf_error];

% --- Affichage graphique ---
figure;
loglog(espace_values, L1_error, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1');
hold on;
loglog(espace_values, L2_error, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2');
loglog(espace_values, Linf_error, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf');
xlabel('Pas de discrétisation');
ylabel('Erreur');
title('Erreurs (L1, L2, Linf) en fonction du maillage');
legend show;
grid on;
