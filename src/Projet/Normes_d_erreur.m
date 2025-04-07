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
espace_values = Longueur ./ (Ntot_values - 1);  % Just used for plotting

L1_error = zeros(length(Ntot_values), 1);
L2_error = zeros(length(Ntot_values), 1);
Linf_error = zeros(length(Ntot_values), 1);

for idx = 1:length(Ntot_values)
    Ntot = Ntot_values(idx);

    [T_numerique, T_analytique] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

    erreur = abs(T_numerique - T_analytique');

    L1_error(idx) = sum(erreur);
    L2_error(idx) = sqrt(sum(erreur.^2));
    Linf_error(idx) = max(erreur);
end

% --- Compilation des erreurs en colonnes ---
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


