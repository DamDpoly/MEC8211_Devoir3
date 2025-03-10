clear;

Deff = 10^(-10);  
S = 2*10^(-8);    
Ce = 20;          
D = 1;            
Rayon = D / 2;   
k = 4 * 10^(-9);  
lambda = 1;

dt_values = 1e8:1e8:1e10;
tmax = 1e11;

C_exact = @(r, t) ((r.^2 / Rayon^2) * (Ce + exp(-lambda * t)) - exp(-lambda * t));

Ntot_values = 5:1:100; 
dr_values = Rayon ./ (Ntot_values - 1);

% --- Erreur en fonction de dr ---
L1_error_dr = [];
L2_error_dr = [];
Linf_error_dr = [];

% Calcul de l'erreur pour différents dr
for Ntot = Ntot_values
    [C_full] = Solution_transitoire_numerique(Ntot, dt_values(1), tmax); 

    r = linspace(0, Rayon, Ntot);
    C_analytic = C_exact(r, tmax);

    error = abs(C_full - C_analytic');

    % Calcul des erreurs L1, L2, Linf pour différents dr
    L1 = sum(error) * dr_values(Ntot == Ntot_values);
    L2 = sqrt(sum(error.^2) * dr_values(Ntot == Ntot_values));
    Linf = max(error);

    % Stockage des erreurs
    L1_error_dr = [L1_error_dr, L1];
    L2_error_dr = [L2_error_dr, L2];
    Linf_error_dr = [Linf_error_dr, Linf];
end

% Tracé de l'erreur en fonction de dr
figure;
loglog(dr_values, L1_error_dr, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1');
hold on;
loglog(dr_values, L2_error_dr, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2');
loglog(dr_values, Linf_error_dr, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf');
xlabel('dr (Échelle Logarithmique)');
ylabel('Erreur (Échelle Logarithmique)');
title('Erreurs (L1, L2, Linf) en fonction de dr');
legend show;
grid on;

% --- Erreur en fonction de dt ---
L1_error_dt = [];
L2_error_dt = [];
Linf_error_dt = [];

Ntot = 100; 

% Calcul de l'erreur pour différents dt
for dt = dt_values
    [C_full] = Solution_transitoire_numerique(Ntot, dt, tmax);

    r = linspace(0, Rayon, Ntot);
    C_analytic = C_exact(r, tmax);

    error = abs(C_full - C_analytic');

    % Calcul des erreurs L1, L2, Linf pour différents dt
    L1 = sum(error) * dr_values(Ntot == Ntot_values);
    L2 = sqrt(sum(error.^2) * dr_values(Ntot == Ntot_values));
    Linf = max(error);

    % Stockage des erreurs
    L1_error_dt = [L1_error_dt, L1];
    L2_error_dt = [L2_error_dt, L2];
    Linf_error_dt = [Linf_error_dt, Linf];
end

% Tracé de l'erreur en fonction de dt
figure;
loglog(dt_values, L1_error_dt, 'o-', 'LineWidth', 2, 'DisplayName', 'Norme L1');
hold on;
loglog(dt_values, L2_error_dt, 's-', 'LineWidth', 2, 'DisplayName', 'Norme L2');
loglog(dt_values, Linf_error_dt, '^-', 'LineWidth', 2, 'DisplayName', 'Norme Linf');
xlabel('dt (Échelle Logarithmique)');
ylabel('Erreur (Échelle Logarithmique)');
title('Erreurs (L1, L2, Linf) en fonction de dt');
legend show;
grid on;




