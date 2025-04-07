% --- Définition des constantes ---
D = 0.05;           % Diamètre du pilier (m)
Longueur = 1;       % Longueur du pilier (m)
k_cuivre = 385;     % Conductivité thermique (W/m.K)
h = 20;             % Coefficient d'échange thermique (W/m².K)
T_inf = 25;         % Température ambiante (°C)
Tm = 100;           % Température à la base (°C)
Ntot = 5000;        % Nombre total de points

% Appel de la fonction pour obtenir les profils de température
[T_numerique, T_analytique] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

% Calcul de dx et Ac
dx = Longueur / (Ntot - 1);
Ac = (pi * D^2) / 4;

% Flux de chaleur numérique en utilisant la méthode des différences finies d'ordre 2
T0 = T_numerique(1);
T1 = T_numerique(2);
T2 = T_numerique(3);
gradT_x0_num = (-3*T0 + 4*T1 - T2) / (2*dx);
q_num = -k_cuivre * Ac * gradT_x0_num;

% Calcul symbolique du flux de chaleur analytique
syms x m L
T_expr = T_inf + (Tm - T_inf) * cosh(m * (L - x)) / cosh(m * L);
dTdx = diff(T_expr, x);

% Calcul de m et évaluation de la dérivée symbolique en x=0
P = pi * D;
m_val = (h * P) / (k_cuivre * Ac);
gradT_x0_ana = double(subs(dTdx, [x, m, L], [0, m_val, Longueur]));
q_ana = -k_cuivre * Ac * gradT_x0_ana;

% Affichage des résultats
fprintf('Flux de chaleur numérique à x = 0 : %.3f W\n', q_num);
fprintf('Flux de chaleur analytique à x = 0 : %.3f W\n', q_ana);

