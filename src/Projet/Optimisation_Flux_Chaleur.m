% --- Définition des constantes ---
D = 0.05;           % Diamètre du pilier (m)
Longueur = 1;       % Longueur du pilier (m)
k_cuivre = 385;     % Conductivité thermique (W/m.K)
h = 20;             % Coefficient d'échange thermique (W/m².K)
T_inf = 25;         % Température ambiante (°C)
Tm = 100;           % Température à la base (°C)

% Call the function to get temperature profiles
[T_numerique, T_analytique] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

% Compute dx and Ac
dx = Longueur / (Ntot - 1);
Ac = (pi * D^2) / 4;

% Numerical heat flux using second-order finite difference
T0 = T_numerique(1);
T1 = T_numerique(2);
T2 = T_numerique(3);
gradT_x0_num = (-3*T0 + 4*T1 - T2) / (2*dx);
q_num = -k_cuivre * Ac * gradT_x0_num;

% Symbolic computation of analytical heat flux
syms x m L
T_expr = T_inf + (Tm - T_inf) * cosh(m * (L - x)) / cosh(m * L);
dTdx = diff(T_expr, x);

% Compute m and evaluate symbolic derivative at x=0
P = pi * D;
m_val = (h * P) / (k_cuivre * Ac);
gradT_x0_ana = double(subs(dTdx, [x, m, L], [0, m_val, Longueur]));
q_ana = -k_cuivre * Ac * gradT_x0_ana;

% Display results
fprintf('Flux de chaleur NUMÉRIQUE à x = 0 : %.3f W\n', q_num);
fprintf('Flux de chaleur ANALYTIQUE à x = 0 : %.3f W\n', q_ana);

% Relative error
erreur_relative = abs(q_num - q_ana) / abs(q_ana) * 100;
fprintf('Erreur relative : %.4f %%\n', erreur_relative);
