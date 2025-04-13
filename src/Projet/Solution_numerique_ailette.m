function [T_numerique, T_analytique, Q_numerique, Q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot)
    % Paramètres géométriques et physiques
    P = pi * D;                       
    Ac = (pi * D^2) / 4;             
    m = sqrt((h * P) / (k * Ac));     

    dx = L / (Ntot - 1);
    x = linspace(0, L, Ntot);

    % Initialisation des matrices A et B
    A = zeros(Ntot);
    B = zeros(Ntot, 1);

    % --- Condition de Dirichlet à x = 0 ---
    A(1, 1) = 1;
    B(1) = Tm;

    % --- Nœuds internes ---
    for i = 2:Ntot-1
        A(i, i-1) = 1;
        A(i, i)   = -2 - m^2 * dx^2;
        A(i, i+1) = 1;
        B(i) = -m^2 * dx^2 * T_inf;
    end

    % --- Condition de Neumann à x = L ---
    % dT/dx = -h/k * (T - T_inf)
    A(Ntot, Ntot-2) = 1 / (2*dx);
    A(Ntot, Ntot-1) = -4 / (2*dx);
    A(Ntot, Ntot)   = 3 / (2*dx) + h / k;
    B(Ntot) = (h / k) * T_inf;

    % Résolution du système linéaire
    T_numerique = A \ B;

    % --- Solution analytique ---
    T_analytique = T_inf + (Tm - T_inf) * ...
        (((h / (k * m)) * sinh(m * (L - x)) + cosh(m * (L - x))) / ...
         ((h / (k * m)) * sinh(m * L) + cosh(m * L)));

    % --- Flux de chaleur total à x = 0 (numérique) ---
    T0 = T_numerique(1);
    T1 = T_numerique(2);
    gradT_x0_num = (T1 - T0) / dx;
    Q_numerique = -k * Ac * gradT_x0_num;

    % --- Flux de chaleur total à x = 0 (analytique) ---
    syms xs ms ys T_infs Tms hs ks real

    T_symb_expr = T_infs + (Tms - T_infs) * ...
        (((hs / (ks * ms)) * sinh(ms * (ys - xs)) + cosh(ms * (ys - xs))) / ...
         ((hs / (ks * ms)) * sinh(ms * ys) + cosh(ms * ys)));

    dTdx = diff(T_symb_expr, xs);

    gradT_x0_ana = double(subs(dTdx, ...
        [xs,    ms,  ys,   T_infs, Tms,  hs,  ks], ...
        [0,     m,   L,    T_inf,  Tm,   h,   k]));

    Q_analytique = -k * Ac * gradT_x0_ana;
end

