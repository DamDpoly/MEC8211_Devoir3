function [T_numerique, T_analytique, Q_numerique, Q_analytique] = Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot)
    % Paramètres géométriques et physiques
    P = pi * D;                     
    Ac = (pi * D^2) / 4;           
    m2 = (h * P) / (k * Ac);        % m² : paramètre de l'équation différentielle

    dx = L / (Ntot - 1);
    x = linspace(0, L, Ntot);

    % Initialisation des matrices A et B
    A = zeros(Ntot);
    B = zeros(Ntot, 1);

    % --- Condition de Dirichlet à x = 0 (température fixée) ---
    A(1, 1) = 1;
    B(1) = Tm;

    % --- Nœuds internes ---
    for i = 2:Ntot-1
        A(i, i-1) = 1;
        A(i, i)   = -2 - m2 * dx^2;
        A(i, i+1) = 1;
        B(i) = -m2 * dx^2 * T_inf;
    end

    % --- Condition de convection à l'extrémité x = L (schéma d'ordre 1) ---
    A(Ntot, Ntot-1) = -1/dx;
    A(Ntot, Ntot)   = 1/dx + h/k;
    B(Ntot) = (h/k) * T_inf;

    % Résolution du système linéaire
    T_numerique = A \ B;

    % --- Solution analytique ---
    m = sqrt(m2);
    T_analytique = T_inf + (Tm - T_inf) * ...
        (((h / (k * m)) * sinh(m * (L - x)) + cosh(m * (L - x))) / ...
         ((h / (k * m)) * sinh(m * L) + cosh(m * L)));

    % --- Puissance thermique totale à x = 0 (numérique) ---
    gradT_x0_num = (T_numerique(2) - T_numerique(1)) / dx;
    Q_numerique = -k * Ac * gradT_x0_num;

    % --- Puissance thermique totale à x = 0 (analytique) ---
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

