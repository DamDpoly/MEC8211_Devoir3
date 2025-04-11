function [T_numerique, T_analytique, q_numerique, q_analytique] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot)
    % Geometry and physical parameters
    P = pi * D;                       
    Ac = (pi * D^2) / 4;             
    m = sqrt((h * P) / (k * Ac));     

    dx = L / (Ntot - 1);
    x = linspace(0, L, Ntot);

    % Initialize matrix A and RHS B
    A = zeros(Ntot);
    B = zeros(Ntot, 1);

    % --- Dirichlet BC at x = 0 ---
    A(1, 1) = 1;
    B(1) = Tm;

    % --- Interior nodes ---
    for i = 2:Ntot-1
        A(i, i-1) = 1;
        A(i, i)   = -2 - m^2 * dx^2;
        A(i, i+1) = 1;
        B(i) = -m^2 * dx^2 * T_inf;
    end

    % --- Neumann BC at x = L ---
    % dT/dx = -h/k * (T - T_inf)
    % Use backward difference: (3T_N - 4T_{N-1} + T_{N-2}) / (2dx) = -h/k (T_N - T_inf)

    A(Ntot, Ntot-2) = 1 / (2*dx);
    A(Ntot, Ntot-1) = -4 / (2*dx);
    A(Ntot, Ntot)   = 3 / (2*dx) + h / k;
    B(Ntot) = (h / k) * T_inf;

    % Solve the linear system
    T_numerique = A \ B;

    % --- Analytical Solution ---
    T_analytique = T_inf + (Tm - T_inf) * ...
        (((h / (k * m)) * sinh(m * (L - x)) + cosh(m * (L - x))) / ...
         ((h / (k * m)) * sinh(m * L) + cosh(m * L)));

    % --- Heat Flux at x = 0 ---
    T0 = T_numerique(1);
    T1 = T_numerique(2);
    gradT_x0_num = (T1 - T0) / dx;
    q_numerique = -k * Ac * gradT_x0_num;

    % --- Analytical heat flux (symbolic) ---
    syms xs ms ys T_infs Tms hs ks real

    % Symbolic expression for temperature
    T_symb_expr = T_infs + (Tms - T_infs) * ...
        (((hs / (ks * ms)) * sinh(ms * (ys - xs)) + cosh(ms * (ys - xs))) / ...
         ((hs / (ks * ms)) * sinh(ms * ys) + cosh(ms * ys)));

    % Derivative of T w.r.t. x
    dTdx = diff(T_symb_expr, xs);

    % Evaluate at x = 0
    gradT_x0_ana = double(subs(dTdx, ...
        [xs,    ms,  ys,   T_infs, Tms,  hs,  ks], ...
        [0,     m,   L,    T_inf,  Tm,   h,   k]));

    % Fourier's law to compute analytical heat flux at x = 0
    q_analytique = -k * Ac * gradT_x0_ana;
    
end
