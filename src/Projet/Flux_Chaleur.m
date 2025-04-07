function q_x0 = Flux_Chaleur(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot)
    % Call the function to get the numerical temperature profile
    [T_numerique, ~] = Solution_numerique_ailette(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);

    % Discretisation (spacing between nodes)
    dx = Longueur / (Ntot - 1);

    % Apply second-order central difference for the temperature gradient at x=0
    % Second-order central difference: (T_2 - T_(-1)) / (2 * dx)
    T0 = T_numerique(1);  % Temperature at x=0 (T_0)
    T1 = T_numerique(2);        % Temperature at x=dx (T_1)
    T2 = T_numerique(3);        % Temperature at x=2*dx (T_2)
    
    gradT_x0 = (-3*T0 + 4*T1 - T2) / (2*dx);  % Central difference approximation

    % Calculate the heat flux at x=0 using Fourier's law
    Ac = (pi() * D^2) / 4;   % Cross-sectional area of the fin
    q_x0 = -k_cuivre * Ac * gradT_x0;  % Heat flux at x=0
    
    % Display the result
    fprintf('Flux de chaleur à x = 0 : %.3f W\n', q_x0);
end
