function [L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur_FEM(D, L, k, h, T_inf, Tm, Ntot, H_values, num_z)

    % Initializing error arrays for temperature (T) and heat flux (q)
    L1_error_T = zeros(length(H_values), 1);
    L2_error_T = zeros(length(H_values), 1);
    Linf_error_T = zeros(length(H_values), 1);
    
    Erreur_q = zeros(length(H_values), 1);

    for idx = 1:length(H_values)
        H = H_values(idx);

        % Appel de la fonction pour obtenir les profils de température et les flux de chaleur
        [T_FEM, q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, H, num_z, false);
        [~, T_analytique, ~, q_ana] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);

        % --- Calcul de l'erreur pour les profils de température ---
        erreur_T = abs(T_FEM - T_analytique');  % Difference between numerical and analytical temperature

        % Calcul des différentes normes d'erreur pour la température
        L1_error_T(idx) = sum(erreur_T) /Ntot;
        L2_error_T(idx) = sqrt(sum(erreur_T.^2) /Ntot);
        Linf_error_T(idx) = max(erreur_T);

        % --- Calcul de l'erreur pour les flux de chaleur ---
        Erreur_q(idx) = abs(q_FEM - q_ana);  % Difference between numerical and analytical heat flux

    end
end