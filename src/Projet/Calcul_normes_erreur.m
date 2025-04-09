function [L1_error_T, L2_error_T, Linf_error_T, Erreur_q] = Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values)
    % --- Erreur en fonction de la discrétisation spatiale ---
    dx_values = L ./ (Ntot_values - 1);

    % Initializing error arrays for temperature (T) and heat flux (q)
    L1_error_T = zeros(length(Ntot_values), 1);
    L2_error_T = zeros(length(Ntot_values), 1);
    Linf_error_T = zeros(length(Ntot_values), 1);
    
    Erreur_q = zeros(length(Ntot_values), 1);

    for idx = 1:length(Ntot_values)
        Ntot = Ntot_values(idx);
        dx = dx_values(idx);

        % Appel de la fonction pour obtenir les profils de température et les flux de chaleur
        [T_numerique, T_analytique, q_num, q_ana] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, Ntot);

        % --- Calcul de l'erreur pour les profils de température ---
        erreur_T = abs(T_numerique - T_analytique');  % Difference between numerical and analytical temperature

        % Calcul des différentes normes d'erreur pour la température
        L1_error_T(idx) = sum(erreur_T) * dx;
        L2_error_T(idx) = sqrt(sum(erreur_T.^2) * dx);
        Linf_error_T(idx) = max(erreur_T);

        % --- Calcul de l'erreur pour les flux de chaleur ---
        Erreur_q(idx) = abs(q_num - q_ana);  % Difference between numerical and analytical heat flux

    end
end
