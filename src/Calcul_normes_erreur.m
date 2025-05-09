%Ce code est une fonction calculant les normes d'erreur pour la simulation
%FDS et la simulation FEM. Les flags sont utilisés pour enlever les calculs
%de l'un ou de l'autre afin de réduire la quantité de calculs tout en
%mettant les calculs pour FDS et FEM dans la même fonction.

function [L1_error_T, L2_error_T, Linf_error_T, Erreur_q, ...
          L1_error_T_FEM, L2_error_T_FEM, Linf_error_T_FEM, Erreur_q_FEM] = ...
          Calcul_normes_erreur(D, L, k, h, T_inf, Tm, Ntot_values, H_values, num_z, flag_FDS, flag_FEM)

    % Initialisation
    L1_error_T = [];
    L2_error_T = [];
    Linf_error_T = [];
    Erreur_q = [];

    L1_error_T_FEM = [];
    L2_error_T_FEM = [];
    Linf_error_T_FEM = [];
    Erreur_q_FEM = [];

    % --- ERREURS POUR LE MODELE FDS --- %
    if flag_FDS
        dx_values = L ./ (Ntot_values - 1);

        L1_error_T = zeros(length(Ntot_values), 1);
        L2_error_T = zeros(length(Ntot_values), 1);
        Linf_error_T = zeros(length(Ntot_values), 1);
        Erreur_q = zeros(length(Ntot_values), 1);

        for idx = 1:length(Ntot_values)
            Ntot = Ntot_values(idx);
            dx = dx_values(idx);

            [T_numerique, T_analytique, q_num, q_ana] = Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot);

            erreur_T = abs(T_numerique - T_analytique');

            L1_error_T(idx) = sum(erreur_T) * dx;
            L2_error_T(idx) = sqrt(sum(erreur_T.^2) * dx);
            Linf_error_T(idx) = max(erreur_T);

            Erreur_q(idx) = abs(q_num - q_ana);
        end
    end

    % --- ERREURS POUR LE MODELE FEM --- %
    if flag_FEM
        L1_error_T_FEM = zeros(length(H_values), 1);
        L2_error_T_FEM = zeros(length(H_values), 1);
        Linf_error_T_FEM = zeros(length(H_values), 1);
        Erreur_q_FEM = zeros(length(H_values), 1);

        for idx = 1:length(H_values)
            H = H_values(idx);

            [T_FEM, q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, H, num_z, false);

            % Use consistent analytic reference (from schema 1)
            [~, T_analytique_FEM, ~, q_ana_FEM] = Solution_numerique_ailette_schema_1(D, L, k, h, T_inf, Tm, Ntot_values(end));

            erreur_T_FEM = abs(T_FEM - T_analytique_FEM');

            L1_error_T_FEM(idx) = sum(erreur_T_FEM) / num_z;
            L2_error_T_FEM(idx) = sqrt(sum(erreur_T_FEM.^2) / num_z);
            Linf_error_T_FEM(idx) = max(erreur_T_FEM);

            Erreur_q_FEM(idx) = abs(q_FEM - q_ana_FEM);
        end
    end
end
