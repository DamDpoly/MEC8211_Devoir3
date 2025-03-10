function [C_full] = Solution_transitoire_numerique(Ntot, dt, tmax)
    % Définition des constantes de l'exercice
    Deff = 10^(-10);  
    S = 2*10^(-8);    
    Ce = 20;          
    D = 1;            
    Rayon = D / 2;   
    k = 4 * 10^(-9);  
    lambda = 1;

    % Discrétisation de l'espace
    dr = Rayon / (Ntot - 1); % Pas d'espace (m)

    % Discrétisation du temps
    t0 = 0;  
    tend = tmax; 

    % Nouvelle condition initiale avec la MMS 
    C = zeros(Ntot-2, 1); 
    for i = 1 : Ntot-2
        C(i) = ((i*dr)/Rayon)^2*(Ce + 1) - 1;
    end

    % Boucle temporelle
    for t = t0:dt:tend
        % Definition de A et B
        A = zeros(Ntot-2, Ntot-2);
        B = zeros(Ntot-2, 1);
        
        % Boucle de remplissage des matrices A et B
        for i = 1:Ntot-2
            % Terme source MMS
            S_i = (Ce + exp(-lambda * t)) * (k * ((i * dr) / Rayon)^2 - Deff * (4 / Rayon^2)) + lambda * exp(-lambda * t) * (1 - ((i * dr) / Rayon)^2 - k / lambda);
            
            if i == 1 % Condition limite de Neumann
                A(i, i) = 1 + Deff * dt * (2 / dr^2) + dt * k;
                A(i, i+1) = -Deff * dt * (1.5 / dr^2);
                B(i) = C(i) + Deff * dt * ((i - 0.5) / (dr^2 * i)) * C(i) + dt * S_i;
            
            elseif i == Ntot-2 % Condition limite de Dirichlet
                A(i, i-1) = -Deff * dt * ((i - 0.5) / (dr^2 * i));
                A(i, i) = 1 + Deff * dt * (2 * i / (dr^2 * i)) + dt * k;
                B(i) = C(i) + Deff * dt * ((i + 0.5) / (dr^2 * i)) * Ce + dt * S_i;
            else
                % Points internes 
                A(i, i-1) = -Deff * dt * ((i - 0.5) / (dr^2 * i));
                A(i, i) = 1 + Deff * dt * (2 * i / (dr^2 * i)) + dt * k;
                A(i, i+1) = -Deff * dt * ((i + 0.5) / (dr^2 * i));
                B(i) = C(i) + dt * S_i;
            end
        end
        
        % Résolution du système A C = B
        C_n_plus_1 = A \ B;

        % Condition de Neumann sur le flux qui impose C0
        C0 = -(C_n_plus_1(2) - 4*C_n_plus_1(1)) / 3;  

        % Concaténation de C0 et Ce pour avoir toutes les concentrations
        C_full = [C0; C_n_plus_1; Ce];

        % Mise à jour de C pour la prochaine itération
        C = C_n_plus_1;
    end
end
