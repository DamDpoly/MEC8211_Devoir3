clear; 
clc; 
close all;

% Définition des constantes de l'exercice
Deff = 10^(-10);  
S = 2*10^(-8);    
Ce = 20;          
D = 1;            
Rayon = D / 2;   
alpha = S / Deff; 
k = 4 * 10^(-9);  
lambda = 1;

% Discrétisation de l'espace
Ntot = 200;            % Nombre total de points n
dr = Rayon / (Ntot - 1); % Pas d'espace (m)

% Discrétisation du temps
t0 = 0; % t départ (s)      
tend = 4000000000; % t final (s)
dt = 20000;  % Pas de temps (s)

% Nouvelle condition initiale avec la MMS 
C = (((1:Ntot-2) * dr / Rayon).^2) * (Ce + 1) - 1;

% Boucle temporelle
for t = t0:dt:tend
    
    % Definition de A et B
    A = zeros(Ntot-2, Ntot-2);
    B = zeros(Ntot-2, 1);
    
    % Boucle de remplissage des matrices A et B
    for i = 1:Ntot-2
        
        % Terme source MMS
        S_i = (Ce + exp(-lambda * t)) * (k * ((i * dr) / Rayon)^2 - Deff * (4 / Rayon^2)) + lambda * exp(-lambda * t) * (1 - ((i * dr) / Rayon)^2 - k / lambda);
        
        if i == 1 % Condition limite pour le cas i = 1 impacté par la condition limite sur C0
            A(i, i) = 1 + Deff * dt * (2 / dr^2) + dt * k;
            A(i, i+1) = -Deff * dt * (1.5 / dr^2);
            B(i) = C(i) + dt * S_i;
        
        elseif i == Ntot-2 % Condition limite au bord
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

    % Affichage de l'évolution
    if mod(t, tend/10) == 0 %Tracer les solutions tous les 10 instants 
        plot(linspace(0, Rayon, Ntot), C_full, 'LineWidth', 2); %Plot des concentrations à l'instant t 
        xlabel('Rayon r (m)');
        ylabel('Concentration C (mol/m³)');
        title(['Évolution de la concentration à t = ' num2str(t) ' s']);
        grid on;
        drawnow; 
    end
end
