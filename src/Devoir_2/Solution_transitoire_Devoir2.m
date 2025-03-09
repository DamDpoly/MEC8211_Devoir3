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

% Discrétisation de l'espace
Ntot = 200;            % Nombre total de points n
dr = Rayon / (Ntot - 1); % Pas d'espace (m)

% Discrétisation du temps
t0 = 0; % t départ (s)      
tend = 4000000000; % t final (s)
dt = 20000;  % Pas de temps (s)

% Initialisation des concentrations à 0 ( C(i,0) = 0 ) 

C = zeros(Ntot-2, 1); 

% Boucle temporelle
for t = t0:dt:tend
    
    % Definition de A et B
    A = zeros(Ntot-2, Ntot-2);
    B = zeros(Ntot-2, 1);
    
    % Boucle de remplissage des matrices A et B
    for i = 1:Ntot-2
        if i == 1 % Condition limite pour le cas i = 1 impacté par la condition limite sur C0

            A(i, i) = 1 + Deff * dt * (2 / dr^2) + dt * k;
            A(i, i+1) = -Deff * dt * (1.5 / dr^2);

            % Dans le second membre on a la concentration calculée au tour
            % précédent + la condition limite C0 avec Neumann
            
            B(i) = C(i) + Deff * dt * ((i - 0.5) / (dr^2 * i)) * C(i);
        
        elseif i == Ntot-2 % Condition limite au bord

            A(i, i-1) = -Deff * dt * ((i - 0.5) / (dr^2 * i));
            A(i, i) = 1 + Deff * dt * (2 * i / (dr^2 * i)) + dt * k;

            % Dans le second membre on a la concetration calculée au tour
            % précédent + la condition limite sur le bord Dirichlet

            B(i) = C(i) + Deff * dt * ((i + 0.5) / (dr^2 * i)) * Ce;
        else
            % Points internes 

            A(i, i-1) = -Deff * dt * ((i - 0.5) / (dr^2 * i));
            A(i, i) = 1 + Deff * dt * (2 * i / (dr^2 * i)) + dt * k;
            A(i, i+1) = -Deff * dt * ((i + 0.5) / (dr^2 * i));

            % Dans le second membre on a la concetration calculée au tour
            % précédent

            B(i) = C(i);
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
    if mod(t, tend/10) == 0 %Tracer les solutions tout les 10 instants 
        plot(linspace(0, Rayon, Ntot), C_full, 'LineWidth', 2); %Plot des concentrations à l'instant t 
        xlabel('Rayon r (m)');
        ylabel('Concentration C (mol/m³)');
        title(['Évolution de la concentration à t = ' num2str(t) ' s']);
        grid on;
        drawnow; 
    end
end