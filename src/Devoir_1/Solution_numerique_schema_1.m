%{
Ce code est une fonction permettant la résolution de l'équation de diffusion 
du sel suivant un modèle stationnaire. La discrétisation spatiale du pilier 
est gérable avec le paramètre Ntot, qualifiant le nombre total de noeuds du 
système. Les paramètres liés aux dimensions du pilier et à la concentration
du sel sont constants. 
%}

function [C, C_analytique, r, dr] = Solution_numerique_schema_1(Ntot)

    % Définition des constantes de l'exercice
    Deff = 10^(-10);  
    S = 2*10^(-8);
    Ce = 20; 
    D = 1; 
    Rayon = D /2;
    alpha = S / Deff;
    
    % Calcul du pas de discrétisation dr
    dr = Rayon / (Ntot - 1);
    
    % Définition des matrices A et B pour pouvoir résoudre A*C = B
    A = zeros(Ntot-2, Ntot-2); 
    B = zeros(Ntot-2, 1);       
    
    % Remplissage de la matrice A
    
    %Condition limite C1 = C0
    A(1, 1) = 1 - (2*1 + 1); 
    A(1, 1+1) = 1 + 1;
    
    %Condition limite CN-1 = CE
    A(Ntot-2, Ntot-2-1) = Ntot-2;
    A(Ntot-2, Ntot-2) = -(2*(Ntot-2) + 1);
    
    % Boucle de remplissage de A
    for i = 2:Ntot-3
        A(i, i-1) = i;
        A(i, i) = -(2*i + 1);
        A(i, i+1) = i + 1;
    end
    
    % Remplissage de la matrice B
    B(Ntot-2, 1) = alpha * dr^2 * (Ntot-2) - (Ntot-1) * Ce;
    for i = 1:Ntot-3
        B(i, 1) = alpha * dr^2 * i;  % Cas global 
    end
    
    % Résolution du système
    C = A \ B; 
    
    % Concaténation des résultats calculés avec ceux des conditions limites
    C = [C(1); C; Ce];
    
    % Définition de la position de chaque noeud dans le vecteur r
    r = linspace(0, Rayon, Ntot);
    
    % Calcul de la concentration avec la solution analytique
    C_analytique = (1/4) * alpha * Rayon^2 * ((r.^2 / Rayon^2) - 1) + Ce;
    
end
