%{
Ce code est une fonction permettant la résolution de l'équation de diffusion 
du sel suivant un modèle stationnaire. La discrétisation spatiale du pilier 
est gérable avec le paramètre Ntot, qualifiant le nombre total de noeuds du 
système. Les paramètres liés aux dimensions du pilier et à la concentration
du sel sont constants. 
%}

% Définition des constantes de l'exercice
D = 0.05; %(m)
Longueur = 1; %(m)
k_cuivre = 385; %(W/m.K)
h = 20; %(W/m^2.K)
P = pi()* D; %(m)
Ac = ( pi() * D^2)/4;%(m^2)
m = (h * P)/(k_cuivre * Ac); %Coefficient

T_inf = 25;
Tm = 100;
Ntot = 5;

% Calcul du pas de discrétisation dx
dx = Longueur / (Ntot - 1);

% Définition des matrices A et B pour résoudre A*C = B
A = zeros(Ntot-2, Ntot-2); 
B = zeros(Ntot-2, 1);       

% Remplissage de la matrice A
% Condition limite T0 = Tm
A(1, 1) = (-m^2 * dx^2) - 2; 
A(1, 2) = 1;

% Condition limite T4 = T3
A(Ntot-2, Ntot-2-1) = 1;
A(Ntot-2, Ntot-2) = (-m^2 * dx^2) - 2 + 1;

% Boucle de remplissage pour les points internes
for i = 2:Ntot-3
    A(i, i-1) = 1;
    A(i, i) = (-m^2 * dx^2) - 2;
    A(i, i+1) = 1;
end

% Remplissage de la matrice B
B(1) = (-m^2 * dx^2 * T_inf) - Tm;
for i = 2:Ntot-2
    B(i) = -m^2 * dx^2 * T_inf;
end

% Résolution du système
T = A \ B;

% Reconstitution du vecteur température complet
T_numerique = [Tm; T; T(Ntot-2)];

% Positions des nœuds
x = linspace(0, Longueur, Ntot);

% Solution analytique
T_analytique = T_inf + (Tm - T_inf) * cosh(m * (Longueur - x)) / cosh(m * Longueur);

% AFFICHAGE
plot(x, T_numerique, 'ro-', 'LineWidth', 1.5)
hold on
plot(x, T_analytique, 'b--', 'LineWidth', 1.5)
xlabel('Position x (m)')
ylabel('Température (°C)')
legend('Numérique', 'Analytique')
title('Comparaison des solutions numérique et analytique')
grid on

