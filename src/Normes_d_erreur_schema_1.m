%{
Ce code appelle la fonction "Solution_numerique_schema_1" et calcule les normes
d'erreur L1, L2 et L∞ pour chaque valeur de Ntot (la quantité totale de
noeuds). Ntot est compris entre 5 et Nmax. Les erreurs sont tracées
graphiquement sur une échelle logarithmique. 
%}

clc;
clear;

Nmax = 100;

% Initialisation des vecteurs
L1_error = zeros(Nmax, 1);
L2_error = zeros(Nmax, 1);
Linf_error = zeros(Nmax, 1);
dr = zeros(Nmax, 1);

% Calcul des erreurs pour chaque valeur de Ntot
for Ntot = 5:Nmax

    % Appel de la fonction Solution_numerique_schema_1 pour obtenir C et dr
    [C, C_analytique, r, dr(Ntot)] = Solution_numerique_schema_1(Ntot);
    
    % Calcul des normes d'erreur
    erreur = abs(C - C_analytique');
	L1_error(Ntot) = sum(erreur) * dr(Ntot);
	L2_error(Ntot) = sqrt(sum(erreur.^2) * dr(Ntot));
	Linf_error(Ntot) = max(erreur);
end

% Tracé des erreurs en fonction du pas de discrétization
figure;
semilogy(dr(5:Nmax), L1_error(5:Nmax), 'r', 'DisplayName', 'Norme L1'); hold on;
semilogy(dr(5:Nmax), L2_error(5:Nmax), 'b', 'DisplayName', 'Norme L2');
semilogy(dr(5:Nmax), Linf_error(5:Nmax), 'g', 'DisplayName', 'Norme L∞');
set(gca, 'XScale', 'log');
legend('show');
xlabel('Pas de grille \Delta x');
ylabel('Erreur');
title("Normes d'erreur en fonction de \Delta x");
grid on;