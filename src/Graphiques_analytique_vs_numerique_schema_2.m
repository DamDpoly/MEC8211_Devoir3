%{
Ce code appelle la fonction "Solution_numerique_schema_2" et crée quatre
graphiques pour comparer la solution numérique et la solution
analytique à diverses valeurs de Ntot (la quantité totale de
noeuds)
%}

clc;
clear;

% Valeurs sélectionnées pour Ntot
Ntot_valeurs = [5, 10, 100, 1000];

figure;

for i = 1:length(Ntot_valeurs)
    Ntot = Ntot_valeurs(i);
    
     % Appel de la fonction Calculate_Num pour obtenir C, C_analytique et r
    [C, C_analytique, r, dr] = Solution_numerique_schema_2(Ntot);
    
    % Tracé de la comparaison des concentrations analytiques et numériques
    subplot(2, 2, i);
    plot(r, C, '-o', 'DisplayName', 'Résultat Numérique');
    hold on;
    plot(r, C_analytique, '-x', 'DisplayName', 'Résultat Analytique');
    xlabel('Rayon r (m)');
    ylabel('Concentration C (mol/m^3)');
    title(['Comparaison entre les résultats numériques et analytiques (N = ' num2str(Ntot) ')']);
    legend('show');
    grid on;
end
