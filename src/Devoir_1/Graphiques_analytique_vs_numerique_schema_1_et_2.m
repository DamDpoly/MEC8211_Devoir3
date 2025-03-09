%{
Ce code appelle les fonctions "Solution_numerique_schema_1" et "Solution_numerique_schema_2" et créer quatre
graphiques pour comparer les solutions numérique obtenues avec les schémas d'ordre 1 et 2 et la solution
analytique à diverses valeurs de Ntot (la quantité totale de noeuds) 
%}

clc;
clear;

% Valeurs sélectionnées pour Ntot
Ntot_valeurs = [5, 10, 100, 1000];

figure;

for i = 1:length(Ntot_valeurs)
    Ntot = Ntot_valeurs(i);
    
    % Appel de la fonction Calculate_Num pour obtenir C, C_analytique et r
    % pour le schéma d'ordre 1 et 2
    [C_1, C_analytique_1, r_1, dr_1] = Solution_numerique_schema_1(Ntot);
    [C_2, C_analytique_2, r_2, dr_2] = Solution_numerique_schema_2(Ntot);

    % Tracé de la comparaison des concentrations analytiques et numériques
    subplot(2, 2, i);
    hold off;
    plot(r_1, C_1, '-o', 'DisplayName', 'Résultat Numérique Ordre 1', 'Color', 'b', 'LineWidth', 1);
    hold on;
    plot(r_1, C_analytique_1, '-x', 'DisplayName', 'Résultat Analytique', 'Color',[1, 0.2, 0.2],'LineWidth', 2);
    hold on;
    plot(r_2, C_2, '-s', 'DisplayName', 'Résultat Numérique Ordre 2', 'Color', [0, 0.5, 0], 'LineWidth', 1);
    
    xlabel('Rayon r (m)');
    ylabel('Concentration C (mol/m^3)');
    title(['Comparaison entre les résultats numériques et analytiques (N = ' num2str(Ntot) ')']);
    legend('show','Location','northwest');
    grid on;
end
