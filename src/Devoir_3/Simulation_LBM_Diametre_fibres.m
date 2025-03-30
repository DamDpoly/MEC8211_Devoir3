% Ce code, basé sur le code donné launch_simulationLBM, permet d'itérer la
% simulation sur un grand nombre d'échantillons de diamètre de fibres et de
% SEED afin de calculer l'incertitude du paramètre d'entrée. Ce code génère
% des PDF et CDF du paramètre d'entrée et de k ainsi qu'une matrice de
% résultats avec toutes les données du système.

% Important: LBM et Generate_sample ont été modifiés pour en extraire les
% résultats, ce code ne fonctionnera pas sans les fichiers modifiés qui
% sont dans le dossier de remise du devoir.

clc; clear; close all;

% Définition des paramètres
seed_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
num_seeds = length(seed_values);
deltaP = 0.1;
poro = 0.9;
nombre_echantillon = 150;

% Fixer dx et NX
dx = 1e-6;
L = 200e-6;
NX = round(L / dx);

% Paramètres pour le diamètre des fibres
mean_fiber_d = 12.5;
std_fiber_d = 2.85;

% Générer des diamètres de fibres aléatoires suivant une distribution normale
seed_fiber_diameter = 123123;
rng(seed_fiber_diameter);
fiber_d_values = mean_fiber_d + std_fiber_d * randn(nombre_echantillon, 1);

% Initialisation des variables de stockage
d_equivalent_results = zeros(nombre_echantillon, num_seeds);
fiber_count_results = zeros(nombre_echantillon, num_seeds);
permeability_results = zeros(nombre_echantillon, num_seeds);
k_in_micron2_results = zeros(nombre_echantillon, num_seeds);
Re_results = zeros(nombre_echantillon, num_seeds);

% Boucle sur les diamètres des fibres
for d = 1:nombre_echantillon
    fiber_d = fiber_d_values(d);
    for s = 1:num_seeds
        seed = seed_values(s);
        [d_equivalent, num_fibers] = Generate_sample(seed, 'fiber_mat.tiff', fiber_d, 0, poro, NX, dx);
        [poro_eff, Re, k_in_micron2] = LBM('fiber_mat.tiff', NX, deltaP, dx, d_equivalent);
        
        % Stocker les résultats en tenant compte du diamètre
        d_equivalent_results(d, s) = d_equivalent;
        fiber_count_results(d, s) = num_fibers;
        permeability_results(d, s) = poro_eff;
        k_in_micron2_results(d, s) = k_in_micron2;
        Re_results(d, s) = Re;
        
        fprintf('Fiber Diameter %.2f, Seed %d - Porosité: %.5f, d_eq: %.2f, Fibres: %d, Perméabilité: %.5f, k: %.5f, Re: %.5f\r', ...
                fiber_d, seed, poro, d_equivalent, num_fibers, poro_eff, k_in_micron2, Re);
    end
end

% Calculer la moyenne de k_in_micron2_results pour chaque diamètre de fibre
mean_k_in_micron2_results = mean(k_in_micron2_results, 2); 

% Convertir les résultats en un vecteur
k_values = squeeze(mean_k_in_micron2_results);  

% Estimer les paramètres log-normaux pour k_in_micron2_results
params_k = fitdist(k_values, 'Lognormal'); 

% Estimer les paramètres normaux pour fiber_d_values
params_fiber_d = fitdist(fiber_d_values, 'Normal');

% Affichage des paramètres estimés pour k_in_micron2_results
disp(['Shape (sigma) de k_in_micron2_results: ', num2str(params_k.sigma)]);
disp(['Scale (mu) de k_in_micron2_results: ', num2str(params_k.mu)]);

% Tracer l'histogramme de k_in_micron2_results
figure;
subplot(2, 2, 1);
histogram(k_values, 20, 'Normalization', 'pdf', 'EdgeColor', 'k', 'FaceColor', [0.7, 0.7, 0.7]);
hold on;
x_k = linspace(min(k_values), max(k_values), 1000);
pdf_k_vals = pdf(params_k, x_k);
plot(x_k, pdf_k_vals, 'r-', 'LineWidth', 2);
title('PDF de la perméabilité');
xlabel('k (micron^2)');

% Tracer l'histogramme de fiber_d_values
subplot(2, 2, 2);
histogram(fiber_d_values, 20, 'Normalization', 'pdf', 'EdgeColor', 'k', 'FaceColor', [0.7, 0.7, 0.7]);
hold on;
x_fiber_d = linspace(min(fiber_d_values), max(fiber_d_values), 1000);
pdf_fiber_d_vals = pdf(params_fiber_d, x_fiber_d);  % Calcul de la PDF pour fiber_d_values
plot(x_fiber_d, pdf_fiber_d_vals, 'g-', 'LineWidth', 2);
title('PDF du diamètre des fibres');
xlabel('Diamètre des fibres (microns)');

% Tracer la CDF de k_in_micron2_results
subplot(2, 2, 3);
cdf_k_vals = cdf(params_k, x_k);
plot(x_k, cdf_k_vals, 'b-', 'LineWidth', 2);
title('CDF de la perméabilité');
xlabel('k (micron^2)');

% Tracer la CDF de fiber_d_values
subplot(2, 2, 4);
cdf_fiber_d_vals = cdf(params_fiber_d, x_fiber_d);
plot(x_fiber_d, cdf_fiber_d_vals, 'm-', 'LineWidth', 2);
title('CDF du diamètre des fibres');
xlabel('Diamètre des fibres (microns)');

% Statistiques supplémentaires pour fiber_d_values
mean_fiber_d = mean(fiber_d_values); 
std_fiber_d = std(fiber_d_values); 
disp(['Moyenne des diamètres des fibres: ', num2str(mean_fiber_d)]);
disp(['Ecart-type des diamètres des fibres: ', num2str(std_fiber_d)]);

% Création d'une matrice des résultats par seed
organized_results = [];

for s = 1:num_seeds
    for d = 1:nombre_echantillon
        d_equivalent = d_equivalent_results(d, s); 
        num_fibers = fiber_count_results(d, s);     
        poro_eff = permeability_results(d, s);  
        k_in_micron2 = k_in_micron2_results(d, s);  
        Re = Re_results(d, s);  
        organized_results = [organized_results; ...
            s, 0.9, d_equivalent, num_fibers, poro_eff, k_in_micron2, Re];
    end
end

% Afficher la matrice des résultats organisés
disp('Matrice des résultats organisés par seed :');
disp('Seed | Porosité | d_equivalent | num_fibers | permeability | k_in_micron2 | Re');
disp(organized_results);
