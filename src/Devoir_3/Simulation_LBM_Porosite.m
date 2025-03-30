% Ce code, basé sur le code donné launch_simulationLBM, permet d'itérer la
% simulation sur un grand nombre d'échantillons de porosité et de
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
mean_fiber_d = 12.5;
std_fiber_d = 0;
poro_moyenne = 0.9;
poro_ecart_type = 7.5e-3;
nombre_echantillon = 150;

% Fixer dx et NX
dx = 1e-6;
L = 200e-6;
NX = round(L / dx);

% Générer des porosités aléatoires
seed_porosity = 123123;
rng(seed_porosity);
poro_values = poro_moyenne + poro_ecart_type * randn(1, nombre_echantillon);  % Variations de porosité


% Initialisation des variables de stockage
d_equivalent_results = zeros(nombre_echantillon, num_seeds);
fiber_count_results = zeros(nombre_echantillon, num_seeds);
porosity_results = zeros(nombre_echantillon, num_seeds);
k_in_micron2_results = zeros(nombre_echantillon, num_seeds);
Re_results = zeros(nombre_echantillon, num_seeds);

% Boucle sur les porosités
for p = 1:nombre_echantillon
    poro = poro_values(p);
    for s = 1:num_seeds
        seed = seed_values(s);
        fiber_d = mean_fiber_d;

        [d_equivalent, num_fibers] = Generate_sample(seed, 'fiber_mat.tiff', mean_fiber_d, std_fiber_d, poro, NX, dx);
        [poro_eff, Re, k_in_micron2] = LBM('fiber_mat.tiff', NX, deltaP, dx, d_equivalent);

        d_equivalent_results(p, s) = d_equivalent;
        fiber_count_results(p, s) = num_fibers;
        porosity_results(p, s) = poro_eff;
        k_in_micron2_results(p, s) = k_in_micron2;
        Re_results(p, s) = Re;

        fprintf('Porosité %.5f, Seed %d - d_eq: %.2f, Fibres: %d, Perméabilité: %.5f, k: %.5f, Re: %.5f\r', ...
                poro, seed, d_equivalent, num_fibers, poro_eff, k_in_micron2, Re);
    end
end

% Calculer la moyenne de k_in_micron2_results pour chaque porosité
mean_k_in_micron2_results = mean(k_in_micron2_results, 2);

% Convertir les résultats en un vecteur
k_values = squeeze(mean_k_in_micron2_results);

% Estimer les paramètres log-normaux pour k_in_micron2_results
params_k = fitdist(k_values, 'Lognormal');

% Calculer la moyenne des porosités mesurées pour chaque échantillon
mean_porosities = mean(porosity_results, 2);

% Estimer les paramètres normaux pour la porosité
params_porosity = fitdist(mean_porosities, 'Normal'); 

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

% Tracer l'histogramme de porosity_results (mesurée)
subplot(2, 2, 2);
histogram(porosity_results(:), 20, 'Normalization', 'pdf', 'EdgeColor', 'k', 'FaceColor', [0.7, 0.7, 0.7]);
hold on;
x_poro = linspace(min(porosity_results(:)), max(porosity_results(:)), 1000);
pdf_poro_vals = pdf(params_porosity, x_poro);
plot(x_poro, pdf_poro_vals, 'g-', 'LineWidth', 2);
title('PDF de la porosité');
xlabel('Porosité');

% Tracer la CDF de k_in_micron2_results
subplot(2, 2, 3);
cdf_k_vals = cdf(params_k, x_k);
plot(x_k, cdf_k_vals, 'b-', 'LineWidth', 2);
title('CDF de la perméabilité');
xlabel('k (micron^2)');

% Tracer la CDF de porosity_results (mesurée)
subplot(2, 2, 4);
cdf_poro_vals = cdf(params_porosity, x_poro);
plot(x_poro, cdf_poro_vals, 'm-', 'LineWidth', 2);
title('CDF de la porosité');
xlabel('Porosité');

% Statistiques supplémentaires pour fiber_d_values
disp(['Moyenne des diamètres des fibres: ', num2str(mean_fiber_d)]);
disp(['Ecart-type des diamètres des fibres: ', num2str(std([mean_fiber_d]))]);

% Création d'une matrice des résultats par seed
organized_results = []; 

for s = 1:num_seeds
    for p = 1:nombre_echantillon
        d_equivalent = d_equivalent_results(p, s); 
        num_fibers = fiber_count_results(p, s);  
        poro_eff = porosity_results(p, s);   
        k_in_micron2 = k_in_micron2_results(p, s); 
        Re = Re_results(p, s);   
        organized_results = [organized_results; ...
            s, poro_values(p), d_equivalent, num_fibers, poro_eff, k_in_micron2, Re];
    end
end

% Afficher la matrice des résultats organisés
disp('Matrice des résultats organisés par seed :');
disp('Seed | Porosité donnée | Diamètre | Nombre de fibres | Porosité mesurée | k | Re');
disp(organized_results);


