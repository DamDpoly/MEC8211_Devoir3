clear
clc

% Données d'entrée
D_mean = 0.050;  % Diamètre du pilier (m)
L_mean = 0.400;  % Longueur du pilier (m)
k_mean = 150;    % Conductivité thermique (W/m.K)
h_mean = 15;     % Coefficient d'échange thermique (W/m².K)
T_inf_mean = 300; % Température ambiante (°C)
Tm_mean = 400;   % Température à la base (°C)

%% Création des données d'entrée
% Incertitudes
std_D = 0.001;
std_L = 0.005;
std_k = 1;
std_h = 2;
std_T_inf = 2;
std_Tm = 10;

num_samples = 3000;

% Échantillonnage positif
D_samples     = sample_positive(D_mean, std_D, num_samples);
L_samples     = sample_positive(L_mean, std_L, num_samples);
k_samples     = sample_positive(k_mean, std_k, num_samples);
h_samples     = sample_positive(h_mean, std_h, num_samples);
T_inf_samples = sample_positive(T_inf_mean, std_T_inf, num_samples);
Tm_samples    = sample_positive(Tm_mean, std_Tm, num_samples);
Biot_samples = (h_samples .* D_samples) ./ (4 .* k_samples);

% Initialisation des résultats
T_FEM_results = zeros(num_samples, 50);
Q_FEM_results = zeros(num_samples, 1);

T_FDS_results = zeros(num_samples, 50);
Q_FDS_results = zeros(num_samples, 1);

% Boucle de simulation
for i = 1:num_samples
    D = D_samples(i);
    L = L_samples(i);
    k = k_samples(i);
    h = h_samples(i);
    T_inf = T_inf_samples(i);
    Tm = Tm_samples(i);

    % FEM
    [T_FEM, Q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, 0.0016, 50, false);
    T_FEM_results(i, :) = T_FEM(:);
    Q_FEM_results(i) = Q_FEM;

    % FDS
    [T_FDS, ~, Q_FDS, ~] = Solution_numerique_ailette(D, L, k, h, T_inf, Tm, 50);
    T_FDS_results(i, :) = T_FDS(:);
    Q_FDS_results(i) = Q_FDS;
end

%% Histogrammes PDF des paramètres d'entrée
param_names = {'D', 'L', 'k', 'h', 'T_{inf}', 'T_m'};
param_samples = {D_samples, L_samples, k_samples, h_samples, T_inf_samples, Tm_samples};

figure('Name', 'PDF des paramètres d''entrée');
for i = 1:length(param_samples)
    subplot(2, 3, i);
    samples = param_samples{i};
    histogram(samples, 20, 'Normalization', 'pdf', 'FaceColor', [0.7, 0.7, 0.9]);
    hold on;
    pd = fitdist(samples, 'Normal');
    x = linspace(min(samples), max(samples), 1000);
    plot(x, pdf(pd, x), 'r', 'LineWidth', 2);
    title(['PDF de ', param_names{i}]);
    xlabel(param_names{i});
    legend('Échantillons', 'Normal Fit');
    grid on;
end

%% CDF des paramètres d'entrée
figure('Name', 'CDF des paramètres d''entrée');
for i = 1:length(param_samples)
    subplot(2, 3, i);
    samples = param_samples{i};
    pd = fitdist(samples, 'Normal');
    x = linspace(min(samples), max(samples), 1000);
    plot(x, cdf(pd, x), 'b', 'LineWidth', 2);
    title(['CDF de ', param_names{i}]);
    xlabel(param_names{i});
    grid on;
end

%% Analyse de Q_FEM et Q_FDS
pd_Q_FEM = fitdist(Q_FEM_results, 'Lognormal');
pd_Q_FDS = fitdist(Q_FDS_results, 'Lognormal');

disp('Paramètres log-normaux de Q_{FEM} :');
disp(['Mu (log): ', num2str(pd_Q_FEM.mu)]);
disp(['Sigma (log): ', num2str(pd_Q_FEM.sigma)]);

disp('Paramètres log-normaux de Q_{FDS} :');
disp(['Mu (log): ', num2str(pd_Q_FDS.mu)]);
disp(['Sigma (log): ', num2str(pd_Q_FDS.sigma)]);

% PDF et CDF comparatifs
figure('Name', 'Comparaison des distributions de Q_{FEM} et Q_{FDS}');

subplot(1,2,1); % PDF
histogram(Q_FEM_results, 20, 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'FaceColor', [0.2 0.6 0.8]);
hold on;
histogram(Q_FDS_results, 20, 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'FaceColor', [0.8 0.5 0.2]);
x = linspace(min([Q_FEM_results; Q_FDS_results]), max([Q_FEM_results; Q_FDS_results]), 1000);
plot(x, pdf(pd_Q_FEM, x), 'b-', 'LineWidth', 2);
plot(x, pdf(pd_Q_FDS, x), 'r-', 'LineWidth', 2);
title('PDF comparée');
xlabel('Q (W)');
ylabel('Densité de probabilité');
legend('FEM Samples','FDS Samples','FEM Fit','FDS Fit');
grid on;

subplot(1,2,2); % CDF
plot(x, cdf(pd_Q_FEM, x), 'b-', 'LineWidth', 2); hold on;
plot(x, cdf(pd_Q_FDS, x), 'r-', 'LineWidth', 2);
title('CDF comparée');
xlabel('Q (W)');
ylabel('Fonction de répartition');
legend('FEM','FDS');
grid on;

%% Construction de la matrice finale pour analyse de sensibilité

% Matrice : [D, L, k, h, T_inf, Tm, Biot, Q_FEM, Q_FDS]
results_matrix = [ ...
    D_samples, ...
    L_samples, ...
    k_samples, ...
    h_samples, ...
    T_inf_samples, ...
    Tm_samples, ...
    Biot_samples, ...
    Q_FEM_results, ...
    Q_FDS_results ...
];

% Noms des colonnes pour référence
results_headers = {'D', 'L', 'k', 'h', 'T_inf', 'Tm', 'Biot', 'Q_FEM', 'Q_FDS'};

% Aperçu
disp('Aperçu des 5 premières lignes du tableau résultats :');
disp(array2table(results_matrix(1:5,:), 'VariableNames', results_headers));

%% Fonction d’échantillonnage positive
function samples = sample_positive(mean_val, std_val, num_samples)
    samples = zeros(num_samples, 1);
    i = 1;
    while i <= num_samples
        val = normrnd(mean_val, std_val);
        if val > 0
            samples(i) = val;
            i = i + 1;
        end
    end
end

