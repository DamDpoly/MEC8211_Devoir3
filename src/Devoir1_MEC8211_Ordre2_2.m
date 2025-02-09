%{ 
Ce code permet la résolution de l'équation de diffusion du sel suivant
un modèle stationnaire. 
La discrétisation spatial du pilier est gérable avec le paramètre Ntot,
qualifiant le nombre total de noeuds du système, par défault il est de 5.
Les paramètres liés aux dimensions du pilier et aux concentrations de sels 
sont constants.
%}

%Ce code utilise des schémas d'approximations des dérivées du second ordre
%pour résoudre l'EDP.
clc

% Definition des constantes de l'exercice
Deff = 10^(-10);  
S = 2*10^(-8);
Ce = 20; 
D = 1; 
Rayon = D /2;
alpha = S / Deff;

% Definition du nombre de noeuds de discrétisation
Ntot = 5;

% Calcul du pas de discrétisation dr
dr = Rayon / (Ntot - 1);

%Definition des matrices A et B pour pouvoir calculer A C = B

A = zeros(Ntot-2,Ntot-2);

B = zeros(Ntot-2,1);
    
%Remplissage de A

%Condition limite C1 = C0

A(1, 1) = (1-0.5) -(2*1); % Coefficient pour A_[1,1]
A(1, 1+1) = 1 + 0.5; % Coefficient pour A_[1,2]

%Condition limite CN-1 = CE

 A(Ntot-2,Ntot-2-1) = (Ntot-2)-0.5; % Coefficient pour A_[Ntot-2,Ntot-3]
        
 A(Ntot-2, Ntot-2) = - (2*(Ntot-2)); % Coefficient pour A_[Ntot-2,Ntot-2]
        
%Boucle de remplissage de A
     
for i = 2:Ntot-3

 
        A(i,i-1) = i - 0.5;       % Coefficient pour A_[i,i-1]
        
        A(i, i) = - (2*i); % Coefficient pour A_[i,i]
        
        A(i, i+1) = i + 0.5; % Coefficient pour A_[i,i+1]
        
end


fprintf('Matrice A :\n');
disp(A);


%Remplissage de B

%Condition limite CN-1 = CE avec CE transféré au second membre

B(Ntot-2,1)= alpha * dr * dr * (Ntot-2)-(Ntot-1.5)*Ce;

% Boucle de remplissage de B

for i=1:Ntot-3
B(i,1)= alpha * dr * dr * (i); % Cas global 
end

fprintf('Matrice B :\n');
disp(B);

% Résolution du système 

C = A\B; 

% Concaténation des résultats calculés avec ceux des conditions limites

C = [C(1); C; Ce];

fprintf('Les concentrations calculées numériquement sont les suivantes :\n ');
for i = 1:length(C)
    fprintf('C%d : %.4f  ', i-1, C(i));
end


% Définition de la position de chaque noeuds dans le vecteur r
r = linspace(0, Rayon, Ntot);

% Calcul de la concentration avec la solution analytique
C_analytique = (1/4) * (alpha) * Rayon^(2) * ((r.^(2) / Rayon^(2)) - 1) + Ce;

fprintf('\n Les concentrations calculées analytiquement sont les suivantes :\n ');
for i = 1:length(C)
    fprintf('C%d : %.4f  ', i-1, C_analytique(i));
end


% Tracé de la comparaison des concentrations analytiques et numériques
figure;
plot(r, C, '-o', 'DisplayName', 'Résultat Numérique');
hold on;
plot(r, C_analytique, '-x', 'DisplayName', 'Résultat Analytique');
xlabel('Rayon r (m)');
ylabel('Concentration C (mol/m^3)');
title('Comparaison entre les résultats numériques et analytiques');
legend('show');
grid on;




