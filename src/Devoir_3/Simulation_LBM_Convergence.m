% MATLAB script to launch a fiber structure generation and the corresponding LBM simulation
%
% INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place). 
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
% 
% NX: domain lateral size in grid cell


% La fonction de ce code est de fournir les données
% nécessaire au tableau Excel de la question A).
% Ce tableau Excel va nous permettre de tracer un graphique de
% convergence de la solution en fonction des maillages Nx choisit.

% Pour y parvenir avec ce code nous utilisons une double boucle :
% Elle permet d'évaluer la perméabilité de plusieurs seeds pour une même
% valeur de Nx afin de pouvoir reporter les valeurs sur Excel et évaluer la
% moyenne des perméabilités en ce Nx. La boucle passe ensuite au maillage Nx suivant
% et fait varier les seeds pour évaluer la perméabilité moyenne à ce Nx. On
% reproduit ce schéma pour chacun des maillages étudiés.

% Finalement on obtient un tableau Excel avec la moyenne des perméabilités
% pour chaque valeur de Nx.

clc;
clear all; 

deltaP= 0.1 ; % pressure drop in Pa
Nx_values = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400]; % Différents maillages testés 
dx_values = 0.0002 ./ Nx_values;  % Adaptation du pas en fonction de Nx
poro= 0.9 ; % Definition de la porosité
mean_fiber_d= 12.5 ; % Diamètre moyen des fibres en microns
std_d= 2.85 ; % Deviation sur les fibres en microns
filename= 'fiber_mat.tiff' ;

% Double boucle 

for i = 1:length(Nx_values) % Boucle sur les différents maillages choisit

    NX = Nx_values(i)
    dx = dx_values(i); 

for j = 100 : 106 % Boucle sur les différents seed choisit

seed = j

% generation of the fiber structure
[d_equivalent]=Generate_sample(seed,filename,mean_fiber_d,std_d,poro,NX,dx);

% calculation of the flow field and the permeability from Darcy Law
LBM(filename,NX,deltaP,dx,d_equivalent);

end 

end




