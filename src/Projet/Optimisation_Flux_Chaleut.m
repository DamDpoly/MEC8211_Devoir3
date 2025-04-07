% Define the parameters
D = 0.05;            % Diameter (m)
Longueur = 0.2;      % Length (m)
k_cuivre = 398;      % Thermal conductivity (W/m.K)
h = 500;             % Heat transfer coefficient (W/m².K)
T_inf = 25;          % Ambient temperature (°C)
Tm = 100;            % Temperature at the base (°C)
Ntot = 50;           % Number of nodes for spatial discretization

% Call the function to compute the heat flux at x=0
Calculate_Heat_Flux(D, Longueur, k_cuivre, h, T_inf, Tm, Ntot);
