function [T_FEM, q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, H, num_z, flag)
    % Définir la géométrie
    g = decsg([3 4 0 0 D/2 D/2 L 0 0 L]');
    
    % Visualiser la géométrie si flag est vrai
    if flag
        figure;
        pdegplot(g, EdgeLabels="on");
        axis equal;
    end
    
    % Créer le modèle FEM
    model = femodel(AnalysisType="thermalSteady", Geometry=g);
    model.PlanarType = "axisymmetric";
    
    % Définir les propriétés des matériaux
    model.MaterialProperties = materialProperties(ThermalConductivity=k);
    model.FaceLoad = faceLoad(Heat=0);
    
    % Conditions aux limites
    model.EdgeBC(2) = edgeBC(Temperature=Tm);  % Condition limite de température sur le bord 2
    model.EdgeLoad(3) = edgeLoad(ConvectionCoefficient=h, AmbientTemperature=T_inf);  % Convection sur le bord 3
    model.EdgeLoad(4) = edgeLoad(ConvectionCoefficient=h, AmbientTemperature=T_inf);  % Convection sur le bord 4

    % Générer le maillage
    model = generateMesh(model, 'Hmax', H);
    
    % Résoudre le modèle
    result = solve(model);
    
    % Extraire les résultats de température
    T = result.Temperature;
    
    % Visualiser la distribution de température si flag est vrai
    if flag
        figure;
        pdeplot(result.Mesh, XYData=T, Contour="on");
        axis equal;
        title("Température à l'état stationnaire");
        figure;
        pdeplot(model.Mesh);
        axis equal;
        title("Maillage FEM");
    end
    
    % Stocker les coordonnées des nœuds et les températures dans le maillage
    Mesh = [result.Mesh.Nodes(1,:)', result.Mesh.Nodes(2,:)', T];
    
    % Définir les positions des tranches axiales (valeurs z pour la moyenne)
    z = linspace(0, L, num_z);
    T_FEM = zeros(length(z), 1);
    tolerance = 5e-4;
    
    % Boucle sur chaque tranche z
    for i = 1:length(z)
        z_val = z(i);
        nodes_at_z = Mesh(abs(Mesh(:, 2) - z_val) < tolerance, :);
        
        % Calcul de la température moyenne pondérée
        r = nodes_at_z(:, 1);  % Coordonnée radiale
        T = nodes_at_z(:, 3);  % Température au nœud
        weighted_sum = sum(r .* T);
        sum_r = sum(r);
        T_FEM(i) = weighted_sum / sum_r;
    end
    
    % Tracer le profil de température si flag est vrai
    if flag
        figure;
        plot(z, T_FEM);
        xlabel('Position axiale (z)');
        ylabel('Température moyenne pondérée');
        title('Profil de température le long du cylindre');
    end

    % Calcul du flux de chaleur
    gradient = (T_FEM(2) - T_FEM(1)) / (z(2) - z(1));
    Ac = (pi * D^2) / 4;
    q_FEM = -k * Ac * gradient;  % Flux de chaleur
    
end


