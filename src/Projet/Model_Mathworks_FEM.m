function [T_FEM, Q_FEM] = Model_Mathworks_FEM(D, L, k, h, T_inf, Tm, H, num_z, flag)
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
    model.EdgeBC(2) = edgeBC(Temperature=Tm);  % Température imposée à la base
    model.EdgeLoad(3) = edgeLoad(ConvectionCoefficient=h, AmbientTemperature=T_inf);  % Convection latérale
    model.EdgeLoad(4) = edgeLoad(ConvectionCoefficient=h, AmbientTemperature=T_inf);  % Convection extrémité

    % Générer le maillage
    model = generateMesh(model, 'Hmax', H);
    
    % Résoudre le modèle
    result = solve(model);
    
    % Extraire la température
    T = result.Temperature;
    
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

    % Stocker coordonnées et températures
    Mesh = [result.Mesh.Nodes(1,:)', result.Mesh.Nodes(2,:)', T];

    % Moyenne pondérée selon r pour obtenir un profil axial 1D
    z = linspace(0, L, num_z);
    T_FEM = zeros(length(z), 1);
    tolerance = 5e-4;

    for i = 1:length(z)
        z_val = z(i);
        nodes_at_z = Mesh(abs(Mesh(:, 2) - z_val) < tolerance, :);
        r = nodes_at_z(:, 1);
        T_vals = nodes_at_z(:, 3);
        weighted_sum = sum(r .* T_vals);
        sum_r = sum(r);
        if sum_r > 0
            T_FEM(i) = weighted_sum / sum_r;
        else
            T_FEM(i) = NaN;
        end
    end

    if flag
        figure;
        plot(z, T_FEM);
        xlabel('Position axiale (z)');
        ylabel('Température moyenne pondérée');
        title('Profil de température le long du cylindre');
    end

    % Calcul du débit de chaleur total Q [W]
    gradient = (T_FEM(2) - T_FEM(1)) / (z(2) - z(1));
    Ac = (pi * D^2) / 4;
    Q_FEM = -k * Ac * gradient;  % Puissance thermique [W]
end

