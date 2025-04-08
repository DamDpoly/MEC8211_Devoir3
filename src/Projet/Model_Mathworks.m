g = decsg([3 4 0 0 .5 .5 1 0 0 1]');

figure
pdegplot(g,EdgeLabels="on")
axis equal

model = femodel(AnalysisType="thermalSteady", ...
                Geometry=g);

model.PlanarType = "axisymmetric";

k = 10; % Thermal conductivity, W/(m*C)

model.MaterialProperties = materialProperties(ThermalConductivity=k);
model.FaceLoad = faceLoad(Heat=0);

model.EdgeBC(2) = edgeBC(Temperature=100);
model.EdgeLoad(3) = edgeLoad(ConvectionCoefficient=100, AmbientTemperature=25);

% Refining the mesh
model = generateMesh(model, 'Hmax', 0.005);

% Visualize mesh
figure
pdemesh(model)
axis equal

% Solve the model
result = solve(model);
T = result.Temperature;

% Visualize the result
figure
pdeplot(result.Mesh, XYData=T, Contour="on")
axis equal
title("Steady-State Temperature")

% Solve the model
result = solve(model);

% Extract temperature values
T = result.Temperature;  % Temperature at each node (should be 1xN)

% Extract mesh node coordinates
coords = result.Mesh.Nodes;  % coords should be 2xN (1st row = x, 2nd row = y)

% Ensure the coordinates and temperature are column vectors
x_coords = coords(1,:)';  % Transpose to make it N x 1
y_coords = coords(2,:)';  % Transpose to make it N x 1
temperature = T;  % Transpose to make it N x 1

% Combine coordinates and temperature into a single matrix
temperatureMatrix = [x_coords, y_coords, temperature];

% Display the matrix
disp('Node Coordinates and Temperature:');
disp(temperatureMatrix);
