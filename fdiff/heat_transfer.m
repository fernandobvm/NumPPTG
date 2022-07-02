thermalmodel = createpde('thermal','transient');
g = @squareg;
geometryFromEdges(thermalmodel,g);
thermalProperties(thermalmodel,'ThermalConductivity',79.5,...
                               'MassDensity',7850,...
                               'SpecificHeat',450,...
                               'Face',1)
internalHeatSource(thermalmodel,25);
thermalBC(thermalmodel,'Edge',[1,3,4],'HeatFlux',0);
thermalBC(thermalmodel,'Edge',2,...
                       'ConvectionCoefficient',5000,...
                       'AmbientTemperature',25);
thermalIC(thermalmodel,25);
thermalIC(thermalmodel,100,'Edge',4);
thermalmodel.StefanBoltzmannConstant = 5.670367e-8;
mesh = generateMesh(thermalmodel);
result = solve(thermalmodel,[0:0.1:10])