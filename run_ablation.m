setup = DE.defaultSetup();
setup.varCount = 4;
setup.goal = 'min';
setup.fitnessFunction = 'DE.kinetic_ablation';
setup.rounds = 1000000; 
setup.populationAmount = 350;
setup.varMin = -1;
setup.varMax = 1;
setup.disturb = 1;
setup.avoidStagnation = 0;
setup.stagnationSteps = 2000;
setup.initializeAroundBest = 1;
setup.crossoverRate = 0.6;
setup.mutationScale = 0.1;
% setup.bestCromo = eliteX;
%%fit 0.00078368
% setup.bestCromo = [-2756.0146991240494571684394031763, 8.2726068408982946777996403397992, 7.8327701190963949784418218769133, 664.76399435408973204175708815455, 1024.1393489834906631585909053683, -14.29944068769194842616343521513];
% setup.bestCromo = [3 22 0.85 650 40 0.15];
%[ -0.72757368577982473034637678210856, 0.013018421302918137261173114893609, 0.87684267235373081206262213527225, 0.93737538208561366204207843111362]
 


% setup.bestCromo = [5 13 0.8 600 40 0.05];
simu = DE(setup);
result = simu.run()