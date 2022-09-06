classdef DE
    
    properties (GetAccess = 'public', SetAccess = 'private')
        setup = 0;
    end
    
    methods
        
        function this = DE(setup)
            this.setup = setup;
        end

        function R = run(this)
                
            rn = 0;
            gc = 1; 
            this.validateSetup();

            task = getCurrentTask;
            taskID = 1; %task.ID;
            
            %Matriz que hospedará a população inicial.
            pop = this.initiPopulation();

            %Marcação de tempo para debug...
            dtIni = now; 
            
            %Inicializa vetor de desempenho...
            D = zeros(this.setup.rounds, 2);
            
            %Executa a função de avaliação de fitness para primeira geradação da população.
            fits = this.evalFitness(pop);
            if this.setup.debug                            
                %fprintf('[');
                fprintf('%d -- [', taskID);
            end

            %Executa as iterações para produção de novas gerações...
            for g = 1:this.setup.rounds
                
                %Itera por todos os indivíduos da população...
                popLen = size(pop, 1);
                for k = 1:popLen
                    %Aplica a mutação...
                    %trial = this.applyMutation(pop, k);
                    trial = applyMutation(this, pop, k);

                    %Gera a novo filho através de crossover...
                    %child = this.applyCrossover(trial, pop(k, :));
                    child = applyCrossover(this, trial, pop(k, :));
                
                    %Compara filho gerado com parent original...
                    fit = this.evalFitness(child);
                    if this.best(fit, fits(k))
                        pop(k, :) = child;
                        fits(k) = fit;
                        
                        %if g > 1 && (fit > D(g-1, 1))
                        %    %fprintf('BINGO!! (%s).\n', fit);
                        %    fprintf('\n %s \n',char(vpa(152.6201399999999921419657766819 - fit)))
                        %end
                    end
                    
                    if mod(g, 1) == 0 && this.setup.debug
                        perc = (k/popLen)*100;
                        if mod(perc, 20) == 0
                            fprintf('%d%%,', perc)
                        end
                    end
                end
                
                %Executa a coleta do elitismo...
                ids = this.evalElitism(fits);
                elite = pop(ids, :);
                D(g, 1) = fits(ids(1)); 
                D(g, 2) = sum(fits)/length(fits);
                %D(g, 3) = elite; %apenas para varCount == 1
                
                %assignin('base', 'D', D);
                %assignin('base', 'pop', pop);
                %assignin('base', 'fits', fits);
                assignin('base', 'eliteX', elite);
                
                %if isfield(this.setup, 'onNewGeneration')
                %    this.setup.onNewGeneration(g, R, elite, D(g, 1), rem(now, 1)-dtIni);
                %end
                
                %Debug..
                if mod(g, 1) == 0 && this.setup.debug
                    gi = DE.iif(g > 100, g-100, 1);
                    disp(strcat(...
                        '] DEBUG-DE ==>  G=', int2str(g), ' P=', int2str(size(pop, 1)),...
                        ' T=', datestr(rem(now, 1)-dtIni, 'HH:MM:SS:FFF'),...
                        ' F=', num2str(D(g, 1)),... %char(vpa(D(g, 1),5)), ...
                        ' M=', num2str(D(g, 2)),...
                        ' D=', num2str(abs(D(g, 1) - D(g, 2))),...
                        ' Dif=', num2str(sum(diff(D(gi:g,1))))...
                    ));
                    
                    %fprintf(char(vpa(152.6201399999999921419657766819 - D(g, 1))));
                    %fprintf('%s \n', char(vpa(152.6201399999999921419657766819 - D(g, 1))));
                    %fprintf(char(vpa(D(g, 1), 32)));

                    %fprintf('[');
                    %pack; %TODO: Remover...
                    fprintf('%d -- [', taskID);
                    dtIni = now;
                    %disp(vpa(elite));
                end                
                %disp(pop);
                
                %{
                if ~isempty(this.setup.eliteFile)
                    %save this.setup.eliteFile  elite -ASCII;
                    %save(this.setup.eliteFile, 'elite');
                    fl = strcat('results\', ...
                        datestr(now, 'yyyy-mm-dd-hh-MM-ss_'),...
                        num2str(D(g, 1)), '.dat');
                    save(fl, 'elite');
                    
                    %datafile = fopen(fl,'a');
                    %fwrite(datafile, elite);
                    %fprintf(datafile, '%s,', num2str(elite));
                    %fclose(datafile);
                end
                %}
                
                %TODO:Parametrizar...
                %%{O ideal é essa distância ser influenciada pela distãncia
                %da solução atual para o ótimo buscado, quanto mais longe
                %tiver do ótimo, menor deve ser o trashhold. Se tiver longe
                %pode-se usar 1e-4 ou 1e-6, quando tiver longe pode-se usar
                %1e-8 ou 1e-12.
                if isfield(this.setup, 'disturb')
                    if this.setup.disturb && (abs(D(g, 1) - D(g, 2)) < this.setup.disturbThreshold)%1e-12) 
                        if rn < 10
                            rn = rn + 1;
                        else
                            rn = 0;
                            
                            opt = randperm(3); opt = opt(1);
                            %opt = randperm(4); opt = opt(1);
                            [pop, fits] = this.regeneratePopulation(pop, fits, elite, opt);
                        end
                    end
                end
                %}
                
                %Avoid Stagination...
                if this.setup.avoidStagnation
                    gR = this.setup.stagnationSteps;
                    gi = DE.iif(g > gR, g-gR, 1);
                    df = sum(diff(D(gi:g,1)));
                    
                    %if (g > gc*gR) && (df < 1e-8)
                    if (g - gc > gR) && (df < 1e-8)
                        if this.setup.debug
                            disp('### --- Opsss, detectada estagnação!!!');
                        end                        
                        opt = randperm(3); opt = opt(1);
                        [pop, fits] = this.regeneratePopulation(pop, fits, elite, opt);

                        gc = g;
                    end
                end
                
                %TODO: Parametrizar....
                if isfield(this.setup, 'limit')
                    if abs(this.setup.limit - D(g,1)) <= 1e-12 %1e-9
                    %if D(g,1) >= this.setup.limit
                    %if this.best(D(g,1), this.setup.limit)
                        disp('Valor limite alcançado, encerrando busca... :)');
                        break;
                    end
                end
            end
            
            R.eliteValue = D(g,1);
            R.elites = D(1:g,1);
            R.bests = D(1:g,1);
            R.means = D(1:g,2);
            %R.elites = D(1:g,3); %apenas para varCount == 1
            R.eliteCromo = elite;
            R.lastPop = pop;

            %if ~isempty(this.setup.eliteFile)
            %    %save this.setup.eliteFile  elite -ASCII;
            %    save(this.setup.eliteFile, 'elite');
            %end
            
            %disp(parcluster);
        end
        
    end
    
    methods (Access = 'private')
        
        function validateSetup(this)
            
            if strcmp(this.setup.fitnessFunction, 'null')
                error 'Deve ser especificada uma função fitness, sorry!';
            end
            
        end

        function [pop, fits] = regeneratePopulation(this, pop, fits, elite, opt)
            switch opt
                case 1
                    if this.setup.debug
                        fprintf('Opsss, refazendo população.... \n1 -- [');
                    end
                    pop = this.initiPopulation();
                case 2
                    if this.setup.debug
                        fprintf('Opsss, pertubando população.... \n1 -- [');
                    end
                    pop = this.disturbPopulation(pop, 0.5); %0.25);
                case 3
                    if this.setup.debug
                        fprintf('Opsss, invertendo população.... \n1 -- [');
                    end
                    pop = this.opposedPopulation(pop);
                case 4
                    if this.setup.debug
                        fprintf('Opsss, refazendo população NAROUND.... \n1 -- [');
                    end
                    pop = this.initiPopulation(true);
                otherwise
                    error 'Unknow option, sorry! :(';
            end
            
            pop(1,:) = elite;
            fits = this.resetFits(fits);
            fits(1) = this.evalFitness(elite);
        end

        function pop = opposedPopulation(this, pop)
            len = size(pop,1);
            %piv = (min(pop) + max(pop));
            piv = 2*pi*ones(1, size(pop, 2));
            pop = repmat(piv, len, 1) - pop;
        end
        
        function f = resetFits(this, f)
            switch this.setup.goal
                case {'max', 'Max', 'MAX'}
                    f(1:end) = -realmax;
                case {'min', 'Min', 'MIN'}
                    f(1:end) = realmax;
            end            
        end        
        
        function u = applyMutation(this, pop, k)
            
            len = size(pop, 1);
            vec = randperm(len);
            all = vec(vec ~= k); all = all(1:3);
            
            x1 = pop(all(1), :); x2 = pop(all(2), :); x3 = pop(all(3), :);
            
            u = x1 + (x2 - x3) * this.setup.mutationScale;
        end

        function u = applyMutationU(this, pop, k)            
            len = size(pop, 1);
            vec = randperm(len);
            n = this.setup.varCount;
            all = vec(vec ~= k); all = all(1:3);
            
            x1 = pop(all(1), :); x2 = pop(all(2), :); x3 = pop(all(3), :);
            
            x1 = vec2mat(x1, sqrt(n))';
            x2 = vec2mat(x2, sqrt(n))';
            x3 = vec2mat(x3, sqrt(n))';
            
            u = x1 + (x2 - x3) * this.setup.mutationScale;
            u = u/norm(u);
            
            %u = x1 * (x2 * x3) * x1';% .* this.setup.mutationScale;
            u = reshape(u, 1, n);
        end
        
        function F = evalFitness(this, pop)
            
            len = size(pop, 1);
            setup_ = this.setup;
            fitnessFunction = this.setup.fitnessFunction;
            
            if (this.setup.vectorizedFitness)
                F = feval(fitnessFunction, setup_, pop);
            else                
                F = zeros(len, 1);
                for k = 1:len
                    F(k) = feval(fitnessFunction, setup_, pop(k,:));
                end                                
            end                                   
        end
        
        function E = evalElitism(this, Fits) 
            
            E = ones(1, this.setup.elitismAmount);

            for l = 1:this.setup.elitismAmount

                switch this.setup.goal
                    case {'max', 'Max', 'MAX'}
                        [~, id] = max(Fits);
                    case {'min', 'Min', 'MIN'}
                        [~, id] = min(Fits);
                end

                E(l) = id;
                Fits(id) = [];
            end
                
        end
        
        function c = applyCrossover(this, trial, current)
            
            len = length(trial);
            vec = rand(1, len) < this.setup.crossoverRate;
            
            %Caso nenhum tenha sido selecionado, será forçado um...
            if sum(vec) == 0
                p = randperm(len);
                vec(p(1)) = 1;
            end
            
            %Faz o crossover do vetor trial com o current, em que vec indica 
            %a probabilidade de usar uma variável do vetor trial...
            c = vec .* trial + ~vec .* current;
            
            %TODO:Remover, caso específico...
            %U = class_vec_to_mat([c 0]);
            %if Kron.isKron2(U) == 1
            %    c = rand(1, 3);
            %end
        end
       
        function c = applyCrossoverU(this, trial, current)

            n = this.setup.varCount;
            vec = rand(sqrt(n)) < this.setup.crossoverRate;
            
            %Caso nenhum tenha sido selecionado, será forçado um...
            if sum(vec) == 0
                p = randperm(n);
                vec(p(1)) = 1;
            end
            
            trial = vec2mat(trial, sqrt(n))';
            current = vec2mat(current, sqrt(n))';

            %Faz o crossover do vetor trial com o current, em que vec indica 
            %a probabilidade de usar uma variável do vetor trial...
            c = vec .* trial + ~vec .* current;
            c = c/norm(c);
            c = reshape(c, 1, n);
        end
        
        function P = initiPopulation(this, naround)
            if nargin == 1
                naround = false;
            end
            varMin = this.setup.varMin;
            varMax = this.setup.varMax;
            P = varMin + (varMax-varMin) .* ...
                rand(this.setup.populationAmount, this.setup.varCount);

            %TODO:Parametrizar...
            if this.setup.complexVar
                P = P + round(rand(this.setup.populationAmount, this.setup.varCount)) .* ...
                    rand(this.setup.populationAmount, this.setup.varCount) * 1i;
            end
            %TODO:Remover, caso específico...
            %P(:,1:40) = real(P(:,1:40));
            
            %Iniciar a população com pertubações diversas do melhor cromossomo.
            if ~naround && this.setup.initializeAroundBest 
                if ~isempty(this.setup.bestCromo)
                    P = repmat(this.setup.bestCromo, this.setup.populationAmount, 1);
                    P = this.disturbPopulation(P, 1); %0.25);
                    P(end,:) = this.setup.bestCromo;
                end
            end
            
            if ~isempty(this.setup.bestCromo)
                P(end,:) = this.setup.bestCromo;
                this.setup.bestCromo = [];
                %P = [P; this.setup.bestCromo];
            end
        end
        
        function P = initiPopulationU(this)
            
            m = this.setup.populationAmount;
            n = this.setup.varCount;
            P = zeros(m, n);

            for k = 1:m                
                U = Telp.genUniRandom(sqrt(n));                
                P(k,:) = reshape(U, 1, n);                
            end
            
            if ~isempty(this.setup.bestCromo)
                P(end,:) = this.setup.bestCromo;
                this.setup.bestCromo = [];
                %P = [P; this.setup.bestCromo];
            end
        end        
        
        function P = disturbPopulation(this, pop, prob)
            
            varCount = this.setup.varCount;
            popAmount = this.setup.populationAmount;
            
            %4 * this.setup.mutationScale * ...
            %rand * 10 * this.setup.mutationScale * ...
            if this.setup.complexVar
                P = pop + ...
                    this.setup.mutationScale * ...
                    ((rand(popAmount, varCount) <= prob) .* ...
                    (rand(popAmount, varCount) * 1i));
            else
                P = pop + ...
                    this.setup.mutationScale * ...
                    ((rand(popAmount, varCount) <= prob) .* ...
                    (rand(popAmount, varCount)));
            end
        end
        
        function r = best(this, i1, i2)
           
            switch this.setup.goal
                case {'max', 'Max', 'MAX'}
                    r = i1 > i2;
                case {'min', 'Min', 'MIN'}
                    r = i1 < i2;
            end            
            
        end
        
    end
    
    methods (Static)
        
        function setup = defaultSetup()
        %DEFAULTSETUP Fornece um conjunto default de parametrizações para a estrutura
        % genérica de otimização baseada em evolução diferencial.
            
            %Quantidade de variáveis integrantes de um cromossomo.
            setup.varCount = 14; 

            %Menor valor para o cromossomo.
            setup.varMin = 0; 

            %Maior valor para o cromossomo.
            setup.varMax = 1; 

            %Maior valor para o cromossomo.
            setup.complexVar = false; 

            %Quantidade de cromossomos da população.
            setup.populationAmount = 150;

            %Percentual de mudança no processo de corssover, em %.
            setup.crossoverRate = 0.65;

            %Fator de incremento da mutação diferencial.
            setup.mutationScale = 0.5; 

            %Inclusão do cromossomo UM... 
%             setup.bestCromo = [0 3.4e-8 0 0.074 1360 0 0.03 0 2e-6 0.03 0.01 0.006 2.5e-6 1e-8];
            setup.bestCromo = [];

            %Quantidade de rounds que serão rodados pelo GA.
            setup.rounds = 30;

            %Indica o objetivo da otimização, se minimizar ou maximizar.
            setup.goal = 'max'; %'max';
            
            %Indica se a função fitness usada é vetorizada ou não. 
            setup.vectorizedFitness = 0; %0;
            
            %Função para avaliação de fitness.
            setup.fitnessFunction = 'DE.ppt';
            
            %Indica a quantidade de indivíduos devem ser usados no elitismo.
            setup.elitismAmount = 1;
            
            %Indica se devem ser emitidas mensagens de depuração.
            setup.debug = true;
            
            %Indica se a população deve ser regerada quando a média tiver a
            %uma distância menor que disturbThreshold em relação ao melhor
            %indivíduo.
            setup.disturb = 1;
            
            %Determina o limiar para pertubação da população em situações 
            %de convergência precoce.
            setup.disturbThreshold = 1e-10;
                    
            %Iniciar a população com pertubações diversas do melhor cromossomo.
            setup.initializeAroundBest = 0;
            
            %Indica se deve pertubar a população quando detectada estagnação.
            setup.avoidStagnation = 1;
            setup.stagnationSteps = 100;
            
            %Handler de função a ser chamada ao final de cada geração...
            %setup.onNewGeneration = null;
            
            setup.MinValues = [];
            
            setup.MaxValues = [];
            
            setup.idx = [];
            
            setup.ppt_geometry = 0;
        end        
        
        function r = iif(c, t, f)
            
            if c
                r = t;
            else
                r = f;
            end
            
        end
                
    end
    
    %TODO:Teste, depois remover...
    methods (Static)
                
        function R = runTest(best)
            %Configurando os parâmetros do SW...
            stp = DE.defaultSetup();
            
            if nargin == 1
                stp.One.bestCromo = best; 
            end
            stp.limit = 0;
            stp.debug = 1;
            stp.disturb = 0;
            stp.goal = 'min';
            stp.varCount = 2;
            stp.rounds = 500;
            stp.varMax = 5.12; %stp.varCount * 1000; %5.12;
            stp.varMin = -5.12; %-stp.varCount * 1000; %-5.12;
            stp.complexVar = false;
            
            stp.crossoverRate = 0.75;
            stp.mutationScale = 0.15;
            
            stp.populationAmount = 100;
            stp.vectorizedFitness = 0;
            %stp.fitnessFunction = 'DE.Kepler';
            %stp.fitnessFunction = 'SW.Perm';
            %stp.fitnessFunction = 'SW.Beale';
            stp.fitnessFunction = 'DE.Rastringin';
            %stp.fitnessFunction = 'DE.LogEqualExp';
            %stp.fitnessFunction = 'DE.SimpleExpression';
            
            %Criando objeto DE para análise da função Rastringin...
            baseDE = DE(stp);
            
            R = baseDE.run();      
        end
        function z = Rastringin(setup, cromo)
            
            x1 = cromo(:,1); x2 = cromo(:,2);
            z = 20 + x1.^2 + x2.^2 - 10*(cos(2*pi*x1)+cos(2*pi*x2));
            
        end
        function z = LogEqualExp(setup, cromo)
            
            x = cromo(:,1);
            z = abs(exp(x) - log(x));
            
        end
        function z = SimpleExpression(setup, cromo)
            
            x = cromo(:,1);
            z = abs(exp(x) - log(x));
            
            %z = abs(x^x^(x^1/4+1/4) - (1/2)^(2^1/2)^(1-2*(2^1/2)));
            %z = abs(x^x^(x^1/4+1/4) - (1/2)^((2^1/2)^(1-2*(2^1/2))));
            %z = abs(x^x^(x^1/4+1/4) - ((1/2)^(2^1/2))^(1-2*(2^1/2)));
            z = abs(x^(x^(x^(1/4)+1/4)) - (1/2)^(sqrt(2)^(1-2*sqrt(2))));
        end
        
        function fit_value = ppt(setup,cromo)
            v = setup.MinValues  + (setup.MaxValues - setup.MinValues).*abs(cos(cromo));
            v(setup.idx) = setup.MinValues(setup.idx);
%             Lc,Le,Re,Rp,V0,Lpe,Rc,Rpe,C,geometry,h,w,l,Te,n_e,t_final,m_bit;
            Lc = v(1);
            Le = v(2);
            Re = v(3);
            Rp = v(4);
            V0 = v(5);
            Lpe = v(6);
            Rc = v(7);
            Rpe = v(8);
            C = v(9);
            h = v(10);
            w = v(11);
            l = v(12);
            t_final = v(13);
            m_bit = v(14);
            geometry = setup.ppt_geometry;
            

%             display(setup.MinValues(14))
%             display(setup.MaxValues(14))
%             display(cromo(14))
%             display(m_bit)
%             display("debug ppt function")
            tspan = [0 t_final];
            opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep', 1e-8);
            mu = 0.0000012566;
            Te = 1.5;
            n_e = 1e21;
            delta = 0;
            R = 8.314;
            V_critic = 13300;
            gamma = 1;
            Req = Re + Rp + Rc + Rpe;
            
            y0 = [0,0,0,0,m_bit];
            %geometry = 0 (parallel)
            %geometry = 1 (coaxial)
            if ~geometry
%                 display("debug ppt parallel begins");
                [t,y] = ode45(@(t,y) ppt_parallel(t,y,V0,C,Rc,Re,Rpe,Rp,mu,h,w,m_bit,Lc,Le,Lpe,delta,Te,n_e,gamma,R,V_critic), tspan, y0, opts);
%                 display("debug ppt parallel end");
                
%                 Ekinetic = y(:,5).*y(:,3).*y(:,3)/2;
%                 m_t = [0; diff(y(:,5))];
%                     
%                 C_m = gamma*(2/(gamma+1))^(gamma+1/(2*(gamma-1)));
%                 C_f = gamma*(2/(gamma+1))^(gamma/(gamma-1)) + (1+(gamma-1)/2)^(-gamma/(gamma-1));
%                 
%                 Eb = (0.5)*(Lc+Le + mu*h.*y(:,1)/w).*y(:,4).^2;
%                 I_bit = zeros(length(t),1);
%                 for i = 2:length(t)
%                     int = trapz(t(1:i),mu*h.*y(1:i,4).*y(1:i,4)./(2.*w) + m_t(1:i).*sqrt(gamma.*R.*Te).*(C_f/C_m));
%                     I_bit(i) = int;
%                 end
%                 
%                 I = y(:,4);
                V = V0 - y(:,2)/C;
                Ec = C*V.^2/2;
%                 Eohm = zeros(length(t),1);  
%                 for i = 2:length(t)
%                     int = trapz(t(1:i),Req*y(1:i,4).*y(1:i,4));
%                     Eohm(i) = int;
%                 end
                
                try
                    index = find(abs(y(:,1)-l) < 0.001,1);
                    t_prime = t(index);
                    exhaust_v = y(index,3);
                    %app.Isp.Value = app.exhaust_v.Value/app.g.Value;
                    %app.I_bit.Value = I_bit(end);
                    n_t= 100*m_bit.*exhaust_v.*exhaust_v./(2*Ec(1));
%                     display("efficiency calculated")
                    if isempty(n_t)
                        fit_value = 0;
                    else
                        fit_value = n_t;
                    end
                catch
                    fit_value = 0;
                end
            else
                r_in = h;
                r_out = w;
                [t,y] = ode45(@(t,y) ppt_coaxial(t,y,V0,C,Rc,Re,Rpe,Rp,mu,m_bit,Lc,Le,Lpe,delta,Te,n_e,r_out,r_in,gamma,R,V_critic), tspan, y0, opts);
                Ekinetic = y(:,5).*y(:,3).*y(:,3)/2;
                Eb = (0.5)*(Lc+Le+mu*log(r_out/r_in)*y0(1)/(4*pi)).*y(:,4).^2;
                I_bit = zeros(length(t),1);
                for i = 2:length(t)
                    int = trapz(t(1:i),mu*log(r_out./r_in).*y(1:i,4).*y(1:i,4)./(4.*pi));
                    I_bit(i) = int;
                end
                I = y(:,4);
                V = V0 - y(:,2)/C;
                Ec = C*V.^2/2;
                Eohm = zeros(length(t),1);

                for i = 2:length(t)
                    int = trapz(t(1:i),Req*y(1:i,4).*y(1:i,4));
                    Eohm(i) = int;
                end
                
                try
                    index = find(abs(y(:,1)-app.l.Value) < 0.001,1);
                    t_prime = t(index);
                    exhaust_v = y(index,3);
                    Isp = exhaust_v/g;
                    I_bit = I_bit(end);
                    n_t = 100*m_bit.*exhaust_v.*exhaust_v./(2*Ec(1));
                    fit_value = n_t;
                catch
                    fit_value = 0;
                end
           
            end            
                
            function y = ppt_parallel(t,y0,V0,C,Rc,Re,Rpe,Rplasma,mu,h,w,m_bit,Lc,Le,Lpe,delta,Te,n_e,gamma,R,V_critic)
                %Lpe = mu*h*y0(1)/w +mu*delta*h/(2*w);
                L = Lpe + Lc + Le;
                %Rplasma = 8.08*h*sqrt(mu*log(1.24e7*((Te*11604.5250061)^3/n_e)^0.5)/t_pulse)/(w*(Te*11604.5250061)^(0.75));
                Req = Rc + Re + Rpe + Rplasma;
                
                m_dot = y0(4)*y0(4)*mu*h/(4.404*w*V_critic);
                C_m = gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
                C_f = gamma*(2/(gamma+1))^(gamma/(gamma-1)) + (1+(gamma-1)/2)^(-gamma/(gamma-1));
                
                y = zeros(5,1);
                y(1) = y0(3);
                y(2) = y0(4);
                y(3) = ((mu*h*y0(4)*y0(4)/(2*w) - m_dot*y0(1) + m_dot*sqrt(gamma*R*Te)*(C_f/C_m)))/y0(5);
                y(4) = (-y0(2)/C - mu*h*y0(3)*y0(4)/w - Req*y0(4) + V0)/(L);
                y(5) = m_dot;
            end
        
            function y = ppt_coaxial(t,y0,V0,C,Rc,Re,Rpe,Rplasma,mu,m_bit,Lc,Le,Lpe,delta,Te,n_e,r_out,r_in,gamma, R,V_critic)
                %Lpe = mu*log(r_out/r_in)*y0(1)/(4*pi) + mu*delta*log(r_out/r_in)/(4*pi);
                L = Lpe + Lc + Le;
                
                %Rplasma = 2.57*(r_out-r_in)*sqrt(mu*log(1.24e7*((Te*11604.5250061)^3/n_e)^0.5)/t_pulse)/((r_out+r_in)*(Te*11604.5250061)^(0.75));
                Req = Rc + Re + Rpe + Rplasma;
                
                m_dot = y0(4)*y0(4)*mu*log(r_out/r_in)/(4.404*V_critic);
                C_m = gamma*(2/(gamma+1))^(gamma+1/(2*(gamma-1)));
                C_f = gamma*(2/(gamma+1))^(gamma/(gamma-1)) + (1+(gamma-1)/2)^(-gamma/(gamma-1));
                
                y = zeros(5,1);
                y(1) = y0(3);
                y(2) = y0(4);
                y(3) = (mu*log(r_out/r_in)*y0(4)*y0(4) - m_dot*y0(1) + m_dot*sqrt(gamma*R*Te)*(C_f/C_m))/(4*y0(5)*pi);
                y(4) = (-y0(2)/C - mu*log(r_out/r_in)*y0(3)*y0(4)/(4*pi) - Req*y0(4) + V0)/(L);
                y(5) = m_dot;
            end
        end
        
        function z = Kepler(setup, cromo)
            E = 0.2;
            M = 0.5;
            x = cromo(1);
            
            %z = ((x - E * sin(x)) - M);
            z = abs((x - E * sin(x)) - M);
        end
    end
    
end