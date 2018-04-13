function [val, index] = addition_chain_genetic(v, fit, popSize, genNum)
%     v = [7, 31];
%     fit = 8;
%     popSize = 50;
%     genNum = 50000;
    v_size = size(v);
    max_number_addition_chain = v(v_size(2));

    generation = 0;
    SelectionProbabilities = rand(popSize, 1);
    
    population = InitialPopulation(popSize, max_number_addition_chain)
    fitness = evaluate(population, v)

    result_check_fitness = check_fitness(fitness, fit);
    while(~result_check_fitness && generation <= genNum)
        parents = select(population, SelectionProbabilities);
        population = reproduce(parents, 50, 0.745);
        fitness = evaluate(population, v);
        result_check_fitness = check_fitness(fitness, fit);
        generation = generation + 1
    end
    
    [val, index] = min(fitness);    
    population(index,:);

end


function out = InitialPopulation(popSize, numberGenes)
    out = randi([0 1],popSize,numberGenes);
%     out = [ ];
end

function parents = select(pop, SelectionProbabilities)
    size_population = size(pop);
    parents = zeros(2,size_population(1),size_population(2));
    
    for i = 1:size_population(1)
        n1 = rand;
        n2 = rand;
        parents(1,i,:) = pop(size_population(1),:);
        parents(2,i,:) = pop(size_population(1),:);
        for j = 1:size_population(1)
            if SelectionProbabilities(j) >= n1
                parents(1,i,:) = pop(j,:);                          
            else if SelectionProbabilities(j) >= n2
                    parents(2,i,:) = pop(j,:);
                end
            end                
        end
    end
end

function result = reproduce(parents, mutationDegree, mutationRate)
    parents_size = size(parents);
    result = zeros(parents_size(2),parents_size(3));    
    for i = 1:parents_size(2)
        result(i,:) = crossover_2(parents(1,i,:), parents(2,i,:));
    end
    result = mutate(result, mutationDegree, mutationRate);
end

function result = crossover_1(parent_1, parent_2)
    parent_size = size(parent_1);
    point = randi([1,parent_size(3)-1],1,1);
    result = cat(3,parent_1(1:point),parent_2(point+1:parent_size(3)));    
end

function result = crossover_2(parent_1, parent_2)
    parent_size = size(parent_1);
    point_1 = randi([1,parent_size(3)-1],1,1);
    point_2 = randi([1,parent_size(3)-1],1,1);
    if (point_1 < point_2)
        result = cat(3,parent_2(1:point_1),parent_1(point_1+1:point_2),parent_2(point_2+1:parent_size(3)));    
    else
        result = cat(3,parent_2(1:point_2),parent_1(point_2+1:point_1),parent_2(point_1+1:parent_size(3)));
    end
end

function pop = mutate(pop, mutationDegree, mutationRate)
    size_population = size(pop);
    if (mutationRate && mutationDegree)
        for i = 1:size_population(1)
            n = rand;
            if n <= mutationRate
                for j = 1:mutationDegree
                    gene = randi([1,size_population(2)],1,1);
                    pop(i,gene) = mod(pop(i,gene) + 1, 2);
                end
            end
        end
    end       
end

function fitness = evaluate(population, v)
    size_population = size(population);
    fitness = zeros(1,size_population(1));
    for i = 1:size_population(1)
        fitness(i) = evaluate_ind(population(i,:), v, size_population(2));
    end
end

function fitness = evaluate_ind(s, v, n)
    largePenalty = inf;
    fitness = 0;
    for i = 1:n
        if (s(i))
            fitness = fitness + 1;
        end
        if ((~isInclude(s,v)) || (s(i) && ~isSumOfPrevous(s, i)))
            fitness = fitness + largePenalty;
            break;
        end
    end
end

function out = check_fitness(fitness, value)
    out = false;
    fitness_size = size(fitness);
    for i = 1:fitness_size(2)
        if (fitness(i) <= value)
            out = true;
            break;
        end
    end
end

function out = isInclude(s,v)
    out = true;
    size_v = size(v);
    for i = 1:size_v(2)
        if (~s(v(i)))
            out = false;
            break;
        end        
    end    
end

function out = isSumOfPrevous(s, i)
    out = false;
    for j = 1:i-1
        for k = 1:i-1
            if ( (j + k == i) && (s(j) == 1) && (s(k) == 1) )
                out = true;
                break;
            end            
        end
    end
end

