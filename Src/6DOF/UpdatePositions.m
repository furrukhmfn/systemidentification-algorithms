function population = UpdatePositions(population, globalBest, lowerBounds, upperBounds)
% Stochastic step toward global best solution with box-constraint projection

[populationSize, parameterCount] = size(population);
stepSize = rand(populationSize, parameterCount) .* (globalBest - population);
population = population + stepSize;
population = max(lowerBounds, min(upperBounds, population));

end