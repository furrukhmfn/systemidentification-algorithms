function population = UpdatePositions(population, globalBest, lowerBounds, upperBounds)
% UpdatePositions  Random step position update toward global best with box-constraint projection.

[populationSize, parameterCount] = size(population);
stepSize = rand(populationSize, parameterCount) .* (globalBest - population);
population = population + stepSize;
population = max(lowerBounds, min(upperBounds, population));

end