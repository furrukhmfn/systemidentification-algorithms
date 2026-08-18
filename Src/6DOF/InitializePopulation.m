function population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds)
% Uniform random initialization across parameter search space
population = lowerBounds + rand(populationSize, parameterCount) .* (upperBounds - lowerBounds);
end