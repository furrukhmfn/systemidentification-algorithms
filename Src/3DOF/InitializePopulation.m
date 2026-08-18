function population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds)
% InitializePopulation  Uniform random initialization within parameter box constraints.

population = lowerBounds + rand(populationSize, parameterCount) .* (upperBounds - lowerBounds);

end