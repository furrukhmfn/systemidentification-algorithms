classdef InputGraphicsData
    % Container for graphics generation and plotting configurations.
    properties
       OutputDataFileName 
       AlgorithmName
       CaseNumber {mustBeNumeric}
    end

    methods
        function obj = InputGraphicsData(outputDataFileName, algorithmName, caseNumber)
            obj.OutputDataFileName = outputDataFileName;
            obj.AlgorithmName = algorithmName;
            obj.CaseNumber = caseNumber;
        end

    end
end
