classdef InputData
    % InputData  Flight telemetry paths and optimization configuration for 6DOF identification.

    properties
        % Identification flight telemetry
        AccelerationFileName 
        OutputStatesFilesName 
        AccutatorsFileName 

        % Out-of-sample validation telemetry (optional)
        TestAccelerationFileName = '';
        TestOutputStatesFilesName = '';
        TestAccutatorsFileName = '';
        TestOffSet {mustBeNumeric} = 1;
        TestPopulationNumber {mustBeNumeric} = 5000;
        HasTestData logical = false;

        % Output workspace path
        WorkspaceFileName 

        % Optimization parameters and physical bounds
        PopulationNumber {mustBeNumeric} = 10100;
        ParameterCount {mustBeNumeric}
        OffSet {mustBeNumeric} = 4000;
        PopulationSize {mustBeNumeric}
        LowerBounds {mustBeMatrix}
        UpperBounds {mustBeMatrix}
        MaxIteration {mustBeNumeric}
        ValueOfGravitationConstant {mustBeFloat} = 9.81
        ValidationMargin {mustBeNumeric} = 1000
    end

    methods
        function obj = InputData(accelerationFileName, outputStatesFilesName, accutatorsFileName, workspaceFileName,...
            lowerBounds, upperBounds, maxIteration, populationSize, parameterCount)
            obj.AccelerationFileName = accelerationFileName;
            obj.OutputStatesFilesName = outputStatesFilesName;
            obj.AccutatorsFileName = accutatorsFileName;
            obj.WorkspaceFileName = workspaceFileName;
            obj.LowerBounds = lowerBounds;
            obj.UpperBounds = upperBounds;
            obj.MaxIteration = maxIteration;
            obj.PopulationSize = populationSize;
            obj.ParameterCount = parameterCount;
        end

        function obj = SetTestData(obj, testAccelFile, testStatesFile, testActuatorFile, testOffset, testLength)
            % Attach separate telemetry dataset for out-of-sample model validation
            obj.TestAccelerationFileName = testAccelFile;
            obj.TestOutputStatesFilesName = testStatesFile;
            obj.TestAccutatorsFileName = testActuatorFile;
            if nargin >= 5 && ~isempty(testOffset), obj.TestOffSet = testOffset; end
            if nargin >= 6 && ~isempty(testLength), obj.TestPopulationNumber = testLength; end
            obj.HasTestData = true;
        end

        function obj = ChangePopulationParameters(obj, populationNumber, offSet)
            obj.PopulationNumber = populationNumber;
            obj.OffSet = offSet;
        end

        function obj = IncreamentOffSet(obj, increament)
            obj.OffSet = obj.OffSet + increament;
        end
    end
end
