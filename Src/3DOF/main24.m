function main24(algorithmFunc, algorithmName)
% Executes the 24-parameter longitudinal identification test suite across all cases in parallel.

if nargin < 1 || isempty(algorithmFunc)
    algorithmFunc = @ALO24;
end
if nargin < 2 || isempty(algorithmName)
    s = functions(algorithmFunc);
    algorithmName = s.function;
end

fprintf("main24: Running algorithm %s across all cases\n", algorithmName);

% Initialize local worker pool if not already active
p = gcp('nocreate');
if isempty(p)
    parpool('local');
end
fprintf("%s parallel script started\n", algorithmName);

algFunc = algorithmFunc;
algName = algorithmName;

% Test matrix: [caseNumber, populationSize, offSetIncreament, maxIterations]
caseSpecs = [
    11, 2000, 10000, 20;
    21, 2000, 10000, 20;
    31, 2000, 10000, 20;
    41, 2000, 10000, 20;
    12, 4000, 8000,  20;
    22, 4000, 8000,  20;
    32, 4000, 8000,  20;
    42, 4000, 8000,  20;
    13, 6000, 6000,  20;
    23, 6000, 6000,  20;
    33, 6000, 6000,  20;
    43, 6000, 6000,  20;
    14, 8000, 6000,  20;
    24, 8000, 6000,  20;
    34, 8000, 6000,  20;
    44, 8000, 6000,  20;
    15, 10000, 4000, 20;
    25, 10000, 4000, 20;
    35, 10000, 4000, 20;
    45, 10000, 4000, 20;
];

nCases = size(caseSpecs, 1);

%% Phase 1: Parallel optimization runs across cases
futures = cell(nCases, 1);
for i = 1:nCases
    cs = caseSpecs(i, :);
    cNum = cs(1); cPop = cs(2); cOff = cs(3); cIter = cs(4);
    futures{i} = parfeval(@() caseFileGeneric(cNum, cPop, cOff, cIter, algFunc, algName), 0);
end

% Collect futures and check for runtime errors
for i = 1:nCases
    try
        wait(futures{i});
        if ~isempty(futures{i}.Error)
            cs = caseSpecs(i, :);
            fprintf('  Case %d failed: %s\n', cs(1), futures{i}.Error.message);
        end
    catch ME
        cs = caseSpecs(i, :);
        fprintf('  Case %d failed: %s\n', cs(1), ME.message);
    end
end

fprintf("All %s optimizations finished. Generating graphics...\n", algorithmName);

% Software OpenGL fallback for headless rendering
try
    opengl software;
catch
end

%% Phase 2: Post-processing and trajectory graphics generation
for i = 1:nCases
    cs = caseSpecs(i, :);
    matFile = fullfile('output', algName + "_Case" + num2str(cs(1)) + "_24.mat");
    if exist(matFile, 'file')
        try
            CreateGraphicsLongitudinalDynamics24(matFile, cs(1), algName);
            fprintf('  Graphics for Case %d done.\n', cs(1));
        catch ME
            fprintf('  Graphics for Case %d failed: %s\n', cs(1), ME.message);
        end
    end
end

fprintf("All %s cases finished.\n", algorithmName);

end
