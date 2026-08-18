function runAllAlgorithms()
% Batch benchmark harness executing the full metaheuristic algorithm suite across all 24-parameter test cases.

% Software OpenGL fallback to prevent driver crashes during batch graphics export
try
    opengl software;
catch
end

% Include additional algorithm implementations
algoDir = fullfile(fileparts(mfilename('fullpath')), 'more_algorithms');
if exist(algoDir, 'dir')
    addpath(algoDir);
end

% Algorithm registry: function handle and identifier
algorithms = {
    @ALO24,  'ALO';
    @GWO24,  'GWO';
    @GHOA24, 'GHOA';
    @WOA,    'WOA';
    @SSA,    'SSA';
    @SCA,    'SCA';
    @WCA,    'WCA';
    @MFO,    'MFO';
    @CMAES,  'CMAES';
    @JADE,   'JADE';
    @DFA,    'DFA';
    @SHADE,  'SHADE';
    @LSHADE, 'LSHADE';
    @HOKALMAN24, 'HOKALMAN'
};

nAlgos = size(algorithms, 1);
totalStart = tic;

fprintf('\n========================================================\n');
fprintf('  runAllAlgorithms: Running %d algorithms\n', nAlgos);
fprintf('========================================================\n');

% Run each algorithm sequentially to avoid worker contention
for a = 1:nAlgos
    func = algorithms{a, 1};
    name = algorithms{a, 2};

    fprintf('\n----------------------------------------\n');
    fprintf('  [%d/%d] Starting %s\n', a, nAlgos, name);
    fprintf('----------------------------------------\n');

    algoStart = tic;
    main24(func, name);
    elapsed = toc(algoStart);

    fprintf('  [%d/%d] %s finished in %.1f s\n', a, nAlgos, name, elapsed);
end

totalElapsed = toc(totalStart);
fprintf('\n========================================================\n');
fprintf('  All %d algorithms completed in %.1f s (%.1f min)\n', ...
    nAlgos, totalElapsed, totalElapsed / 60);
fprintf('========================================================\n');

end
