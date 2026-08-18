function caseStr = caseNaming(num)
% Converts a two-digit case identifier (e.g. 11) to formatted grid label 'Case (1,1)'.
    x = floor(num / 10);
    y = mod(num, 10);
    caseStr = sprintf('Case (%d,%d)', x, y);
end