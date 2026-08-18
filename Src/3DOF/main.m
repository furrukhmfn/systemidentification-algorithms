function main()
% Runs the 15-parameter ALO optimization test suite across all case configurations in parallel.

parpool('local');  
disp("script started");

% Dispatch parallel asynchronous evaluations
futures(1) = parfeval(@case11, 0);  
futures(2) = parfeval(@case21, 0);
futures(3) = parfeval(@case31, 0);
futures(4) = parfeval(@case41, 0);

futures(5) = parfeval(@case12, 0);
futures(6) = parfeval(@case22, 0);
futures(7) = parfeval(@case32, 0);
futures(8) = parfeval(@case42, 0);

futures(9) = parfeval(@case13, 0);
futures(10) = parfeval(@case23, 0);
futures(11) = parfeval(@case33, 0);
futures(12) = parfeval(@case43, 0);

futures(13) = parfeval(@case14, 0);
futures(14) = parfeval(@case24, 0);
futures(15) = parfeval(@case34, 0);
futures(16) = parfeval(@case44, 0);

futures(17) = parfeval(@case15, 0);
futures(18) = parfeval(@case25, 0);
futures(19) = parfeval(@case35, 0);
futures(20) = parfeval(@case45, 0);

% Wait for pool completion
wait(futures);
disp("All scripts finished.");

%% Case configurations: Population = 2000
function case11()
    caseFile(11,2000,10000,100);
end
function case21()
    caseFile(21,2000,10000,250);
end
function case31()
    caseFile(31,2000,10000,500);
end
function case41()
    caseFile(41,2000,10000,1000);
end

%% Case configurations: Population = 4000
function case12()
    caseFile(12,4000,8000,100);
end
function case22()
    caseFile(22,4000,8000,250);
end
function case32()
    caseFile(32,4000,8000,500);
end
function case42()
    caseFile(42,4000,8000,1000);
end

%% Case configurations: Population = 6000
function case13()
    caseFile(13,6000,6000,100);
end
function case23()
    caseFile(23,6000,6000,250);
end
function case33()
    caseFile(33,6000,6000,500);
end
function case43()
    caseFile(43,6000,6000,1000);
end

%% Case configurations: Population = 8000
function case14()
    caseFile(14,8000,6000,100);
end
function case24()
    caseFile(24,8000,6000,250);
end
function case34()
    caseFile(34,8000,6000,500);
end
function case44()
    caseFile(44,8000,6000,1000);
end

%% Case configurations: Population = 10000
function case15()
    caseFile(15,10000,4000,100);
end
function case25()
    caseFile(25,10000,4000,250);
end
function case35()
    caseFile(35,10000,4000,500);
end
function case45()
    caseFile(45,10000,4000,1000);
end

end