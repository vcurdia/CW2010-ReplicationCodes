function setpath

% setpath
%
% Set path to codes needed.

pathBase = '../../Matlab/';
pathList = {...
    'ACR-LQ',...
    'Sims-Gensys',...
    'Sims-Optimize',...
    };
for j=1:length(pathList)
    pathAdd{j} = [pathBase,pathList{j}];
end
addpath(pathAdd{:})
