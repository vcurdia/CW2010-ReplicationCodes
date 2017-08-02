% RunAllModels
%
% Commands the execution of all models under all specifications
%
% .........................................................................
%
% Copyright 2009 by Vasco Curdia and Michael Woodford
% Created: March 24, 2009
% Updated: January 16, 2010

%% ------------------------------------------------------------------------

%% Preamble
clear all
tic

%% setup options

% Models = {'FF','NoFF','RepHH'};
Models = {'FF'};

PersList = {'PersLevel_0','PersLevel_50','PersLevel_90','PersLevel_99'};
SigmaRatioList = {'SigmaRatio_2','SigmaRatio_5'};
etaList = {'eta_1','eta_5','eta_50'};

isNatVars = 1;

NoSpread = 0;
NoDist = 0;
GDebt = 0;

%% Prepare exercise names
nExercise = 0;
for jM=1:length(Models)
    for jPers=1:length(PersList)
        for jSigma=1:length(SigmaRatioList)
            for jeta=1:length(etaList)
                for isNoSpread=NoSpread
                    for isNoDist=NoDist
                        for isGDebt=GDebt
                            nExercise = nExercise+1;
                            Options = {PersList{jPers},SigmaRatioList{jSigma},etaList{jeta}};
                            if isNatVars, Options{end+1} = 'NatVars'; end
                            if isNoSpread, Options{end+1} = 'NoSpread'; end
                            if isNoDist, Options{end+1} = 'NoDist'; end
                            if isGDebt, Options{end+1} = 'GDebt'; end
                            ExerciseCmd{nExercise} = sprintf('IntModel%s',Models{jM});
                            ExerciseOptions{nExercise} = Options;
                        end
                    end
                end
            end
        end
    end
end

%% Run exercises
for jE=1:nExercise
	feval(ExerciseCmd{jE},ExerciseOptions{jE}{:})
end

%% ------------------------------------------------------------------------

%% Elapsed time
% fprintf('\n\n%s\n\n',vctoc)

%% ------------------------------------------------------------------------
