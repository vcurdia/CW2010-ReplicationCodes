function IRFPlotCompare(Models,Policies,FileNameSuffix)

% IRFPlotCompare
%
% Generate plots of IRFs in Curdia and Woodford (2010) comparing different
% specifications, or policies (not both)
%
% Usage:
%   IRFPlotCompare(Models)
%   IRFPlotCompare(Models,Policies)
%   IRFPlotCompare(Models,Policies,ExerciseNames)
%
% Arguments:
%   Models
%     List with model names
%   Policies [optional]
%     List with policy names. Default: 'LQ'
%   FileNameSuffix [optional]
%     File name suffix, defining parametrization different from benchmark.
%     Default: '' - sets parametrization to benchmark of '_Pers_Exo'
%
% Specifications available:
%   FF - Financial Frictions model
%   NoFF - model with heterogenous agents but financial frictions shut down
%   RepHH - model reduced to standard Neo-Keynesian model
%
% Policies available:
%   LQ - Optimal policy in timeless perspective
%   Taylor - Simple Taylor rule
%   Nat - natural variables (if on its own or if with more than one rule,
%         use Taylor rule Nat
%
% .........................................................................
%
% Copyright 2009 by Vasco Curdia and Michael Woodford
% Created: March 24, 2009
% Created: January 13, 2010

%% ------------------------------------------------------------------------

%% preamble
nsteps = 16;
yMaxSlack = []; % 0.001
yMinScale = 1e-2; % 0.01
FigShow = 1;
FigPrint = 0;
FigShape = {3,2}; %{3,2},{4,2},{4,3},{3,4},{1,2}
ShowLegend = 1;
SaveFigData = 0;
UseBW = 1;

%% designate and label the variables to plot and scale
if all([FigShape{:}]==[4,2])
    var_plot = {'Y','Pi','Rd','Rb','Rrd','omegatil','b'};
    var_plot_nat = {'Yn','','','','Rrdn','omegatiln','bn'};
    var_label = {'Y','\pi','i^d','i^b','r^d','\omega','b'};
    scale = [1,4,4,4,4,4,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[3,2])
    var_plot = {'Y','Pi','Rd','omegatil','b'};
    var_plot_nat = {'Yn','','','omegatiln','bn'};
    var_label = {'Y','\pi','i^d','\omega','b'};
    scale = [1,4,4,4,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[2,3])
    var_plot = {'Y','Pi','Rd','b','omegatil'};
    var_plot_nat = {'Yn','','','bn','omegatiln'};
    var_label = {'Y','\pi','i^d','b','\omega'};
    scale = [1,4,4,1,4]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[4,3])
    var_plot = {'cs','cb','Y','Rd','Rb','Pi','Rrd','omegatil','b','ws','wb'};
    var_plot_nat = {'csn','cbn','Yn','','','','Rrdn','omegatiln','bn','wsn','wbn'};
    var_label = {'c^s','c^b','Y','i^d','i^b','\pi','r^d','\omega','b','w^s','w^b'};
    scale = [1,1,1,4,4,4,4,4,1,1,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[3,4])
    var_plot = {'cs','Rd','Rrd','ws','cb','Rb','omegatil','wb','Y','Pi','b'};
    var_plot_nat = {'csn','','Rrdn','wsn','cbn','','omegatiln','wbn','Yn','','bn'};
    var_label = {'c^s','i^d','r^d','w^s','c^b','i^b','\omega','w^b','Y','\pi','b'};
    scale = [1,4,4,1,1,4,4,1,1,4,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[1,2])
    var_plot = {'Y','Pi'};
    var_plot_nat = {'Yn',''};
    var_label = {'Output','Inflation'};
    scale = [1,4]; % annualize inflation and interest rates
end
nPlots = length(var_plot);

%% check args
if nargin==1
    fprintf('\nNo policies supplied. LQ optimal policy is assumed.\n')
    Policies = {'LQ'};
end
if nargin<3
    fprintf('\nNo file name suffix provided. Benchmark assumed.\n\n')
    FileNameSuffix = '_PersLevel_90_SigmaRatio_05_eta_05'; 
elseif nargin>3
    error('Invalid number of arguments')
end
fprintf('\nFile name suffix provided: %s',FileNameSuffix)
fprintf('\nPolicies provided: %s\n\n',[Policies{:}])
nM = length(Models);
nP = length(Policies);
% nE = length(ExerciseNames);
if nM>1&&nP>1
    error('FcnError:nMnP','Cannot specify multiple models and policies.\n%s',...
        'Specify either one model and multiple policies, or multiple models and a single policy.')
end
if nP~=2
    PolicyNat = Policies{1};
    if ismember(PolicyNat,{'Nat'})
        PolicyNat = 'Taylor';
    end
else
    PolicyNat = Policies{1};
end

%% Plotting options
if UseBW
    if all(ismember({'LQ','Taylor','FlexTarget'},Policies))
        LineStyle = {'-','--','--+'};
    else
        LineStyle = {'-','--','-.','--+','--x','--o','--s'};
    end
    LineColor = {'k','k','k','k','k','k'};
    MarkerSize = {1,1,5,5,5,3,3};
else
    LineStyle = {'-','--','--+','--x','--o','--s'};
    LineColor = {'b','r','k',[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
    MarkerSize = {1,1,5,5,3,3};
end
LineWidth = 1;

%% load mat file
for j=1:nM
    M.(Models{j}) = load(sprintf('Output_%s%s',Models{j},FileNameSuffix),...
        'zz','csi','IRF');
%     M.(Models{j}).zz = {M.(Models{j}).y(:).name};
    if all(ismember(Policies,{'LQ','PiStab','FlexTarget','Nat'}))
        [tf,idx] = ismember('xi_i',M.(Models{j}).csi);
        M.(Models{j}).csi(idx) = [];
        M.(Models{j}).IRF.LQ(:,:,idx) = [];
        M.(Models{j}).IRF.PiStab(:,:,idx) = [];
        M.(Models{j}).IRF.FlexTarget(:,:,idx) = [];
    end
end

%% Manipulate shocks
if nM>1
    if ismember('FF',Models)
        csi = M.FF.csi(1:end-2);
        M.FF.IRF.(Policies{1})(:,:,end-1:end) = [];
    elseif ismember('NoFF',Models)
        csi = M.NoFF.csi;
    else
        error('Either FF or NoFF should be in the list of models')
    end
else
    csi = M.(Models{1}).csi;
end
ncsi = length(csi);
[tf,idx] = ismember('RepHH',Models);
if tf&&nM>1
    % this section assumes that after eliminating the FF shocks, if
    % appropriate, the last shocks are the hCbar ones...
    csi(ncsi+1) = {'hCbar'};
    jj = 1:nM;
    for j=jj(jj~=idx)
        M.(Models{j}).IRF.(Policies{1})(:,:,ncsi+1)=...
            M.(Models{j}).IRF.(Policies{1})(:,:,ncsi-1)+...
            M.(Models{j}).IRF.(Policies{1})(:,:,ncsi);
    end
    M.RepHH.IRF.(Policies{1})(:,:,ncsi:ncsi+1) = ...
        repmat(M.RepHH.IRF.(Policies{1})(:,:,end),[1,1,2]);
    ncsi = ncsi+1;
    s = load(sprintf('Output_NoFF%s',FileNameSuffix),'s_c','s_b','s_s','pi_b');
    s_c = s.s_c;
    s_b = s.s_b;
    s_s = s.s_s;
    pi_b = s.pi_b;
    M.RepHH.IRF.(Policies{1})(:,:,ncsi-2) = M.RepHH.IRF.(Policies{1})(:,:,ncsi-2)*pi_b*s_b/s_c;
    M.RepHH.IRF.(Policies{1})(:,:,ncsi-1) = M.RepHH.IRF.(Policies{1})(:,:,ncsi-1)*(1-pi_b)*s_s/s_c;
end

%% Prepare legend(s)
if nM==1
    if all(ismember({'LQ','Taylor','TaylorSP25','TaylorSP50','TaylorSP75',...
            'TaylorSP100'},Policies)) || all(ismember({'LQ','TaylorNoNat',...
            'TaylorNoNatSP25','TaylorNoNatSP50','TaylorNoNatSP75','TaylorNoNatSP100'},Policies))
        IRFLegends = {'Optimal','\phi_\omega = 0','\phi_\omega = 0.25',...
            '\phi_\omega = 0.50','\phi_\omega = 0.75','\phi_\omega = 1'};
    elseif all(ismember({'LQ','TaylorBm100','TaylorBm50','Taylor','TaylorBp50',...
            'TaylorBp100'},Policies)) || all(ismember({'LQ','TaylorNoNatBm100',...
            'TaylorNoNatBm50','TaylorNoNat','TaylorNoNatBp50','TaylorNoNatBp100'},Policies))
        IRFLegends = {'Optimal','\phi_b = -1','\phi_b = -0.5',...
            '\phi_b = 0','\phi_b = 0.5','\phi_b = 1'};
    elseif all(ismember({'LQ','TaylorBm50','Taylor','TaylorBp50'},Policies))...
            || all(ismember({'LQ','TaylorNoNatBm50','TaylorNoNat','TaylorNoNatBp50'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.5','\phi_b = 0','\phi_b = 0.5'};
    elseif all(ismember({'LQ','TaylorBm1','TaylorBm05','Taylor','TaylorBp05','TaylorBp1'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.01','\phi_b = -0.005','\phi_b = 0','\phi_b = 0.005','\phi_b = 0.01'};
    elseif all(ismember({'LQ','TaylorBm10','TaylorBm5','Taylor','TaylorBp5','TaylorBp10'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.1','\phi_b = -0.05','\phi_b = 0','\phi_b = 0.05','\phi_b = 0.1'};
    elseif all(ismember({'LQ','TaylorBm20','TaylorBm10','Taylor','TaylorBp10','TaylorBp20'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.2','\phi_b = -0.1','\phi_b = 0','\phi_b = 0.1','\phi_b = 0.2'};
    elseif all(ismember({'LQ','TaylorBm30','TaylorBm15','Taylor','TaylorBp15','TaylorBp30'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.3','\phi_b = -0.15','\phi_b = 0','\phi_b = 0.15','\phi_b = 0.3'};
    elseif all(ismember({'LQ','TaylorBm40','TaylorBm20','Taylor','TaylorBp20','TaylorBp40'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.4','\phi_b = -0.2','\phi_b = 0','\phi_b = 0.2','\phi_b = 0.4'};
    elseif all(ismember({'LQ','TaylorBm50','TaylorBm25','Taylor','TaylorBp25','TaylorBp50'},Policies))
        IRFLegends = {'Optimal','\phi_b = -0.5','\phi_b = -0.25','\phi_b = 0','\phi_b = 0.25','\phi_b = 0.5'};
    elseif all(ismember({'LQ','Taylor','TaylorSP300','TaylorSP600','TaylorSP900','TaylorSP1200'},Policies))
        IRFLegends = {'Optimal','\phi_\omega = 0','\phi_\omega = 3','\phi_\omega = 6',...
            '\phi_\omega = 9','\phi_\omega = 12'};
    elseif all(ismember({'LQ','Taylor','TaylorNoNat'},Policies))
        IRFLegends = {'Optimal','Taylor','Taylor (basic)'};
    else
        for jP=1:nP
            jPolicy = Policies{jP};
            jPolicy = strrep(jPolicy,'LQ','Optimal');
    %         jPolicy = strrep(jPolicy,'Yn','');
            if ismember('SP',jPolicy)
                jPolicy = strrep(jPolicy,'SP','+');
            elseif ismember('Bm',jPolicy)
                jPolicy = strrep(jPolicy,'Bm','-');
            elseif ismember('Bp',jPolicy)
                jPolicy = strrep(jPolicy,'Bp','+');
            elseif ismember('New',jPolicy)
                jPolicy = strrep(jPolicy,'New','New+');
            end
            IRFLegends{1+jP} = jPolicy;
        end
    end
elseif nP==1
    for jM=1:nM
        IRFLegends{1+jM} = Models{jM};
    end
end

%% Plot IRFs
tid = 0:1:nsteps; ntid = length(tid);
nIRF = nM*nP;
for j=1:ncsi
    if FigShow
        figure('Name',sprintf('Responses to a shock in %s',csi{j}))
    else
        figure('Visible','off')
    end
    FigData = cell(nPlots,1);
    for jj=1:nPlots
        hsubp(jj) = subplot(FigShape{:},jj);
        IRF = NaN(nIRF,ntid);
        for jM=1:nM
%             [tf,var_pos] = ismember(var_plot{jj},M.(Models{jM}).zz);
%             if tf
%                 for jP=1:nP
%                     IRF((jM-1)+(jP-1)+1,:) = scale(jj)*...
%                         M.(Models{jM}).IRF.(Policies{jP})(var_pos,1:nsteps+1,j); 
%                 end
%             end
            for jP=1:nP
                if ~ismember(Policies{jP},{'Nat'})
                    [tf,var_pos] = ismember(var_plot{jj},M.(Models{jM}).zz);
                    if tf
                        IRF((jM-1)+(jP-1)+1,:) = scale(jj)*...
                            M.(Models{jM}).IRF.(Policies{jP})(var_pos,1:nsteps+1,j);
                    end
                else
                    [tf,var_pos] = ismember(var_plot_nat{jj},M.(Models{jM}).zz);
                    if tf
                        IRF((jM-1)+(jP-1)+1,:) = scale(jj)*...
                            M.(Models{jM}).IRF.(PolicyNat)(var_pos,1:nsteps+1,j);
                    end
                end
            end
        end
        FigData{jj} = IRF';
        for jM=1:nM
            for jP=1:nP
                jR = (jM-1)+(jP-1)+1;
                plot(tid,IRF(jR,:),LineStyle{jR},...
                    'Color',LineColor{jR},'LineWidth',LineWidth,...
                    'MarkerSize',MarkerSize{jR},'MarkerFaceColor',LineColor{jR})
                hold on
            end
        end
        title(var_label{jj})
        plot(tid,zeros(size(tid)),'k:')
        xlim([0 nsteps])
        set(gca,'XTick',0:4:nsteps)
        yMax = max([0,max(IRF(:))]);
        yMin = min([0,min(IRF(:))]);
        ySlack = max([0.05*(yMax-yMin),yMaxSlack]);
        ylim([min([yMin-ySlack,-yMinScale]) max([yMax+ySlack,yMinScale])])
    end
    if ShowLegend
        if FigShape{1}>1
            hleg = legend(IRFLegends{:});
            legPos = get(hleg,'Position');
            xUL = get(hsubp((FigShape{1}-1)*FigShape{2}-1),'Position');
            xU = get(hsubp((FigShape{1}-1)*FigShape{2}),'Position');
            xL = get(hsubp(FigShape{1}*FigShape{2}-1),'Position');
            legPos(1) = xL(1)+xU(1)-xUL(1)+(xU(3)-legPos(3))/2;
            legPos(2) = xL(2)+(xL(4)-legPos(4))/2;
            set(hleg,'Position',legPos)
        else
            legend(IRFLegends{:},'Location','SE')
        end
%     legend('boxoff')
    end
    % convert plot to an eps file
    if FigPrint
        FigName = sprintf('IRF_%s_%s_%s.eps',[Models{:}],[Policies{:}],csi{j});
        FigName = strrep(FigName,'Taylor','T');
        FigName = strrep(FigName,'FlexTarget','FT');
        FigName = strrep(FigName,'PiStab','Pi');
        FigName = strrep(FigName,'NoNat','NN');
        if UseBW
            print('-deps',FigName)
        else
            print('-depsc2',FigName)
        end
    end
    if SaveFigData
        FigName = sprintf('IRF_%s_%s_%s',[Models{:}],[Policies{:}],csi{j});
        FigName = strrep(FigName,'Taylor','T');
        FigName = strrep(FigName,'FlexTarget','FT');
        FigName = strrep(FigName,'PiStab','Pi');
        FigName = strrep(FigName,'NoNat','NN');
        save(FigName,'FigData')
    end
end

%% ------------------------------------------------------------------------
