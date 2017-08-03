% function IRFPlotCompareZLB(FileNameSuffix)
%
% IRFPlotCompareZLB
%
% Usage:
%   IRFPlotCompareExercise6(FileNameSuffix)
%
% Arguments:
%   FileNameSuffix [optional] [string]
%     File name suffix, defining parametrization different from benchmark.
%     Default: '_dSP_4_Pers_99'
%
% .........................................................................
%
% Created: April 6, 2010
% Updated: April 24, 2010
% by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
nsteps = 16;
yMaxSlack = []; % 0.001
yMinScale = 1e-2; % 0.01
FigShow = 1;
FigPrint = 0;
SaveFigData = 0;
FigShape = {3,2}; %{3,2},{4,2},{4,3},{3,4},{1,2}
UseRdLevel = 1;
UseBW = 1;

ShowLegend = 1;
XTickStep = 4; %(nsteps-1)/10; %16;

FileNameSuffix = '_PersLevel_90_SigmaRatio_05_eta_50_dSP_12'

Shocks2Plot = {'xi_i','hZ','hmu_w','htau','hG','hb_g','hHbar','hCbar_b','hCbar_s','hchitil','hXitil'};
Shocks2Plot = {'hchitil'};

Pol2Plot = {'OptPol','TaylorSP0','TaylorSP25','TaylorSP50','TaylorSP75','TaylorSP100'};
Pol2PlotPretty = {'Optimal','\phi_\omega = 0','\phi_\omega = 0.25','\phi_\omega = 0.5',...
    '\phi_\omega = 0.75','\phi_\omega = 1'};
% Pol2Plot = {'OptPol'};
% Pol2PlotPretty = {'Optimal'};

nPol = length(Pol2Plot);

% LineStyle = {'-','--','--+','--x','--o','--s'};
% % LineStyle = {'-','--','-.',':','--.','--'};
% MarkerSize = {1,1,4,4,2,3};
% % LineColor = {'b','r','k',[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
% LineColor = {'b','r',[0.87,0.49,0],[0,0.5,0],[0,0.5,0.5],'k'};
% LineWidth = {1.5,1.5,1,1,1,1};

% LineStyle = {'-','--','-','--','-','--'};
% MarkerSize = {1,1,5,5,3,3};
% LineColor = {'b','r',[0,0.8,0],[0,0.6,0],[0,0.4,0],[0,0.2,0],[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
% LineWidth = {2,1.5,1,1,1,1};

% LineStyle = {'-','--','--+','--x','--o','--s'};
% MarkerSize = {1,1,5,5,3,3};
% LineColor = {'b','r','k',[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
if UseBW
    if all(ismember({'LQ','Taylor','FlexTarget'},Pol2Plot))
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
LineWidth = {1,1,1,1,1,1};

%% designate and label the variables to plot and scale
if all([FigShape{:}]==[4,2])
    var_plot = {'Y','Pi','Rd','Rb','Rrd','omegatil','b'};
    var_plot_nat = {'Yn','','','','Rrdn','omegatiln','bn'};
    var_label = {'Y','\pi','i^d','i^b','r^d','\omega','b'};
    scale = [1,4,4,4,4,4,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[3,2])
    if UseRdLevel
        var_plot = {'Y','Pi','RdLevel','omegatil','b'};
        var_label = {'Y','\pi','i^d (level)','\omega','b'};
        scale = [1,4,1,4,1]; % annualize inflation and interest rates
    else
        var_plot = {'Y','Pi','Rd','omegatil','b'};
        var_label = {'Y','\pi','i^d','\omega','b'};
        scale = [1,4,4,4,1]; % annualize inflation and interest rates
    end
%     var_plot = {'Y','Pi','RdLevel','omegatil','Upsilon'};
%     var_label = {'Y','\pi','i^d (level)','\omega','Upsilon'};
%     scale = [1,4,1,4,1]; % annualize inflation and interest rates
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

%% load mat file
load(['Output_ZLB',FileNameSuffix],'zz','nzz','IRF','csi','Rd_ss');

%% Plot IRFs
tid = 0:1:nsteps-1; ntid = length(tid);
nS = length(Shocks2Plot);
nIRF = nPol;
for j=1:nS, Sj = Shocks2Plot{j};
    if ~ismember(Sj,csi), continue, end
    if FigShow
        figure('Name',sprintf('Responses to a shock in %s',Sj))
    else
        figure('Visible','off')
    end
    FigData = cell(nPlots,1);
    for jj=1:nPlots
        hsubp(jj) = subplot(FigShape{:},jj);
        IRF2Plot = NaN(nIRF,ntid);
        for jP=1:nPol
            [tf,var_pos] = ismember(var_plot{jj},zz);
            if tf
                IRF2Plot(jP,:) = scale(jj)*IRF.(Sj).(Pol2Plot{jP})(var_pos,1:nsteps);
            end
        end
        for jP=1:nPol
            plot(tid,IRF2Plot(jP,:),LineStyle{jP},...
                'Color',LineColor{jP},'LineWidth',LineWidth{jP},...
                'MarkerSize',MarkerSize{jP},'MarkerFaceColor',LineColor{jP})
            hold on
        end
        if ismember('\varphi',var_label{jj})
            h=title('');
            set(h,'Interpreter','latex');
            set(h,'String',strrep(['$',var_label{jj},'$'],' (level)$','$ (level)'));
        else
            title(var_label{jj})
        end
        xlim([0 nsteps-1])
        set(gca,'XTick',0:XTickStep:nsteps-1)
%         if strcmp(var_plot{jj},'RdLevel')
%             RefLevel = (Rd_ss^4-1)*100;
%         else
%             RefLevel = 0;
%         end
        RefLevel = 0;
        plot(tid,RefLevel*ones(size(tid)),'k:')
        hold off
        FigData{jj} = [tid;IRF2Plot;RefLevel*ones(size(tid))]';
        yMax = max([RefLevel,max(IRF2Plot(:))]);
        yMin = min([RefLevel,min(IRF2Plot(:))]);
        ySlack = max([0.05*(yMax-yMin),yMaxSlack]);
        ylim([min([yMin-ySlack,RefLevel-yMinScale]) max([yMax+ySlack,RefLevel+yMinScale])])
    end
    if ShowLegend
        if FigShape{1}>1
            hleg = legend(Pol2PlotPretty{:});
            legPos = get(hleg,'Position');
            xUL = get(hsubp((FigShape{1}-1)*FigShape{2}-1),'Position');
            xU = get(hsubp((FigShape{1}-1)*FigShape{2}),'Position');
            xL = get(hsubp(FigShape{1}*FigShape{2}-1),'Position');
            legPos(1) = xL(1)+xU(1)-xUL(1)+(xU(3)-legPos(3))/2;
            legPos(2) = xL(2)+(xL(4)-legPos(4))/2;
            set(hleg,'Position',legPos)
        else
            legend(Pol2PlotPretty{:},'Location','SE')
        end
%     legend('boxoff')
    end
    % convert plot to an eps file
    if FigPrint
        if ismember('RdLevel',var_plot)
            FigName = strrep(sprintf('IRF_ZLB%s_%s_RdLevel.eps',FileNameSuffix,Shocks2Plot{j}),'_h','_');
        else
            FigName = strrep(sprintf('IRF_ZLB%s_%s.eps',FileNameSuffix,Shocks2Plot{j}),'_h','_');
        end
        if UseBW
            print('-deps',FigName)
        else
            print('-depsc2',FigName)
        end
    end
    if SaveFigData
        if ismember('RdLevel',var_plot)
            FigName = strrep(sprintf('FigData_IRF_ZLB%s_%s_RdLevel',FileNameSuffix,Shocks2Plot{j}),'_h','_');
        else
            FigName = strrep(sprintf('FigData_IRF_ZLB%s_%s',FileNameSuffix,Shocks2Plot{j}),'_h','_');
        end
        save(FigName,'FigData')
    end
end

%% ------------------------------------------------------------------------
