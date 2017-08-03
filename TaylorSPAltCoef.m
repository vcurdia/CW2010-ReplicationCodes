% TaylorSPAltCoef
%
% optimizes a simple rule w.r.t. response to spreads
%
% Options:
%   Model
%     Model name: 'FF'.
%   ShockUsed
%     name of shock used
%   FileNameSuffix
%     File name suffix, defining parametrization different from benchmark.
%
% .........................................................................
%
% Copyright 2009 by Vasco Curdia and Michael Woodford
% Created: April 9, 2009
% Updated: January 16, 2010

%% ------------------------------------------------------------------------

%% preamble
clear all
tic

%% Options
Model = 'FF';
ShockUsed = 'hchitil';
FileNameSuffix = '_PersLevel_90_SigmaRatio_05_eta_05';
Policy = 'TaylorYnSP';
RefPol = 'TaylorYn';
phi_pi = 1.5;
phi_y = 0.5;
Alt_phi_omega = .25:.05:1;
nRulesMax = 5000;
ShowGapRed = 0;
ShowImprovement = 0;
ShowConsEquiv = 1;
ShowBest = 1;
ShowRel = 0;
isTight = 0;
isPrint = 0;
isListIncrease = 0;
NumFactor = 10^5;

%% Parallel options
UseParallel = 0;
% nMaxWorkers = 20;
% nRules = length(Alt_phi_omega);
% nMaxWorkers = min(45,round(nRules/25))
% UseParallel = (nMaxWorkers>=2);

%% ------------------------------------------------------------------------

%% load  environment
FileName = sprintf('Output_%s%s',Model,FileNameSuffix);
load(FileName,'csi','ncsi','LQmat','REE','Rule','LQz_t','W','beta','pi_b',...
    'lambda_b_ss','ctil_b_ss','lambda_s_ss','ctil_s_ss')

%% ConsEquiv factor
ConsEquivFactor = (1-beta)/(pi_b*lambda_b_ss*ctil_b_ss+(1-pi_b)*lambda_s_ss*ctil_s_ss);

%% normalize shocks
ShockSize = eye(ncsi);
for j=1:ncsi
    Vj = diag(ShockSize(:,j));
    idx_Y = ismember(REE.Taylor.z_t,sym('hY_t'));
    Phi1j = REE.Taylor.Phi1;
    Phi2j = REE.Taylor.Phi2;
    VYj = idx_Y*real(lyapcsd(Phi1j,Phi2j*Vj*Phi2j'))*idx_Y';
    ShockSize(j,j) = 0.000292747*ShockSize(j,j)/VYj;
end

%% prepare shock variance matrix
if ismember(ShockUsed,{'','All','all'})
    % ignore xi_i
    VarShocks = diag([0;diag(ShockSize(2:end,2:end))]);
    ShockUsed = 'All';
elseif ismember(ShockUsed,{'Mix'})
    sigma2_mu_w = 1.09;
    sigma2_G = 35.5;
    V_Z = diag(ismember(csi,'hZ'));
    V_mu_w = diag(ismember(csi,'hmu_w'))*sigma2_mu_w;
    V_G = diag(ismember(csi,'hG'))*sigma2_G;
    VarShocks = V_Z+V_mu_w+V_G;
elseif ismember(ShockUsed,{'Mix1'})
    idx_Z = ismember(csi,'hZ');
    idx_G = ismember(csi,'hG');
    idx_chitil = ismember(csi,'hchitil');
    sigma2_Z = 0.5*ShockSize(idx_Z,idx_Z);
    sigma2_G = 0.25*ShockSize(idx_G,idx_G);
    sigma2_chitil = 0.25*ShockSize(idx_chitil,idx_chitil);
    V_Z = diag(idx_Z)*sigma2_Z;
    V_G = diag(idx_G)*sigma2_G;
    V_chitil = diag(idx_chitil)*sigma2_chitil;
    VarShocks = V_Z+V_G+V_chitil;
elseif ismember(ShockUsed,{'Mix2'})
    idx_Z = ismember(csi,'hZ');
    idx_Cbar_s = ismember(csi,'hCbar_s');
    sigma2_Z = 0.5*ShockSize(idx_Z,idx_Z);
    sigma2_Cbar_s = 0.5*ShockSize(idx_Cbar_s,idx_Cbar_s);
    V_Z = diag(idx_Z)*sigma2_Z;
    V_Cbar_s = diag(idx_Cbar_s)*sigma2_Cbar_s;
    VarShocks = V_Z+V_Cbar_s;
elseif ismember(ShockUsed,{'Mix3'})
    idx_Z = ismember(csi,'hZ');
    idx_mu_w = ismember(csi,'hmu_w');
    idx_Hbar = ismember(csi,'hHbar');
    sigma2_Z = 0.5*ShockSize(idx_Z,idx_Z);
    sigma2_mu_w = 0.25*ShockSize(idx_mu_w,idx_mu_w);
    sigma2_Hbar = 0.25*ShockSize(idx_Hbar,idx_Hbar);
    V_Z = diag(idx_Z)*sigma2_Z;
    V_mu_w = diag(idx_mu_w)*sigma2_mu_w;
    V_Hbar = diag(idx_Hbar)*sigma2_Hbar;
    VarShocks = V_Z+V_mu_w+V_Hbar;
elseif ismember(Model,{'FF','NoFF'}) && ismember(ShockUsed,{'hCbar'})
    load(FileName,'s_c','s_b','s_s','pi_b')
    VarShocks = ShockSize*(...
        diag(ismember(csi,'hCbar_b'))*pi_b*s_b/s_c+...
        diag(ismember(csi,'hCbar_s'))*(1-pi_b)*s_s/s_c);
else
    VarShocks = ShockSize*diag(ismember(csi,ShockUsed));
end
if all(VarShocks(:)==0), error('Shock does not match available list!'), end

%% ------------------------------------------------------------------------

%% prepare grid
nCoef_omega = length(Alt_phi_omega);
nRules = nCoef_omega;
fprintf('\nNumber of rules: %.0f\n\n',nRules)
% safeguard agains too long executions
if nRules>nRulesMax
    error('Too many policy alternatives! Choose smaller number of policies or increase NMax.')
end

%% Run through the grid
disp('Evaluating alternative rules...')
if ~UseParallel
%     WAlt = TaylorSPAltCoefFcn(Alt_phi_omega(:)',LQmat,LQz_t,REE.LQ,VarShocks,...
%         Policy,phi_pi,phi_y);
    for j=1:nRules
        WAlt(:,j) = TaylorSPAltCoefFcn(Alt_phi_omega(j)',LQmat,LQz_t,REE.LQ,VarShocks,Policy,phi_pi,phi_y);
    end
else
    for j=1:nMaxWorkers
        idx = ceil(nRules/nMaxWorkers)*(j-1)+(1:ceil(nRules/nMaxWorkers));
        idx = idx(idx<=nRules);
        if isempty(idx)
            break
        end
        JobOptions{j} = {Alt_phi_omega(idx),LQmat,LQz_t,REE.LQ,VarShocks,...
            Policy,phi_pi,phi_y};
    end
    JobOut = vcSubmitQueuedJobs(length(JobOptions),@TaylorSPAltCoefFcn,1,JobOptions,...
        ListPathDependencies,nMaxWorkers);
    WAlt = cat(2,JobOut{:});
    clear JobOut JobOptions
end
WAlt = reshape(WAlt,size(Alt_phi_omega));
WLQ = LQWEval(LQmat,REE.LQ,REE.LQ,VarShocks);
% WTaylor = LQWEval(LQmat,REE.LQ,REE.Taylor,VarShocks);
% WRef = min(WAlt(:));
WRef = LQWEval(LQmat,REE.LQ,REE.(RefPol),VarShocks);
WGapRed = (WAlt-WRef)./(WLQ-WRef);
WAltConsEquiv = ConsEquivFactor.*(WAlt-WRef)*NumFactor;
WLQConsEquiv = ConsEquivFactor.*(WLQ-WRef)*NumFactor;
WRefConsEquiv = ConsEquivFactor.*(WRef-WRef)*NumFactor;

%% list them
[WAltSorted,idxSort] = sort(WAlt(:));
if isListIncrease
    idxUsed = idxSort;
    fprintf('\nWelfare for each policy in increasing order:\n')
else
    idxUsed = 1:nRules;
    fprintf('\nWelfare for each policy:\n')
end
for jRule=1:nRules
    if ShowGapRed
        fprintf('W=%9.6f, WGapRed=%.1f, phi_omega = %.4f\n',...
            WAlt(idxUsed(jRule)),WGapRed(idxUsed(jRule)),Alt_phi_omega(idxUsed(jRule)))
    elseif ShowImprovement
        fprintf('W=%9.6f, W-WRef=%.1f, phi_omega = %.4f\n',...
            WAlt(idxUsed(jRule)),WAlt(idxUsed(jRule))-WRef,Alt_phi_omega(idxUsed(jRule)))
    elseif ShowConsEquiv
        fprintf('W=%9.6f, ConsEquiv=%11.3f, phi_omega = %.4f\n',...
            WAlt(idxUsed(jRule)),WAltConsEquiv(idxUsed(jRule)),Alt_phi_omega(idxUsed(jRule)))
    else
        fprintf('W=%9.6f, phi_omega = %.4f\n',WAlt(idxUsed(jRule)),Alt_phi_omega(idxUsed(jRule)))
    end
end
% if ShowGapRed
%     fprintf('W=%.6f, WGapRed=%.6f, optimal policy\n',WLQ,1)
% elseif ShowImprovement
%     fprintf('W=%.6f, W-WRef=%.6f,optimal policy\n',WLQ,WLQ-WRef)
% else
%     fprintf('W=%.6f, optimal policy\n',WLQ)
% end

%% plot
% figure
Bounds_w = [min(WAlt(:)),WLQ];
Bounds_w(2) = Bounds_w(2)+0.01*(Bounds_w(2)-Bounds_w(1));
Bounds_omega = [min(Alt_phi_omega(:)),max(Alt_phi_omega(:))];
if ShowGapRed
    plot(Alt_phi_omega,WGapRed)
    ylim([0,1])
else
    plot(Alt_phi_omega,WAlt)
    if ShowBest
        hold on
        plot(Alt_phi_omega,ones(1,nCoef_omega)*WLQ,'-r')
        hold off
    end
    ylim(Bounds_w)
end
xlim(Bounds_omega)
if isTight
    axis tight
end
xlabel('\phi_\omega')
title(sprintf('Welfare (%s, \\phi_\\pi=%.2f, \\phi_y=%.2f, %s)',Policy,phi_pi,phi_y,ShockUsed))
if isPrint
    print('-depsc',['Opt_',Policy,'_',ShockUsed,'.eps'])
end

%% identify maximum
disp(' ')
disp('Optimal policy:')
disp('---------------')
fprintf('(%s, phi_pi=%.2f, phi_y=%.2f, %s)\n\n',Policy,phi_pi,phi_y,ShockUsed)
[WOpt,idxOpt] = max(WAlt(:));
if ShowGapRed
    fprintf(' Optimal: W=%9.6f, WGapRed=1.000\n',WLQ)
    fprintf('     Ref: W=%9.6f, WGapRed=0.000\n',WRef)
    fprintf('BestRule: W=%9.6f, WGapRed=%.3f, phi_omega=%.3f\n',...
        WOpt,WGapRed(idxOpt),Alt_phi_omega(idxOpt))
elseif ShowImprovement
    fprintf(' Optimal: W=%9.6f, W-WRef=%.3f\n',WLQ,WLQ-WRef)
    fprintf('     Ref: W=%9.6f, W-WRef=0.000\n',WRef)
    fprintf('BestRule: W=%9.6f, W-WRef=%.3f, phi_omega=%.3f\n',...
        WOpt,WOpt-WRef,Alt_phi_omega(idxOpt))
elseif ShowConsEquiv
    fprintf(' Optimal: W=%9.6f, ConsEquiv=%11.7f\n',WLQ,WLQConsEquiv)
    fprintf('     Ref: W=%9.6f, ConsEquiv=%11.7f\n',WRef,WRefConsEquiv)
    fprintf('BestRule: W=%9.6f, ConsEquiv=%11.7f, phi_omega=%.3f\n',...
        WOpt,WAltConsEquiv(idxOpt),Alt_phi_omega(idxOpt))
else
    fprintf(' Optimal: W=%9.6f\n',WLQ)
    fprintf('BestRule: W=%9.6f, phi_omega=%.3f\n',WOpt,Alt_phi_omega(idxOpt))
end

%% ------------------------------------------------------------------------

%% VD
% if ismember(ShockUsed,{'Mix'})
%     idx_Y = ismember(REE.(RuleName).z_t,sym('hY_t'));
%     Phi1j = REE.(RuleName).Phi1;
%     Phi2j = REE.(RuleName).Phi2;
%     VY_Z = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_Z*Phi2j'))*idx_Y';
%     VY_mu_w = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_mu_w*Phi2j'))*idx_Y';
%     VY_G = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_G*Phi2j'))*idx_Y';
%     VY = VY_Z+VY_mu_w+VY_G;
%     VDY_Z = VY_Z/VY;
%     VDY_mu_w = VY_mu_w/VY;
%     VDY_G = VY_G/VY;
%     fprintf('\nVariance decomposition of output under optimized rule:\n')
%     fprintf('     Z: %6.4f\n',VDY_Z)
%     fprintf('  mu_w: %6.4f\n',VDY_mu_w)
%     fprintf('     G: %6.4f\n',VDY_G)
% elseif  ismember(ShockUsed,{'Mix1'})
%     Phi1j = REE.(RuleName).Phi1;
%     Phi2j = REE.(RuleName).Phi2;
%     VY_Z = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_Z*Phi2j'))*idx_Y';
%     VY_G = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_G*Phi2j'))*idx_Y';
%     VY_chitil = idx_Y*real(lyapcsd(Phi1j,Phi2j*V_chitil*Phi2j'))*idx_Y';
%     VY = VY_Z+VY_G+VY_chitil;
%     VDY_Z = VY_Z/VY;
%     VDY_G = VY_G/VY;
%     VDY_chitil = VY_chitil/VY;
%     fprintf('\nVariance decomposition of output under optimized rule:\n')
%     fprintf('       Z: %6.4f\n',VDY_Z)
%     fprintf('       G: %6.4f\n',VDY_G)
%     fprintf('  chitil: %6.4f\n',VDY_chitil)
% end

%% ------------------------------------------------------------------------

%% Elapsed time
% disp(' '), vctoc, disp(' ')

%% ------------------------------------------------------------------------

