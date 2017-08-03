% ExerciseZLB
%
% Solution of model with alternative policy rules with and without response
% to the spreads, accounting for the ZLB.
%
% Convention for refering to variables in equations:
%   x_t   refers to x{t}
%   x_tF  refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% Required m-files:
%   - symbolic toolbox
%   - LQ package
%   - csolve.m, available in Chris Sims's website
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar  
%
% .........................................................................
%
% Created: April 6, 2010
% Updated: April 24, 2010
% by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
clear all
tic
ttic=toc();
format short g

nsteps = 51;
SolveIterations = 1000;
NumPrecision = 1e-10;

%% Options
% isNatVars = 1;
isNoSpread = 0;
isNoDist = 0;
isGDebt = 0;
isIgnoreZLB = 0;

PersLevel = 90;
% PersLevelAlt = 99;

SigmaRatio = 5;
eta = 50;
dSP = 12;

nMaxZLB = 10;

phi_omega = 0:25:100;

%% Exercise type
ExerciseName = sprintf('ZLB_PersLevel_%02.0f_SigmaRatio_%02.0f_eta_%02.0f_dSP_%.0f',...
    PersLevel,SigmaRatio,eta,dSP);
if exist('PersLevelAlt','var')
    ExerciseName = [ExerciseName,'_PersLevelAlt_',int2str(PersLevelAlt)];
end
if isNoSpread
    ExerciseName = [ExerciseName,'_NoSpread'];
end
if isNoDist
    ExerciseName = [ExerciseName,'_NoDist'];
end
if isGDebt
    ExerciseName = [ExerciseName,'_GDebt'];
end
if isIgnoreZLB
    ExerciseName = [ExerciseName,'_IgnoreZLB'];
end
fprintf('\n****%s****',repmat('*',1,length(ExerciseName)))
fprintf('\n*   %s   *',ExerciseName)
fprintf('\n****%s****\n\n',repmat('*',1,length(ExerciseName)))

%% define parameters
rd_ss = 0.01;
sigmabari = 0.16; %0.16
sigmabar = 1/sigmabari;
phi = 1/0.75;
alpha = 0.66;
omega_y = 0.473;
nu = (omega_y+1)/phi-1;
mu_p = 1.15;
theta = 1+1/(mu_p-1);
spread_ss = (1.02)^(1/4)-1;
omega_ss = ~isNoSpread*spread_ss; 
delta = 0.975;
pi_b = 0.5;
rho_b = 3.2;
elast_omega_b = (eta-1)*spread_ss/(1+spread_ss);
elast_omega_b_TGT = elast_omega_b*4;
isNoRes = 0;
varkappa=0;
Yss = 1;
s_c = 0.7;
sigma_bs = SigmaRatio;
psi = 1;
mu_w_ss = 1;
tau_ss = ~isNoDist*0.2+isNoDist*(1-mu_p*mu_w_ss);
rho_bg = isGDebt*0.1;
Hbar_ss = 1;
Z_ss = 1;
rho_csi = PersLevel/100;
rho_m = (PersLevel~=0)*0.6;
phi_pi = 1.5; %1.5; 2
phi_y = 0.5; %0.5; 1

%% Endogenous variables
y = {'Rd','Pi','Y','lambda_b','lambda_s','b','Delta','K','F','omegatil','b_gy',...
    'Rrdn','Yn','lambda_bn','lambda_sn'};
ny = length(y);
y_ss = []; y_t = []; y_tF = []; y_tL = [];
for j=1:ny
    eval(['syms ',y{j},'_ss ',y{j},'_t ',y{j},'_tF ',y{j},'_tL'])
    y_ss  = [y_ss, eval([y{j},'_ss'])];
    y_t   = [y_t,  eval([y{j},'_t'])];
    y_tF  = [y_tF, eval([y{j},'_tF'])];
    y_tL  = [y_tL, eval([y{j},'_tL'])];
end

%% Identify which variables are to be log-linearized
% yLogIdx = ~ismember(y(:,1),'b');
% yLogIdx = true(1,size(y,1));
yLogIdx = ~ismember(y,{'b_gy'});

%% Exogenous variables
csi = {'xi_i','hZ','hmu_w','htau','hG','hb_g','hHbar','hCbar_b','hCbar_s','hchitil','hXitil'};
ncsi = length(csi);
S = diag([rho_m,rho_csi*ones(1,ncsi-1)]); 
if exist('PersLevelAlt','var')
    S(find(ismember(csi,'hCbar_b')),find(ismember(csi,'hCbar_b')))=PersLevelAlt/100;
    S(find(ismember(csi,'hCbar_s')),find(ismember(csi,'hCbar_s')))=PersLevelAlt/100;
end
csi_t = []; csi_tF = []; csi_tL = []; eps_t = [];
for j=1:ncsi
    eval(['syms ',csi{j},'_t ',csi{j},'_tF ',csi{j},'_tL'])
    csi_t  = [csi_t,  eval([csi{j},'_t'])];
    csi_tF = [csi_tF, eval([csi{j},'_tF'])];
    csi_tL = [csi_tL, eval([csi{j},'_tL'])];
    eval(['syms eps_',csi{j},'_t'])
    eps_t = [eps_t, eval(['eps_',csi{j},'_t'])];
end
csi_ss = zeros(1,ncsi);

%% Lagrange multipliers for F constraints
nF = 8; for jF=1:nF,FLM{jF}=sprintf('FLM%.0f',jF);end
nF = nF+1; FLM{nF}='Upsilon';
FLM_ss = []; FLM_t = []; FLM_tF = [];
for j=1:nF
    eval(['syms ',FLM{j},'_ss ',FLM{j},'_t ',FLM{j},'_tF'])
    FLM_ss  = [FLM_ss, eval([FLM{j},'_ss'])];
    FLM_t   = [FLM_t,  eval([FLM{j},'_t'])];
    FLM_tF  = [FLM_tF, eval([FLM{j},'_tF'])];
end

%% Lagrange multipliers for G constraints
nG = 6; for jG=1:nG,GLM{jG}=sprintf('GLM%.0f',jG);end
GLM_ss = []; GLM_t = []; GLM_tL = [];
for j=1:nG
    eval(['syms ',GLM{j},'_ss ',GLM{j},'_t ',GLM{j},'_tL'])
    GLM_ss  = [GLM_ss, eval([GLM{j},'_ss'])];
    GLM_t   = [GLM_t,  eval([GLM{j},'_t'])];
    GLM_tL  = [GLM_tL, eval([GLM{j},'_tL'])];
end

%% Extended vector z
nz = ny+ncsi;
z = [y,csi];
z_t = [y_t,csi_t];
z_tF = [y_tF,csi_tF];
z_tL = [y_tL,csi_tL];

%% ------------------------------------------------------------------------

%% Auxiliary definitions, used in later sections

%% parameter values and steady state values
Rd_ss = 1+rd_ss;
Pi_ss = 1;
Delta_ss = 1;
Y_ss = Yss;
if omega_ss>0
    beta = (delta+1+omega_ss*(delta+(1-delta)*pi_b)-...
        (((delta+1)+omega_ss*(delta+(1-delta)*pi_b))^2-4*delta*(1+omega_ss))^(1/2))...
        /2/delta/(1+omega_ss)/Rd_ss;
elseif omega_ss==0
	beta = 1/Rd_ss;
else
    error('omega_ss needs to be non-negative')
end
Omega_ss = (1-Rd_ss*beta*(delta+(1-delta)*(1-pi_b)))/Rd_ss/beta/(1-delta)/pi_b;
psi_bs = Omega_ss;
psi_s = psi*(pi_b*psi_bs^(-1/nu)+(1-pi_b))^nu;
psi_b = psi_bs*psi_s;
lambda_s_ss = mu_p*(1+omega_y)*mu_w_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/...
    ((1-tau_ss)*Y_ss*(pi_b*Omega_ss^(1/nu)*psi_b^(-1/nu)+(1-pi_b)*psi_s^(-1/nu))^nu);
lambda_b_ss = Omega_ss*lambda_s_ss;
Lambda_ss = pi_b*lambda_b_ss+(1-pi_b)*lambda_s_ss;
lambdatil_ss = psi*(pi_b*(lambda_b_ss/psi_b)^(1/nu)+(1-pi_b)*(lambda_s_ss/psi_s)^(1/nu))^nu;
Lambdatil_ss = psi^(1/(1+nu))*(pi_b*psi_b^(-1/nu)*lambda_b_ss^((1+nu)/nu)+...
    (1-pi_b)*psi_s^(-1/nu)*lambda_s_ss^((1+nu)/nu))^(nu/(1+nu));
b_ss = rho_b*Y_ss;
sb_ss = (1+pi_b*omega_ss-delta*(1+omega_ss)*Rd_ss)/pi_b/(1-pi_b)*rho_b+...
    (1-delta*Rd_ss)/(1-pi_b)*rho_bg+...
    ((Omega_ss/psi_bs)^(1/nu)-1)/(pi_b*(Omega_ss/psi_bs)^(1/nu)+(1-pi_b))*(1-tau_ss)/mu_p/(1+omega_y);
s_b = s_c+(1-pi_b)*sb_ss;
s_s = s_c-pi_b*sb_ss;
s_bs = s_b/s_s;
s_FI = ~isNoRes*omega_ss/eta*rho_b;
s_g = 1-s_c-(~isNoRes)*omega_ss/eta*rho_b;
sigma_s = sigmabar/(pi_b*s_b*sigma_bs+(1-pi_b)*s_s);
sigma_b = sigma_bs*sigma_s;
sigma = sigmabar/s_c;
Cbar_b_ss = s_b*Y_ss*lambda_b_ss^sigma_b;
Cbar_s_ss = s_s*Y_ss*lambda_s_ss^sigma_s;
ctil_b_ss = s_b*Y_ss;
ctil_s_ss = s_s*Y_ss;
omegatil_ss = 1+omega_ss;
b_gy_ss = rho_bg*Y_ss;
K_ss = Lambda_ss*mu_p*(1+omega_y)*psi*mu_w_ss/lambdatil_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/...
    (1-alpha*beta);
F_ss = Lambda_ss*(1-tau_ss)*Y_ss/(1-alpha*beta);
chitil_ss = isNoRes*omega_ss/(1+varkappa)/b_ss^varkappa;
Xitil_ss = (~isNoRes)*omega_ss/eta/b_ss^(eta-1);
omega_b = elast_omega_b;
omega_chi = 1/(1+omega_ss);
omega_Xi = 1/(1+omega_ss)*eta/rho_b;

%% natural variables
Rrdn_ss = Rd_ss/Pi_ss;
Yn_ss = Y_ss;
lambda_bn_ss = lambda_b_ss;
lambda_sn_ss = lambda_s_ss;
YnL_ss = Y_ss;
YL_ss = Y_ss;

%% a couple more ratios
s_hb = (psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu);
s_hs = (psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu);
s_hbs = s_hb/s_hs;

%% shocks
Z_t = Z_ss*exp(hZ_t);
mu_w_t = mu_w_ss*exp(hmu_w_t);
tau_t = tau_ss+tau_ss*htau_t;
G_t = s_g*Y_ss+Y_ss*hG_t;
b_g_t = rho_bg*Y_ss+Y_ss*hb_g_t;
Hbar_t = Hbar_ss*exp(hHbar_t);
Cbar_b_t = Cbar_b_ss*exp(hCbar_b_t);
Cbar_s_t = Cbar_s_ss*exp(hCbar_s_t);
chitil_t = chitil_ss+1/b_ss^varkappa/(1+varkappa)*hchitil_t;
Xitil_t = Xitil_ss+Y_ss/b_ss^eta*hXitil_t;

%% definitions needed for BW variables
Phi = 1-(theta-1)/theta*(1-tau_ss)/mu_w_ss;
xi = (1-alpha)/alpha*(1-alpha*beta)/(1+omega_y*theta);
kappa = xi*(omega_y+sigmabari);
omega_tau = tau_ss/(1-tau_ss);
q_pi = theta/kappa*(omega_y+sigmabari+Phi*(1-sigmabari));
q_y = omega_y+sigmabari+Phi*(1-sigmabari)-Phi/sigmabar*(1/s_c-1)/(omega_y+sigmabari);
omega_1 = 1/q_y*(omega_y+sigmabari+Phi*(1-sigmabari));
omega_2 = Phi/s_c*sigmabari...
    /((omega_y+sigmabari)^2+Phi*((1-sigmabari)*(omega_y+sigmabari)-(1/s_c-1)*sigmabari));
omega_3 = (1-Phi)...
    /(omega_y+sigmabari+Phi*(1-sigmabari-(1/s_c-1)*sigmabari/(omega_y+sigmabari)));
omega_4 = omega_tau...
    /(omega_y+sigmabari+Phi*(1-sigmabari-(1/s_c-1)*sigmabari/(omega_y+sigmabari)));
q_t = ((1+omega_y)*hZ_t+nu*hHbar_t)/omega_y;
q_tL = ((1+omega_y)*hZ_tL+nu*hHbar_tL)/omega_y;
g_t = pi_b*s_b*hCbar_b_t+(1-pi_b)*s_s*hCbar_s_t+hG_t;
g_tL = pi_b*s_b*hCbar_b_tL+(1-pi_b)*s_s*hCbar_s_tL+hG_tL;
hYnBW_t = (sigmabari*(g_t+hXitil_t)+omega_y*q_t-hmu_w_t-omega_tau*htau_t)/(omega_y+sigmabari);
hYnBW_tL = (sigmabari*(g_tL+hXitil_tL)+omega_y*q_tL-hmu_w_tL-...
    omega_tau*htau_tL)/(omega_y+sigmabari);
hYsBW_t = omega_1*hYnBW_t-omega_2*hG_t+omega_3*hmu_w_t+omega_4*htau_t;
hYsBW_tL = omega_1*hYnBW_tL-omega_2*hG_tL+omega_3*hmu_w_tL+omega_4*htau_tL;
YsBW_t = Y_ss+Y_ss*hYsBW_t;
YsBW_tL = Y_ss+Y_ss*hYsBW_tL;
lambda_y = q_y/q_pi;
phi_u = kappa/(sigmabar*(kappa^2+lambda_y));
phi_x = lambda_y/(sigmabar*(kappa^2+lambda_y));
s_Omega = pi_b*(1-pi_b)*(s_b*sigma_b-s_s*sigma_s)/sigmabar;
phi_Omega = xi*phi_u*s_Omega;

%% auxiliary definitions
Lambda_t = pi_b*lambda_b_t+(1-pi_b)*lambda_s_t;
lambdatil_t = psi*(pi_b*(lambda_b_t/psi_b)^(1/nu)+(1-pi_b)*(lambda_s_t/psi_s)^(1/nu))^nu;
Lambdatil_t = psi^(1/(1+nu))*(pi_b*psi_b^(-1/nu)*lambda_b_t^((1+nu)/nu)+...
    (1-pi_b)*psi_s^(-1/nu)*lambda_s_t^((1+nu)/nu))^(nu/(1+nu));
ctil_b_t = Cbar_b_t*lambda_b_t^(-sigma_b);
ctil_s_t = Cbar_s_t*lambda_s_t^(-sigma_s);
B_t = ctil_b_t-ctil_s_t-((lambda_b_t/psi_b)^(1/nu)-(lambda_s_t/psi_s)^(1/nu))...
    *(lambdatil_t/psi)^(-(1+nu)/nu)*mu_w_t*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)*Delta_t;
lambdatiln_t = psi*(pi_b*(lambda_bn_t/psi_b)^(1/nu)+(1-pi_b)*(lambda_sn_t/psi_s)^(1/nu))^nu;

%% Monetary Policy rules
ZLBRule.NoZLB = Rd_t-(Rd_ss*Pi_t^phi_pi*(Y_t/Y_ss)^(phi_y/4));
ZLBRule.ZLB = Rd_t-1;
ZLBList = fieldnames(ZLBRule);
nZLB = length(ZLBList);

%% Utility function
U = vpa(pi_b*ctil_b_t^(1-1/sigma_b)*Cbar_b_t^(1/sigma_b)/(1-1/sigma_b)...
    +(1-pi_b)*ctil_s_t^(1-1/sigma_s)*Cbar_s_t^(1/sigma_s)/(1-1/sigma_s)...
    -psi/(1+nu)*(lambdatil_t/Lambdatil_t)^(-(1+nu)/nu)*Hbar_t^(-nu)...
    *(Y_t/Z_t)^(1+omega_y)*Delta_t);

%% G constraints
G = vpa([...
        Rd_t*omegatil_t*beta*...
            ((delta+(1-delta)*pi_b)*lambda_b_tF/Pi_tF+(1-delta)*(1-pi_b)*lambda_s_tF/Pi_tF)-lambda_b_t;
        Rd_t*beta*((1-delta)*pi_b*lambda_b_tF/Pi_tF+...
            (delta+(1-delta)*(1-pi_b))*lambda_s_tF/Pi_tF)-lambda_s_t;
        Lambda_t/lambdatil_t*mu_p*(1+omega_y)*psi*mu_w_t*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)+...
            alpha*beta*Pi_tF^(theta*(1+omega_y))*K_tF-K_t;
        Lambda_t*(1-tau_t)*Y_t+alpha*beta*Pi_tF^(theta-1)*F_tF-F_t;
        Rrdn_t*omegatil_ss*beta*((delta+(1-delta)*pi_b)*lambda_bn_tF+(1-delta)*(1-pi_b)*lambda_sn_tF)-lambda_bn_t;
        Rrdn_t*beta*((1-delta)*pi_b*lambda_bn_tF+(delta+(1-delta)*(1-pi_b))*lambda_sn_tF)-lambda_sn_t;
    ]);
nG = length(G);

%% F constraints
F = vpa([...
    pi_b*(1-pi_b)*B_t-pi_b*b_g_t+...
        delta*(b_tL*omegatil_tL+pi_b*b_gy_tL)*Rd_tL/Pi_t-...
        (1+pi_b*(omegatil_t-1))*b_t;
    pi_b*Cbar_b_t*lambda_b_t^(-sigma_b)+(1-pi_b)*Cbar_s_t*lambda_s_t^(-sigma_s)+...
        Xitil_t*b_t^eta+G_t-Y_t;
    alpha*Delta_tL*Pi_t^(theta*(1+omega_y))+...
        (1-alpha)*((1-alpha*Pi_t^(theta-1))/(1-alpha))^(theta*(1+omega_y)/(theta-1))-...
        Delta_t;
    (F_t/K_t)^((theta-1)/(1+omega_y*theta))-(1-alpha*Pi_t^(theta-1))/(1-alpha);
    1+(1+varkappa)*chitil_t*b_t^varkappa+eta*Xitil_t*b_t^(eta-1)-omegatil_t;
    b_g_t-b_gy_t;
    pi_b*Cbar_b_t*lambda_bn_t^(-sigma_b)+(1-pi_b)*Cbar_s_t*lambda_sn_t^(-sigma_s)+...
        Xitil_ss*b_ss^eta+G_t-Yn_t;
    1/lambdatiln_t*mu_p*(1+omega_y)*psi*mu_w_ss*Hbar_t^(-nu)*(Yn_t/Z_t)^(1+omega_y)-(1-tau_ss)*Yn_t;
    Rd_t-1;
    ]);
nF = length(F);

%% ------------------------------------------------------------------------

%% Compute steady state
fprintf('\nSolving for steady state...\n')

%% generate the steady state system
ssSys = vpa(jacobian(U,y_t).'+...
        (FLM_ss*jacobian(F,y_t)).'+beta*(FLM_ss*jacobian(F,y_tL)).'+...
        (GLM_ss*jacobian(G,y_t)).'+beta^(-1)*(GLM_ss*jacobian(G,y_tF)).');
ssSys = [ssSys; G];
ssSys = [ssSys; F];
ssSys(end) = [];
% plug in the steady state values in symbolic form
ssSys = subs(ssSys,[y_t,y_tF,y_tL,csi_t],[y_ss,y_ss,y_ss,csi_ss]);

%% prepare variables
ssSys1 = eval(ssSys(1:ny));
nSys = length(ssSys1);
x = [FLM,GLM];
nx = length(x);

%% generate function for csolve
SolveFileName = sprintf('ssSys%s',ExerciseName);
fid=fopen([SolveFileName,'.m'],'w');
fprintf(fid,'function f=%s(x) \n',SolveFileName);
fprintf(fid,'f = ones(size(x));\n');
fprintf(fid,'for j=1:size(x,2)  \n');
for j=1:nSys
    fprintf(fid,'%s_ss = x(%.0f,j);  \n',x{j},j);
end
for j=1:nSys
    fprintf(fid,['f(',int2str(j),',j) = ',char(ssSys1(j)),';  \n']);
end
fprintf(fid,'end  \n');
fclose(fid);

%% Solve
[x1,rc] = csolve(SolveFileName,zeros(nx,1),[],NumPrecision,1000);
if rc~=0, error(['Solution of steady state system is not normal, rc = ', int2str(rc)]), end
% [(1:nx)' feval(SolverFileName,x1)] % check system solution
delete([SolveFileName,'.m'])
% x1 = round(x1/NumPrecision)*NumPrecision;

%% evaluate variables
for j=1:nx
    eval(sprintf('%s_ss = %.16f;',x{j},x1(j)))
end

%% Check steady state
ssSys = eval(ssSys);
if ~all(abs(ssSys)<1e-6)
    fprintf('\nWARNING: system solution is not precise\n')
    [(1:length(ssSys))' ssSys]
end

%% Present steady state
fprintf('\nSteady state results:')
fprintf('\n=====================\n\n')
for j=1:ny
    fprintf('%15s_ss = %12.6f\n',y{j},eval([y{j},'_ss']))
end
for j=1:nF
    fprintf('%15s_ss = %12.6f\n',FLM{j},eval([FLM{j},'_ss']))
end
for j=1:nG
    fprintf('%15s_ss = %12.6f\n',GLM{j},eval([GLM{j},'_ss']))
end
disp(' ')
y_ss = eval(y_ss);
FLM_ss = eval(FLM_ss);
GLM_ss = eval(GLM_ss);

%% present some ratios
fprintf('\nSome ratios:')
fprintf('\n============\n\n')
RatioList = {...
    'omega_ss','delta','pi_b','beta',...
    'rho_b','s_g','rho_bg','s_c','s_b','s_s','s_bs','s_hb','s_hs','s_hbs',...
    'sigmabar','sigma','sigma_b','sigma_s','sigma_bs',...
    'Cbar_b_ss','Cbar_s_ss','ctil_b_ss','ctil_s_ss',...
    'Omega_ss','psi','psi_b','psi_s','psi_bs',...
    'lambdatil_ss','Lambda_ss','Lambdatil_ss',...
    'phi_pi','phi_y',...
    'Phi','xi','kappa','s_Omega','phi_Omega',...
    'eta','varkappa','chitil_ss','Xitil_ss','omega_chi','omega_Xi',...
    's_FI','elast_omega_b',...
    };
xdisp(RatioList)
disp(' ')

%% Check that steady state satisfies the constraints imposed in model

% check that delta<beta
fprintf('Checking if delta<beta...')
if delta<beta
    fprintf(' Passed\n')
else
    fprintf(' Failed\n')
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

% Check that lambda_b_ss>lambda_s_ss
fprintf('Checking if lambda_b_ss>lambda_s_ss...')
if lambda_b_ss>lambda_s_ss
    fprintf(' Passed\n')
else
    fprintf(' Failed\n')
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

% Check that one of the following holds:
%   beta*(delta+(1-delta)*pi_b)-delta>=0
% or
%   omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta)
fprintf('Checking if beta*(delta+(1-delta)*pi_b)-delta>=0...')
Check1 = (beta*(delta+(1-delta)*pi_b)-delta>=0);
if Check1, fprintf(' Passed\n'), else fprintf(' Failed\n'), end
fprintf('Checking if omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta)...')
Check2 = (omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta));
if Check2, fprintf(' Passed\n\n'), else fprintf(' Failed\n'), end
if ~Check1&&~Check2
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

%% ------------------------------------------------------------------------

%% Solve for FOC
fprintf('Solving for FOC...\n')

%% FOC
FOC = vpa(...
    jacobian(U,y_t)+...
    FLM_t*jacobian(F,y_t)+beta*FLM_tF*subs(jacobian(F,y_tL),[z_t,z_tL],[z_tF,z_t])+...
    GLM_t*jacobian(G,y_t)+beta^(-1)*GLM_tL*subs(jacobian(G,y_tF),[z_t,z_tF],[z_tL,z_t])...
    ).';

%% ------------------------------------------------------------------------

%% Log Linearize everything
fprintf('Log-linearizing...\n')

%% Substitute variables with logs
hy_t = y_t;
hy_tF = y_tF;
hy_tL = y_tL;
hy_ss = y_ss;
for j=1:ny
    if yLogIdx(j)
        hy_t = subs(hy_t,y_t(j),['h' char(y_t(j))]);
        hy_tF = subs(hy_tF,y_tF(j),['h' char(y_tF(j))]);
        hy_tL = subs(hy_tL,y_tL(j),['h' char(y_tL(j))]);
        hy_ss(j) = log(y_ss(j));
        % Plug in transformed variables into the model
        G = subs(G, [y_t(j),y_tF(j)], exp([hy_t(j),hy_tF(j)]));
        F = subs(F,[y_t(j),y_tL(j)],exp([hy_t(j),hy_tL(j)]));
        FOC = subs(FOC,[y_t(j),y_tL(j),y_tF(j)],exp([hy_t(j),hy_tL(j),hy_tF(j)]));
    end
end
hz_t = [hy_t,csi_t];
hz_tF = [hy_tF,csi_tF];
hz_tL = [hy_tL,csi_tL];
hz_ss = [hy_ss,csi_ss];

%% Log-linearize expressions
FnD = {...
    'G_z',G,hz_t;
    'G_zF',G,hz_tF;
	'F_z',F,hz_t;
	'F_zL',F,hz_tL;
	'FOC_z',FOC,hz_t;
	'FOC_zL',FOC,hz_tL;
	'FOC_zF',FOC,hz_tF;
	'FOC_FLM',FOC,FLM_t;
	'FOC_FLMF',FOC,FLM_tF;
	'FOC_GLM',FOC,GLM_t;
	'FOC_GLML',FOC,GLM_tL;
    };
for j=1:size(FnD,1)
    FnDj = vpa(jacobian(FnD{j,2},FnD{j,3}));
    idxSubs = find(FnDj~=0);
    FnDj(idxSubs) = subs(FnDj(idxSubs),...
        [hz_t,hz_tL,hz_tF,FLM_t,FLM_tF,GLM_t,GLM_tL],...
        [hz_ss,hz_ss,hz_ss,FLM_ss,FLM_ss,GLM_ss,GLM_ss]);
    eval([FnD{j,1},' = eval(FnDj);'])
end
clear FnD

%% Find constants
FnD = {...
    'C_G',G;
    'C_F',F;
    'C_FOC',FOC;
    };
for j=1:size(FnD,1)
    FnDj = FnD{j,2};
    idxSubs = find(FnDj~=0);
    FnDj(idxSubs) = subs(FnDj(idxSubs),...
        [hz_t,hz_tL,hz_tF,FLM_t,FLM_tF,GLM_t,GLM_tL],...
        [hz_ss,hz_ss,hz_ss,FLM_ss,FLM_ss,GLM_ss,GLM_ss]);
    eval([FnD{j,1},' = eval(FnDj);'])
end
clear FnD

%% ------------------------------------------------------------------------

%% Optimal policy
fprintf('Solving for Optimal policy state space matrices...\n')

%% some variables needed
Lhz_t = hz_t;
Lhz_tL = hz_tL;
LGLM_t = GLM_t;
LGLM_tL = GLM_tL;
for j=1:nz
    Lhz_t = subs(Lhz_t,hz_t(j),['L' char(hz_t(j))]);
end
for j=1:nG
    LGLM_t = subs(LGLM_t,GLM_t(j),['L' char(GLM_t(j))]);
end
k_t = [hz_t,FLM_t,GLM_t,Lhz_t,LGLM_t];
nk = length(k_t);

%% OptPol, ZLB
G0j = [...
    -FOC_zF,-FOC_FLMF,zeros(ny,nz+2*nG);
    -G_zF,zeros(nG,nz+nF+2*nG);
    -F_z,zeros(nF,nz+nF+2*nG);
    zeros(ncsi,ny),eye(ncsi),zeros(ncsi,nz+nF+2*nG);
    zeros(nz+nG,nz+nF+nG),eye(nz+nG);
    ];
G1j = [...
    FOC_z,FOC_FLM,FOC_GLM,FOC_zL,FOC_GLML;
    G_z,zeros(nG,nz+nF+2*nG);
    F_zL,zeros(nF,nz+nF+2*nG);
    zeros(ncsi,ny),S,zeros(ncsi,nz+nF+2*nG);
    eye(nz),zeros(nz,nz+nF+2*nG);
    zeros(nG,nz+nF),eye(nG),zeros(nG,nz+nG);
    ];
Cj = [...
    C_FOC;
    C_G;
    C_F;
    zeros(ncsi+nz+nG,1)
    ];
G2j = [...
    zeros(ny+nG+nF,ncsi);
    eye(ncsi);
    zeros(nz+nG,ncsi);
    ];
G3j = eye(nk,ny+nG);
cv = find(all(G0j(1:ny+nG,:)==0,2)==1);
G0j(cv,:) = -G1j(cv,:);
G1j(cv,:) = 0;
G3j(:,cv) = [];
cv = find(~all(G3j==0,2));
if ~all(all(G2j(cv,:)==0))
    error('Elements of G2j in forward looking equations are non-zero!')
end
Mat.OptPol.ZLB.G0 = G0j;
Mat.OptPol.ZLB.G1 = G1j;
Mat.OptPol.ZLB.C = Cj;
Mat.OptPol.ZLB.G2 = G2j;
Mat.OptPol.ZLB.G3 = G3j;

%% OptPol, NoZLB
G0j(ny+nG+nF,:) = 0;
idxUpsilon = find(double(jacobian(k_t,Upsilon_t)));
G0j(ny+nG+nF,idxUpsilon) = 1;
Cj(ny+nG+nF) = -Upsilon_ss;
cv = find(all(G0j(1:ny+nG,:)==0,2)==1);
G0j(cv,:) = -G1j(cv,:);
G1j(cv,:) = 0;
G3j(:,cv) = [];
cv = find(~all(G3j==0,2));
if ~all(all(G2j(cv,:)==0))
    error('Elements of G2j in forward looking equations are non-zero!')
end
Mat.OptPol.NoZLB.G0 = G0j;
Mat.OptPol.NoZLB.G1 = G1j;
Mat.OptPol.NoZLB.C = Cj;
Mat.OptPol.NoZLB.G2 = G2j;
Mat.OptPol.NoZLB.G3 = G3j;
[Phi1j,Cj,Phi2j,fmat,fwt,ywt,gev,eu] = gensys(G0j,G1j,Cj,G2j,G3j);
if any(eu~=1),fprintf('WARNING: eu = (%.0f,%.0f)\n',eu),end
REE.OptPol.NoZLB.Phi1 = Phi1j;
REE.OptPol.NoZLB.C = Cj;
REE.OptPol.NoZLB.Phi2 = Phi2j;

%% Spread adjusted policy rules
for jAP=1:length(phi_omega)
    APj = ['TaylorSP',int2str(phi_omega(jAP))];
    fprintf('Solving for %s state space matrices...\n',APj)
    Rulej = sym(['hRrdn_t+phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)-',...
        int2str(phi_omega(jAP)),'/100*homegatil_t+xi_i_t-hRd_t']);
%     RulejG0 = eval(jacobian(Rulej,hz_t));
%     RulejG1 = zeros(1,nz);
%     RulejC = 0;
%     G0j = [...
%         -FOC_zF,-FOC_FLMF,zeros(ny,nz+2*nG);
%         -G_zF,zeros(nG,nz+nF+2*nG);
%         -F_z(1:end-1,:),zeros(nF-1,nz+nF+2*nG);
%         RulejG0,zeros(1,nz+nF+2*nG);
%         zeros(ncsi,ny),eye(ncsi),zeros(ncsi,nz+nF+2*nG);
%         zeros(nz+nG,nz+nF+nG),eye(nz+nG);
%         ];
%     G1j = [...
%         FOC_z,FOC_FLM,FOC_GLM,FOC_zL,FOC_GLML;
%         FOC_z,FOC_FLM,FOC_GLM,FOC_zL,FOC_GLML;
%         G_z,zeros(nG,nz+nF+2*nG);
%         F_zL(1:end-1,:),zeros(nF-1,nz+nF+2*nG);
%         RulejG1,zeros(1,nz+nF+2*nG);
%         zeros(ncsi,ny),S,zeros(ncsi,nz+nF+2*nG);
%         eye(nz),zeros(nz,nz+nF+2*nG);
%         zeros(nG,nz+nF),eye(nG),zeros(nG,nz+nG);
%         ];
%     Cj = [...
%         C_FOC;
%         C_G;
%         C_F(1:end-1,:);
%         RulejC;
%         zeros(ncsi+nz+nG,1)
%         ];
%     G2j = [...
%         zeros(ny+nG+nF,ncsi);
%         eye(ncsi);
%         zeros(nz+nG,ncsi);
%         ];
%     G3j = eye(nk,ny+nG);
%     cv = find(all(G0j(1:ny+nG,:)==0,2)==1);
    RulejG0.NoZLB = eval(jacobian(Rulej,hz_t));
    RulejG1.NoZLB = zeros(1,nz);
    RulejC.NoZLB = 0;
    RulejG0.ZLB = -F_z(end,:);
    RulejG1.ZLB = F_zL(end,:);
    RulejC.ZLB = C_F(end);
    for jZLB=1:nZLB,ZLBj=ZLBList{jZLB};
        G0j = [...
            -G_zF;
            -F_z(1:end-1,:);
            RulejG0.(ZLBj);
            zeros(ncsi,ny),eye(ncsi);
            ];
        G1j = [...
            G_z;
            F_zL(1:end-1,:);
            RulejG1.(ZLBj);
            zeros(ncsi,ny),S;
            ];
        Cj = [...
            C_G;
            C_F(1:end-1,:);
            RulejC.(ZLBj);
            zeros(ncsi,1)
            ];
        G2j = [...
            zeros(nG+nF,ncsi);
            eye(ncsi);
            ];
        G3j = eye(nz,nG);
        cv = find(all(G0j(1:nG,:)==0,2)==1);
        G0j(cv,:) = -G1j(cv,:);
        G1j(cv,:) = 0;
        G3j(:,cv) = [];
        cv = find(~all(G3j==0,2));
        if ~all(all(G2j(cv,:)==0))
            error('Elements of G2j in forward looking equations are non-zero!')
        end
        Mat.(APj).(ZLBj).G0 = G0j;
        Mat.(APj).(ZLBj).G1 = G1j;
        Mat.(APj).(ZLBj).C = Cj;
        Mat.(APj).(ZLBj).G2 = G2j;
        Mat.(APj).(ZLBj).G3 = G3j;
        if strcmp(ZLBj,'NoZLB')
            [Phi1j,Cj,Phi2j,fmat,fwt,ywt,gev,eu] = gensys(G0j,G1j,Cj,G2j,G3j);
            if any(eu~=1),fprintf('WARNING: eu = (%.0f,%.0f)\n',eu),end
            REE.(APj).NoZLB.Phi1 = Phi1j;
            REE.(APj).NoZLB.C = Cj;
            REE.(APj).NoZLB.Phi2 = Phi2j;
        end
    end
% % k_t = [hz_t,FLM_t,GLM_t,Lhz_t,LGLM_t];
%     REE.(APj).Phi1 = [Phi1j,zeros(nz,nk-nz);zeros(nk-nz,nk)];
%     REE.(APj).C = [Cj;zeros(nk-nz,1)];
%     REE.(APj).Phi2 = [Phi2j;zeros(nk-nz,ncsi)];
end

%% List all policies and count them
PolList = fieldnames(REE);
nPol = length(PolList);

% vctoc, return

%% ------------------------------------------------------------------------

%% Generate IRFs
fprintf('\nGenerating IRFs...\n')
% zz = {z{:},FLM{:},...
zz = {z{:},...
    'RdLevel','Rb','Rrd','Rrb','cb','cs','w','wb','ws',...
    'Rrbn','cbn','csn','wn','wbn','wsn',...
    };
nzz = length(zz);
for j=1:nzz
    eval(sprintf('[tf,idx%1$s] = ismember(''%1$s'',zz);',zz{j}))
end
ShockSize = ones(ncsi,1)/100;
ShockSize(ismember(csi,'htau')) = 1/tau_ss/100;
ShockSize(ismember(csi,'xi_i')) = 1/4/100;
ShockSize(ismember(csi,'hb_g')) = 4/100;
ShockSize(ismember(csi,'hchitil')) = dSP/400/omega_chi;
ShockSize(ismember(csi,'hXitil')) = dSP/400/omega_Xi;
ShockSize = diag(ShockSize);
for jS=1:ncsi,Sj = csi{jS};
    fprintf('%s...\n',Sj)
    % OptPol
    for TZLBj=0:nMaxZLB
        % Generate REE
        REEj(nsteps) = REE.OptPol.NoZLB;
        for t=nsteps-1:-1:1
            if t<=TZLBj,ZLBj = 'ZLB';else ZLBj = 'NoZLB';end
            G0j = Mat.OptPol.(ZLBj).G0;
            G1j = Mat.OptPol.(ZLBj).G1;
            Cj  = Mat.OptPol.(ZLBj).C;
            G2j = Mat.OptPol.(ZLBj).G2;
            G3j = Mat.OptPol.(ZLBj).G3;
            cv = find(~all(G3j==0,2));
            Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
            G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
            G1j(cv,:) = 0;
            G0ji = rbinv(G0j);
            REEj(t).Phi1 = G0ji*G1j;
            REEj(t).Phi2 = G0ji*G2j;
            REEj(t).C = G0ji*Cj;
        end
        % Generate IRF
        irfj = REEj(1).C+REEj(1).Phi2*ShockSize(:,jS);
        for t=2:nsteps
            irfj(:,t) = REEj(t).C+REEj(t).Phi1*irfj(:,t-1);
        end
        clear REEj
        % check solution
        CheckZLBj = all((irfj(idxRd,:)>(log(1/Rd_ss)-NumPrecision)));
        CheckUpsilonj = all((irfj(idxUpsilon,:)>-NumPrecision));
        isFoundSequence = CheckZLBj && CheckUpsilonj;
        if isFoundSequence, break, end
    end
    if ~isFoundSequence
        fprintf('WARNING: Did not find sequence for OptCP!\n')
    end
    IRF.(Sj).OptPol = irfj(1:nz+nF,:);
    TZLB.(Sj).OptPol = TZLBj;
    CheckZLB.(Sj).OptPol = CheckZLBj;
    CheckUpsilon.(Sj).OptPol = CheckUpsilonj;

    %% Other policies
    for jPol=2:nPol,Polj=PolList{jPol};
        for TZLBj=0:nMaxZLB
            % Generate REE
            REEj(nsteps) = REE.(Polj).NoZLB;
            for t=nsteps-1:-1:1
                if t<=TZLBj,ZLBj = 'ZLB';else ZLBj = 'NoZLB';end
                G0j = Mat.(Polj).(ZLBj).G0;
                G1j = Mat.(Polj).(ZLBj).G1;
                Cj  = Mat.(Polj).(ZLBj).C;
                G2j = Mat.(Polj).(ZLBj).G2;
                G3j = Mat.(Polj).(ZLBj).G3;
                cv = find(~all(G3j==0,2));
                Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
                G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
                G1j(cv,:) = 0;
                G0ji = rbinv(G0j);
                REEj(t).Phi1 = G0ji*G1j;
                REEj(t).Phi2 = G0ji*G2j;
                REEj(t).C = G0ji*Cj;
%                 if t<=TZLBj
%                     G0j = Mat.OptPol.ZLB.G0;
%                     G1j = Mat.OptPol.ZLB.G1;
%                     Cj  = Mat.OptPol.ZLB.C;
%                     G2j = Mat.OptPol.ZLB.G2;
%                     G3j = Mat.OptPol.ZLB.G3;
%                     cv = find(~all(G3j==0,2));
%                     Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
%                     G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
%                     G1j(cv,:) = 0;
%                     G0ji = rbinv(G0j);
%                     REEj(t).Phi1 = G0ji*G1j;
%                     REEj(t).Phi2 = G0ji*G2j;
%                     REEj(t).C = G0ji*Cj;
%                 else
%                     G0j = Mat.(Polj).G0;
%                     G1j = Mat.(Polj).G1;
%                     Cj  = Mat.(Polj).C;
%                     G2j = Mat.(Polj).G2;
%                     G3j = Mat.(Polj).G3;
%                     cv = find(~all(G3j==0,2));
%                     Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
%                     G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
%                     G1j(cv,:) = 0;
%                     G0ji = rbinv(G0j);
%                     REEj(t).Phi1 = G0ji*G1j;
%                     REEj(t).Phi2 = G0ji*G2j;
%                     REEj(t).C = G0ji*Cj;
%                 end
            end
            % Generate IRF
            irfj = REEj(1).C+REEj(1).Phi2*ShockSize(:,jS);
            for t=2:nsteps
                irfj(:,t) = REEj(t).C+REEj(t).Phi1*irfj(:,t-1);
            end
            clear REEj
            % check solution
            CheckZLBj = all((irfj(idxRd,:)>(log(1/Rd_ss)-NumPrecision)));
            isFoundSequence = CheckZLBj;
            if isFoundSequence, break, end    
        end
        if ~isFoundSequence
            fprintf('WARNING: Did not find sequence for %s!\n',Polj)
        end
        IRF.(Sj).(Polj) = [irfj;NaN(nF,nsteps)];
        TZLB.(Sj).(Polj) = TZLBj;
        CheckZLB.(Sj).(Polj) = CheckZLBj;
    end
end
clear irfj

%% Add more variables
for jS=1:ncsi,Sj = csi{jS};
    for jPol=1:nPol,Polj = PolList{jPol};
        irf = IRF.(Sj).(Polj)(1:ny+ncsi,:);
        % add deposit rate level (annualized pct)
        irf(idxRdLevel,:) = (Rd_ss*exp(irf(idxRd,:))).^4-1;
        % add borrowing rate
        irf(idxRb,:) = irf(idxRd,:)+irf(idxomegatil,:);
        % add real deposit rate
        irf(idxRrd,:) = irf(idxRd,:)-cat(2,irf(idxPi,2:end),NaN);
        % add real lending rate
        irf(idxRrb,:) = irf(idxRb,:)-cat(2,irf(idxPi,2:end),NaN);
        % add cb
        irf(idxcb,:) = irf(idxhCbar_b,:)-sigma_b*irf(idxlambda_b,:);
        % add cs
        irf(idxcs,:) = irf(idxhCbar_s,:)-sigma_s*irf(idxlambda_s,:);
        % add w
        irf(idxw,:) = irf(idxhmu_w,:)-nu*irf(idxhHbar,:)-...
            (pi_b*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_b,:)+...
            (1-pi_b)*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_s,:))+...
            (1+omega_y)*(irf(idxY,:)-irf(idxhZ,:,:))+irf(idxDelta,:);
        % add wb
        irf(idxwb,:) = irf(idxw,:)+...
            (1-pi_b)/nu*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_s,:)-irf(idxlambda_b,:));
        % add ws
        irf(idxws,:) = irf(idxw,:)-...
            pi_b/nu*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_s,:)-irf(idxlambda_b,:));
        % save the irf
        IRF.(Sj).(Polj) = 100*irf;
        clear irf
    end
end

%% display checks and number of periods
for jS=1:ncsi,Sj = csi{jS};
    fprintf('\nShock: %s\n',Sj)
    fprintf('\n  CheckZLB:\n')
    disp(CheckZLB.(Sj))
    fprintf('\n  CheckUpsilon:\n')
    disp(CheckUpsilon.(Sj))
    fprintf('  TZLB:\n')
    disp(TZLB.(Sj))
end

%% ------------------------------------------------------------------------

%% save information
ExerciseName = ['Output_',ExerciseName];
fprintf('\nMAT file: %s\n',ExerciseName)
save(ExerciseName)

%% Elapsed time
% disp(' '), vctoc(ttic), disp(' ')

%% ------------------------------------------------------------------------

