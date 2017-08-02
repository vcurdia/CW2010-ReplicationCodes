function W=TaylorSPAltCoefFcn(x,LQmat,LQz_t,REELQ,VarShocks,Policy,phi_pi,phi_y)

for j=1:size(x,2)
    phi_omega = x(1,j);
    syms hPi_t hY_t xi_i_t hRd_t homegatil_t hRrdn_t hYn_t
    if strcmp(Policy,'TaylorNoNatSP')
        RuleAlt = phi_pi*hPi_t+phi_y/4*hY_t-phi_omega*homegatil_t+xi_i_t-hRd_t;
    elseif strcmp(Policy,'TaylorRnSP')
        RuleAlt = phi_pi*hPi_t+phi_y/4*hY_t-phi_omega*homegatil_t+xi_i_t-(hRd_t-hRrdnNoDFF_t);
    elseif strcmp(Policy,'TaylorYnSP')
        RuleAlt = phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)-phi_omega*homegatil_t+xi_i_t-hRd_t;
    elseif strcmp(Policy,'TaylorSP')
        RuleAlt = phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)-phi_omega*homegatil_t+xi_i_t-(hRd_t-hRrdn_t);
    end
    REERule = LQAltRuleLinear(RuleAlt,LQmat,LQmat.S,LQz_t);
    W(1,j) = LQWEval(LQmat,REELQ,REERule,VarShocks);
end
