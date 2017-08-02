function W=TaylorBAltCoefFcn(x,LQmat,LQz_t,REELQ,VarShocks,Policy,phi_pi,phi_y)

for j=1:size(x,2)
    phi_b = x(1,j);
    syms hPi_t hY_t xi_i_t hRd_t hb_t hRrdn_t hYn_t
    if strcmp(Policy,'TaylorNoNatB')
        RuleAlt = phi_pi*hPi_t+phi_y/4*hY_t+phi_b/4*hb_t+xi_i_t-hRd_t;
    elseif strcmp(Policy,'TaylorRnB')
        RuleAlt = phi_pi*hPi_t+phi_y/4*hY_t+phi_b/4*hb_t+xi_i_t-(hRd_t-hRrdnNoDFF_t);
    elseif strcmp(Policy,'TaylorYnB')
        RuleAlt = phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)+phi_b/4*hb_t+xi_i_t-hRd_t;
    elseif strcmp(Policy,'TaylorB')
        RuleAlt = phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)+phi_b/4*hb_t+xi_i_t-(hRd_t-hRrdn_t);
    end
    REERule = LQAltRuleLinear(RuleAlt,LQmat,LQmat.S,LQz_t);
    W(1,j) = LQWEval(LQmat,REELQ,REERule,VarShocks);
end
