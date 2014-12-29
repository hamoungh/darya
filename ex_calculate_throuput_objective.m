function [ X_sla ] = ex_calculate_throuput_objective...
    ( session,N1, Z, RT_sla, T_plus_nsteps )
  X_sla = session.optimizer.calculate_throuput_objective(...
                N1, Z', RT_sla', T_plus_nsteps);

end

