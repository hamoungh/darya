function [ output_args ] = ex_solve_lqm( theta_ , N(:,t),  Z )
    [x_q, r_q] =  session.optimizer.solve_lqm( theta_ , N(:,t),  Z)

end

