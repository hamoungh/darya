function [theta_,thetaV_0]=ex_mpc_step(session,X_target_now,T,thetaV_0)
                % this is the perfect information case
                % X_target_now=X_sla(:,t:t+T-1); 
                cap = session.optimizer.cap*ones(1,T);
                % uses: o.num_host , o.num_service, o.num_classes, o.speed_ratio
                [u, thetaV_round, gammaV, X_nfm,  U_nfm, theta,thetaV_orig] = ...
                      session.optimizer.solve_nfm_over_time(...
                       X_target_now, T, thetaV_0, cap, session.r, 'abs');
                theta_=theta(:,:,1); 
                thetaV_0=thetaV_orig(: , :, 2);
                % assume here i solve the real sys to obtain thetaV_0 or x
                % thetaV_0=thetaV_orig(: , :, 2);
end


