function [inner_iteration_counts,execution_time] = GmresUnpreconditionedDouble(A,b,x0,restart_m,tol)
    disp('start unpreconditioned GMRES in double precision.');

    inner_iteration_counts = 0;
    tic;% starting point for execution time
    while true


        % Arnoldi iterative process.
        % Input: restart_m is the restart parameter
        [Vm,Hm_bar] = Arnoldi(A,b,x0,restart_m);
        inner_iteration_counts = inner_iteration_counts + size(Hm_bar,2)

        r0 = b-A*x0;
        beta = norm(r0);
        [~,real_m] = size(Hm_bar);
        beta_e1 = zeros(real_m+1,1);beta_e1(1) = beta;
        [Rm_bar,gm_bar] = Givens( Hm_bar,beta_e1 );
        % resize Rm_bar and gm_bar
        Rm = Rm_bar(1:real_m,1:real_m);
        gm = gm_bar(1:real_m);

        % y is the shift from the initial x0.
        % ym = inv(Rm)*gm;% solve directly
        ym = BackwardUpperTriangular( Rm,gm );% solve backward

        % x = x0 + V * y
        % real_solution = inv(A)*b;
        % xm = x0+Vm*ym;
        xm = x0+Vm*ym;
        % judge whether to restart
        rm = norm(b-A*xm);
        if rm <= tol% condition of convergence
            break;
        end

        % The specified accuary was not achieved, meet the condition to restarting
        x0 = xm;
    end
    toc;% end for execution time
    execution_time = toc;
end