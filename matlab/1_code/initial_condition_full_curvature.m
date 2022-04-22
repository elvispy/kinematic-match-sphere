function eta = initial_condition_full_curvature(Fr, dr, Ntot)
    
    
    %% Solving the system
    NEWTON_METHOD_TRESHOLD = 1e-10;
    Eta_k = zeros(Ntot, 1);
    W = Eta_k; % we cut the first equations
    gravityInfluence = Fr * ones(Ntot, 1); 
    %independent_term = [Eta_k((newCPoints+1):end); u_k; z_k;v_k];
    extract = @(A) A((Ntot + 1):(2*Ntot), 1:Ntot);
    subJacobian = @(Eta_k) extract(JacobianCurvature(Eta_k, 0, dr, 0));
    extract2 = @(V) V((Ntot + 1):(2*Ntot));
    subCurvatures = @(Eta_k) extract2(curvature(Eta_k, 0, dr, 0));
    next_iteration = W - subJacobian(Eta_k)\(subCurvatures(Eta_k) - gravityInfluence);
    
    while norm(next_iteration - W) >= NEWTON_METHOD_TRESHOLD
        W = next_iteration;
        next_iteration = W -  subJacobian(W(1:Ntot))  ...
            \(subCurvatures(W(1:Ntot)) - gravityInfluence );
    end
    eta = next_iteration(1:Ntot);
end
