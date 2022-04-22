function eta = initial_condition_linearized(Fr, dr, Ntot)
    %Lets build A'
    A_prime = diag( 1/(dr^2) * ones(Ntot, 1)); %diagonal terms are equal
        A_prime(1, 1) =  (2*1)/(dr^2);
        A_prime(1, 2) = -(2*1)/(dr^2);
        for i = 2:(Ntot-1)
             A_prime(i, i-1)= -(1)/(2*dr^2) * (2*i-3)/(2*i-2);
             A_prime(i, i+1)= -(1)/(2*dr^2) * (2*i-1)/(2*i-2);
        end
        A_prime(Ntot, Ntot - 1) = -(1)/(2*dr^2) * (2*Ntot-3)/(2*Ntot-2);

    R = -Fr * ones(Ntot, 1);
    eta = (A_prime\R)/2;
end
