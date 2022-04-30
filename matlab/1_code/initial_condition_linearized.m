function eta = initial_condition_linearized(Fr, dr, Ntot)
    %Lets build A'
    outputname = sprintf("../2_pipeline/%s/out", mfilename);
    file_name = fullfile(outputname, sprintf("Ntot%gdr%gFr%g.mat", Ntot, dr, Fr));
    
    if isfile(file_name) == true
        load(file_name, 'eta');
    else
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

        if isfolder(outputname) == false; mkdir(outputname); end

        save(file_name, 'eta');
    end
    assert(exist('eta', 'var') ~= 0);
end
