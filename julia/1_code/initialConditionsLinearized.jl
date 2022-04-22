using LinearAlgebra;
# DEBUGGING
# using CSV, Tables;
# END DEBUGGING

function initialConditionsLinearized(Fr::Float64, dr::Float64, 
    Ntot::Int64)::Vector{Float64}
    A_prime = diagm(0 => (1/dr^2) * ones(Ntot, ));
    A_prime[1, 1] =  (2*1)/(dr^2);
    A_prime[1, 2] = -(2*1)/(dr^2);
    for ii = 2:(Ntot - 1)
        A_prime[ii, ii-1] = -(1)/(2*dr^2) * (2*ii-3)/(2*ii-2);
        A_prime[ii, ii+1] = -(1)/(2*dr^2) * (2*ii-1)/(2*ii-2);
    end
    A_prime[Ntot, Ntot - 1] =  -(1)/(2*dr^2) * (2*Ntot-3)/(2*Ntot-2);
    R = -Fr * ones(Ntot, 1);
    # DEBUGGING
    # CSV.write("./2_pipeline/Ntot$(Ntot)dr$(dr)Fr$(Fr).csv", Tables.table(vec((A_prime\R)/2)), writeheader = false);
    # END DEBUGGING

    return vec((A_prime\R)/2);
end
