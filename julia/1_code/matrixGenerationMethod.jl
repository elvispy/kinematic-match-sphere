using .problemStructsAndFunctions: JacobianPiece
using LinearAlgebra: I
using SparseArrays;
# DEBUGGING
# using CSV, Tables
# END DEBUGGING

function matrixGenerationMethod(method::String)::Function
    if occursin(r"linear"i, method)
        return returnLinearizedMatrices;
    else
        return returnFullMatrices;
    end
end

function returnLinearizedMatrices(previous_jacobian_pieces::Array{JacobianPiece}, 
    Ntot::Int64, dr::Float64, Mr::Float64, maximum_time_step::Float64, Vu::Float64)::Array{JacobianPiece}
    first_empty_index = 0;
    while isassigned(previous_jacobian_pieces, first_empty_index + 1) == true && first_empty_index < 10010
        first_empty_index  = first_empty_index + 1;
    end
    @assert first_empty_index < 10000;
    I_M = spdiagm(0 => ones(Float64, Ntot)); #Sparse diagonal matrix
    Z_M = spzeros(Ntot, Ntot);
    eye(N::Int64)::Matrix{Float64} = Matrix(1.0I, N, N); #Redefinition
    A_prime = 2 * eye(Ntot);
    A_prime[1, 1] =  4;
    A_prime[1, 2] = -4;

    for ii = 2:Ntot
        A_prime[ii, ii-1] = -(2*ii-3)/(2*ii-2);
        if ii < Ntot
            A_prime[ii, ii+1] = -(2*ii-1)/(2*ii-2);
        end
    end
    A_prime = (1/ dr^2) * A_prime;

    A_2 = [-I_M; Vu * deepcopy(A_prime) ; spzeros(2, Ntot)];
    A_3 = [ Z_M;           I_M; spzeros(2, Ntot)];

    for new_nb_contact_points = first_empty_index:(first_empty_index + 10)
        rows = zeros(Int64, (2 * Ntot - new_nb_contact_points + 2, ));
        cols = zeros(Int64, (2 * Ntot - new_nb_contact_points + 2, ));

        myA_prime = A_prime;
        myA_prime[1:new_nb_contact_points, :] = zeros(new_nb_contact_points, Ntot);
        myA_prime = sparse(myA_prime);

        myA_1 = [Z_M; myA_prime; spzeros(2, Ntot)];
        residual_1 = myA_1[:, 1:new_nb_contact_points] + 
            [eye(new_nb_contact_points); 
            zeros(2*Ntot - new_nb_contact_points + 2, new_nb_contact_points)];

        myA_1 = myA_1[:, (new_nb_contact_points+1):end];

        rows[1:(Ntot - new_nb_contact_points)] = 1:(Ntot - new_nb_contact_points);
        cols[1:(Ntot - new_nb_contact_points)] = 1:(Ntot - new_nb_contact_points);

        myA_2 = A_2;
        residual_2 = myA_2[:, 1:new_nb_contact_points];
        myA_2 = myA_2[:, (new_nb_contact_points+1):end];

        rows[(Ntot - new_nb_contact_points + 1):(2*Ntot - 2*new_nb_contact_points)] = 
            (Ntot + 1):(2*Ntot - new_nb_contact_points);
        cols[(Ntot - new_nb_contact_points + 1):(2*Ntot - 2*new_nb_contact_points)] = 
            (Ntot -new_nb_contact_points + 1):(2*Ntot - 2*new_nb_contact_points);
        
        rows[(2*Ntot-2*new_nb_contact_points + 1):(2*Ntot - new_nb_contact_points)] = 
            (Ntot - new_nb_contact_points + 1):Ntot;
        cols[(2*Ntot-2*new_nb_contact_points + 1):(2*Ntot - new_nb_contact_points)] = 
            (2*Ntot - new_nb_contact_points + 2) * ones(new_nb_contact_points, 1);

        myA_3 = A_3[:, 1:new_nb_contact_points];
        myA_3[end, :] = -dr^2 * Mr * sparse(int_vector(new_nb_contact_points));

        myA_4 = hcat(sum(residual_1, dims=2), sum(residual_2, dims=2)); 
        myA_4[end-1, end] = -1.0;
        myA_4 = sparse(myA_4);
        myA = hcat(myA_1, myA_2, myA_3, myA_4);
        myA = myA[(new_nb_contact_points+1):end, :];

        rows[end - 1] = 2 * Ntot - new_nb_contact_points + 1;
        cols[end - 1] = 2 * Ntot - new_nb_contact_points + 1;
        rows[end] = 2 * Ntot - new_nb_contact_points + 2;
        cols[end] = 2 * Ntot - new_nb_contact_points + 2;

        time_dependent_part = myA;
        constant_part = sparse(rows, cols, ones((length(cols), )));
        full_jacobian = constant_part + maximum_time_step * time_dependent_part;
        # DEBUGGING
        #CSV.write("./2_pipeline/matrixCP$(new_nb_contact_points)_full.csv",  Tables.table(Matrix(full_jacobian)), writeheader = false);
        #CSV.write("./2_pipeline/matrixCP$(new_nb_contact_points)_const.csv", Tables.table(Matrix(constant_part)), writeheader = false);
        #CSV.write("./2_pipeline/matrixCP$(new_nb_contact_points)_time.csv",  Tables.table(Matrix(time_dependent_part)), writeheader = false);
        # END DEBUGGING
        previous_jacobian_pieces[new_nb_contact_points + 1] = JacobianPiece(
            full_jacobian,
            constant_part,
            time_dependent_part, 
            residual_1
        )
    end
    return previous_jacobian_pieces;   

end

function returnFullMatrices(previous_jacobian_pieces::Array{JacobianPiece}, 
    Ntot::Int64, dr::Float64, Mr::Float64)::Array{JacobianPiece}
    return a;
end

function int_vector(n::Int64)::Vector{Float64}
    if n == 0
        return zeros((0, ));
    elseif n == 1;
        return Float64[π/12];
    else
        S = π * ones(Float64, (n, ));
        S[1] = S[1] / 3;
        for ii = 2:(n-1)
            S[ii] = 2 * (ii - 1) * S[ii];
        end
        S[n] = (3/2 * (n-1) - 1/4) * S[n];
        return S;
    end
end

