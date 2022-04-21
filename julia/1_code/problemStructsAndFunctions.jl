module problemStructsAndFunctions
export ProblemConditions, JacobianPiece, initialConditionsLinearized, 
    returnMethod, matrixGenerationMethod

# Struct definitions
struct ProblemConditions
    η_k::Array{Float64}
    u_k::Array{Float64}
    z_k::Float64
    v_k::Float64
    P_k::Array{Float64}
    dt ::Float64
end
using SparseArrays
struct JacobianPiece
       full_matrix::SparseArrays.SparseMatrixCSC{Float64, Int64};
          constant::SparseArrays.SparseMatrixCSC{Float64, Int64};
    time_dependent::SparseArrays.SparseMatrixCSC{Float64, Int64};
          residual::Matrix{Float64};
end

include("./initialConditionsLinearized.jl");
include("./matrixGenerationMethod.jl");
include("./returnMethod.jl")

end