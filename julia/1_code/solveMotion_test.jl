include("solveMotion.jl")
using .kinematic_match

if true
    
    R_f = 4.0;
    N = 45;
    Ntot = ceil(Int64, R_f * N);
    X = LinRange(0, 5.0, Ntot);
    η_k = -5*exp.(-(X.^2)) ;#.* cos.(2*pi*X);

    solveMotion(plotter = true, vu = 0.1); #, R_f = R_f, N = N);
end