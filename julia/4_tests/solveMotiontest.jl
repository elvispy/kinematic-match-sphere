include("./solveMotion.jl");

# Test 1: Visual tests
TEST_1 = true;
z_k = 1.5;
v_k = 0.0;
N = 50;
R_f = 10.0;
rS = 1.0;
dt = 1e-2;
Tm = 1.0;
mu = 1.68e-2;

X = LinRange(0, R_f*rS, ceil(Int64, R_f * N));
Eta_k = @. exp(-(X^2)/R_f) * cos(X/rS);

# Sphere Falling with initial condition on Eta.
if TEST_1 == true
    solve_motion(N =N, 
        η_k = Eta_k,
        R_f = 10.0,
        rS = rS,
        Tm = Tm,
        z_k= 3.5,
        v_k=0.0,
        plotter = true,
        export_data = false,
        simul_time = 50.0);
end