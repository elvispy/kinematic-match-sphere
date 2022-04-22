#= # Testing Plots
using Plots; 
function circleShape(h, k, r)
    θ = LinRange(0, 2*π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

#anim = Animation();
plt = plot([0,0.1], Any[rand(2),sin]);
xlabel!("HOla gente");
for x in 0.2:0.1:π
    plot(circleShape(0, x, 1), xlims=[-2, 2], ylims=[-10, 10], size = (1080, 960), position = (5, 5))
    plot!(1:5, 1:5)
    #push!(plt, 2, x, sin(x))
    gui(); sleep(0.1)
end

#= #gif(anim) =#
using FileIO;
function test()
    mech = NaN;
    a = 1;
    while a < 10
        a = a+ 1;
        if a == 2
            mech = 2;
        elseif a == 7
            println(mech)
        end
    end    
    println(pwd())
    save("1_code/test.jld2", "mech", mech);
end


test() =#


include("dummy_test.jl");

errortan = 0
A(3)
println(errortan)