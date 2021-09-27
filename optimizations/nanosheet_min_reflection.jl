using simulator.fresnelltools
using simulator.materials.spesifics

using JuMP
import Ipopt


function Reflectance(λ, θ, d)
    n_silica = 1.4
    n_gold = Au
    n_air = 1

    _, S_tm = ThreeLayerSystem(n_silica, n_gold, n_air, λ, θ, d)

    abs((S_tm * [1, 0])[2])^2
end


function Optimize()
    model = Model(Ipopt.Optimizer)
    register(model, :Reflectance, 3, Reflectance; autodiff = true)
    
    @variable(model, 400e-9 <= λ <= 1200e-9, start = 750e-9)
    @variable(model, 15e-9 <= d <= 150e-9, start = 0.8)
    @variable(model, 0 <= θ <= π/2, start = 55e-9)

    @NLobjective(model, Min, Reflectance(λ, θ, d))

    optimize!(model)

    println(model)
    println("λ = ", value(λ))
    println("θ = ", value(θ))
    println("d = ", value(d))
    println("R = ", objective_value(model))
end