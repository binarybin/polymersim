include("Space.jl")
include("EndMove.jl")
include("SnakeMove.jl")
include("CornerMove.jl")

type App
    running_scale :: Int
    total_run :: Int
    NSim :: Int
    LSim :: Int
    NSumo :: Int
    LSumo :: Int
    Lx :: Int
    Ly :: Int
    
    cm_succ :: Int
    cm_succ_list :: Array{Float64, 1}
    sm_succ :: Int
    sm_succ_list :: Array{Float64, 1}
    em_succ :: Int
    em_succ_list :: Array{Float64, 1}
    energy :: Int
    energy_list :: Array{Int, 1}
    App(running_scale::Int, total_run::Int, NSim::Int, LSim::Int, NSumo::Int, LSumo::Int, Lx::Int, Ly::Int) = 
    new(running_scale, total_run, NSim, LSim, NSumo, LSumo, Lx, Ly, 0, [], 0, [], 0, [], 0, [])
end


function move(app::App, themove::SnakeMove, space::Space, polyid::Int, polytyp::ASCIIString)
    (success, bond_change) = move(themove, space, polyid, polytyp)
    if success
        app.sm_succ += 1
        app.energy -= bond_change
    end
end

function move(app::App, themove::EndMove, space::Space, polyid::Int, polytyp::ASCIIString)
    (success, bond_change) = move(themove, space, polyid, polytyp)
    if success
        app.em_succ += 1
        app.energy -= bond_change
    end
end

function move(app::App, themove::CornerMove, space::Space, polyid::Int, polytyp::ASCIIString)
    (success, bond_change) = move(themove, space, polyid, polytyp)
    if success
        app.cm_succ += 1
        app.energy -= bond_change
    end
end

function stochastic_move(app::App, themove, space::Space)
    obj = rand(["sim", "sumo"])
    if obj == "sim"
        id = rand(1:space.NSim)
    else
        id = rand(1:space.NSumo)
    end    
    move(app, themove, space, id, obj)
end

function reset_succ(app::App)
    push!(app.cm_succ_list, app.cm_succ/(2*app.running_scale))
    push!(app.sm_succ_list, app.sm_succ/(2*app.running_scale))
    push!(app.em_succ_list, app.em_succ/(2*app.running_scale))
    sm_succ = 0
    em_succ = 0
    cm_succ = 0
end