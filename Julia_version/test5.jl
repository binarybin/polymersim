include("App.jl")

beta = 3

function run(app::App, space::Space)
    initialize(space)
    
    em = EndMove()
    cm = CornerMove()
    sm = SnakeMove()
    
    one_percent = div(app.total_run, 100)
    
    for run_idx in 1:app.total_run
        if mod(run_idx, app.running_scale) == 1
            reset_succ(app)
        end
        
        if mod(run_idx, one_percent) == 0
            println(div(run_idx, one_percent), " percent finished")
        end
        
        stochastic_move(app, em, space)
        stochastic_move(app, sm, space)
        stochastic_move(app, cm, space)
        
        if run_idx == div(app.total_run, 2)
            sm.beta = beta
            cm.beta = beta
            em.beta = beta
        end
        push!(app.energy_list, app.energy)
    end
end

NSim = 150
LSim = 20
NSumo = 250
LSumo = 20
Lx = 150
Ly = 150

app = App(10000, 100000, NSim, LSim, NSumo, LSumo, Lx, Ly)
space = Space(app.NSim, app.LSim, app.NSumo, app.LSumo, app.Lx, app.Ly)

@time run(app, space)
@time open("app", "w") do file
    serialize(file, app)
end

@time open("data", "w") do file
    serialize(file, space)
end
