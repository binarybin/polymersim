
include("Space.jl")
include("EndMove.jl")
include("SnakeMove.jl")
include("CornerMove.jl")

using SpaceModule
using EndMoveModule
using SnakeMoveModule
using CornerMoveModule

em = EndMove()
cm = CornerMove()
sm = SnakeMove()

nsim = 150
nsumo = 250
beta = 3
total_run = 200000000

space = Space(nsim, 20, nsumo, 20, 250, 250)

initialize(space)

move(move::SnakeMove, space::Space, polyid::Int, polytyp::ASCIIString) =
    SnakeMoveModule.move(move, space, polyid, polytyp)

move(move::EndMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    EndMoveModule.move(move, space, polyid, polytyp)

move(move::CornerMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    CornerMoveModule.move(move, space, polyid, polytyp)

function run()
    running_scale = 10000 # scale for running average for success rates
    
    cm_succ = 0
    cm_succ_list = []
    sm_succ = 0
    sm_succ_list = []
    em_succ = 0
    em_succ_list = []
    
    energylist = []
    energy = 0
    one_percent = div(total_run, 100)
    for run_idx in 1:total_run
        if mod(run_idx, running_scale) == 1
            push!(cm_succ_list, cm_succ/(2*running_scale))
            push!(sm_succ_list, sm_succ/(2*running_scale))
            push!(em_succ_list, em_succ/(2*running_scale))
            sm_succ = 0
            em_succ = 0
            cm_succ = 0
        end
        
        if mod(run_idx, one_percent) == 0
            println(div(run_idx, one_percent), " percent finished")
        end
        
        (success, bond_change) = move(sm, space, rand(1:nsim), "sim")
        if success
            sm_succ += 1
            energy -= bond_change
        end
        
        (success, bond_change) = move(sm, space, rand(1:nsumo), "sumo")
        if success
            sm_succ += 1
            energy -= bond_change
        end
        
        (success, bond_change) = move(em, space, rand(1:nsim), "sim")
        if success
            em_succ += 1
            energy -= bond_change
        end
        
        (success, bond_change) = move(em, space, rand(1:nsumo), "sumo")
        if success
            em_succ += 1
            energy -= bond_change
        end
        
        (success, bond_change) = move(cm, space, rand(1:nsim), "sim")
        if success
            cm_succ += 1
            energy -= bond_change
        end
        
        (success, bond_change) = move(cm, space, rand(1:nsumo), "sumo")
        if success
            cm_succ += 1
            energy -= bond_change
        end
        
        if run_idx == div(total_run, 2)
            sm.beta = beta
            cm.beta = beta
            em.beta = beta
        end
        push!(energylist, energy)
    end
#    println("Energy scale: ", energylist)
#    println("Snake Move Success rate: ", sm_succ_list)
#    println("Corner Move Success rate: ", cm_succ_list)
#    println("End Move Success rate: ", em_succ_list)
    return (cm_succ_list, sm_succ_list, em_succ_list, energylist)
end

@time result = run()
(cm_succ_list, sm_succ_list, em_succ_list, energy_list) = result
#println(space.space)
using PyCall
@pyimport pylab
pylab.figure()
pylab.plot(cm_succ_list)
pylab.plot(sm_succ_list)
pylab.plot(em_succ_list)
pylab.show()

pylab.figure()
pylab.plot(energy_list)
pylab.show()

pylab.figure()
pylab.imshow(space.space)
pylab.show()



