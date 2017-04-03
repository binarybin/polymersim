
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

nsim = 100
nsumo = 200
beta = 3

space = Space(nsim, 20, nsumo, 20, 2*(nsim+nsumo), 100)

initialize(space)

move(move::SnakeMove, space::Space, polyid::Int, polytyp::ASCIIString) =
    SnakeMoveModule.move(move, space, polyid, polytyp)

move(move::EndMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    EndMoveModule.move(move, space, polyid, polytyp)

move(move::CornerMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    CornerMoveModule.move(move, space, polyid, polytyp)

function run()
    for idx in 1:100000000
        move(sm, space, rand(1:nsim), "sim")
        move(sm, space, rand(1:nsumo), "sumo")
        move(cm, space, rand(1:nsim), "sim")
        move(cm, space, rand(1:nsumo), "sumo")
        move(em, space, rand(1:nsim), "sim")
        move(em, space, rand(1:nsumo), "sumo")
    end
    sm.beta = beta
    cm.beta = beta
    em.beta = beta
    for idx in 1:100000000
        move(sm, space, rand(1:nsim), "sim")
        move(sm, space, rand(1:nsumo), "sumo")
        move(cm, space, rand(1:nsim), "sim")
        move(cm, space, rand(1:nsumo), "sumo")
        move(em, space, rand(1:nsim), "sim")
        move(em, space, rand(1:nsumo), "sumo")
    end
end

@time run()

using PyCall
@pyimport pylab
pylab.imshow(space.space)
pylab.show()



