
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

np = 100

space = Space(np, 20, np, 20, 3*np, 100)

initialize(space)

move(move::SnakeMove, space::Space, polyid::Int, polytyp::ASCIIString) =
    SnakeMoveModule.move(move, space, polyid, polytyp)

move(move::EndMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    EndMoveModule.move(move, space, polyid, polytyp)

move(move::CornerMove, space::Space, polyid::Int, polytyp::ASCIIString) = 
    CornerMoveModule.move(move, space, polyid, polytyp)

function run()
    for idx in 1:100000000
        move(sm, space, rand(1:np), rand(["sim", "sumo"]))
        move(cm, space, rand(1:np), rand(["sim", "sumo"]))
        move(em, space, rand(1:np), rand(["sim", "sumo"]))
    end
    sm.beta = 3
    cm.beta = 3
    em.beta = 3
    for idx in 1:100000000
        move(sm, space, rand(1:np), rand(["sim", "sumo"]))
        move(cm, space, rand(1:np), rand(["sim", "sumo"]))
        move(em, space, rand(1:np), rand(["sim", "sumo"]))
    end
end

@time run()

using PyCall
@pyimport pylab
pylab.imshow(space.space)
pylab.show()



