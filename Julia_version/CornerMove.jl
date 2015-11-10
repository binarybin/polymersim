module CornerMoveModule
export CornerMove, move

include("Space.jl")
using SpaceModule

function weight(nbr_bond_inc, beta)
    exp(nbr_bond_inc*beta)
end

type CornerMove
    beta :: Float64
    movetype :: ASCIIString
    CornerMove() = new(0, "CornerMove")
end

function get_possible_moves(move::CornerMove, space::Space, poly::Polymer)
    # This functions needs to be created for each space and each move
    possible_moves = []
    for pos in 1:length(poly.locs) - 2 # I have to use specific geometry here
        x1, y1 = poly.locs[pos]
        x2, y2 = poly.locs[pos+1]
        x3, y3 = poly.locs[pos+2]
        if x1 != x3 && y1 != y3
            if x1 == x2 # (x1, y1), (x1, y3), (x3, y3) shape
                assert(y1 != y2 && y2 == y3 && x2 != x3)
                if space.space[x3, y1] == 0
                    push!(possible_moves, (pos+1, (x3, y1)))
                end
            else # (x1, y1), (x3, y1), (x3, y3) shape
                assert(y1 == y2 && x2 == x3 && y2 != y3)
                if space.space[x1, y3] == 0
                    push!(possible_moves, (pos+1, (x1,y3)))
                end
            end
        end
    end
    possible_moves
end

function update_reverse_checking_space(move::CornerMove, space::Space, oldpoint, newpoint)
    xnew, ynew = newpoint
    xold, yold = oldpoint
    assert(space.rspace[xnew,ynew,1] == 0 && space.rspace[xold,yold,1] != 0)
    space.rspace[xnew,ynew,1] = space.rspace[xold,yold,1]
    space.rspace[xnew,ynew,2] = space.rspace[xold,yold,2]
    space.rspace[xold,yold,1], space.rspace[xold,yold,2] = 0, 0
end

function move(move::CornerMove, space::Space, polyid::Int, polytyp::ASCIIString)
    if polytyp == "sim"
        poly = space.Sims[polyid]
        sitevalue = 1 # code for free sim
    elseif polytyp == "sumo"
        poly = space.Sumos[polyid]
        sitevalue = -1 # code for free sumo
    else
        error("polymer type undefined")
    end
    
    possible_moves = get_possible_moves(move, space, poly)

    if isempty(possible_moves)
        return
    end
    
    (pointid, newpoint) = possible_moves[rand(1:end)]
    oldpoint = poly.locs[pointid]

    nbr_bond_inc = 0
    
    if in_a_bond(space, oldpoint) # (xold, yold) is in a bond
        nbr_bond_inc -= 1 # will lose a bond
    end
    
    bond_filter(bpoint) = can_build_bond_tentative(space, sitevalue, bpoint)
    bond_choice = filter(bond_filter, neighbor(space, newpoint))

    if !isempty(bond_choice) # (x, y) can build a bond
        nbr_bond_inc += 1 # will create a bond
    end
    
    if rand() < weight(nbr_bond_inc, move.beta) 
        # removes the monomer on oldpoint
        # handles everything except in its polymer.locs and rspace
        safe_remove(space, oldpoint) 

        # create a monomer on newpoint
        # handles everything except in its polymer.locs and rspace
        safe_create(space, newpoint, sitevalue)

        # make changes on polymer.locs
        poly.locs[pointid] = newpoint

        # handle the reverse checking space (This part needs to use specific geometry)
        update_reverse_checking_space(move, space, oldpoint, newpoint)

        # try building new bonds
        if !isempty(bond_choice)
            bpoint = bond_choice[rand(1:end)]
            create_bond(space, newpoint, bpoint)
        end
    end

end
end

