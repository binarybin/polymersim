module EndMoveModule
export EndMove, move


include("Space.jl")
using SpaceModule

function weight(nbr_bond_inc, beta)
    exp(nbr_bond_inc*beta)
end

type EndMove
    beta :: Float64
    movetype :: ASCIIString
    EndMove() = new(0, "EndMove")
end

function avail_filter(space::Space, point::Tuple{Int64,Int64})
    space.space[point[1], point[2]] == 0
end

function get_possible_moves(move::EndMove, space::Space, poly::Polymer)
    # This functions is almost identical in different spaces except for the avail_filter
    endloc_head, nextloc_head = 1, 2
    endloc_tail, nextloc_tail = length(poly.locs), length(poly.locs)-1
    
    possible_moves_head = neighbor(space, poly.locs[nextloc_head])
    possible_moves_tail = neighbor(space, poly.locs[nextloc_tail])

    # this automatically chooses the filter for that geometry
    avail_f(point) = avail_filter(space, point) 
    
    moves_head = filter(avail_f, possible_moves_head)
    moves_tail = filter(avail_f, possible_moves_tail)

    append!([(endloc_head, point) for point in moves_head] , 
    [(endloc_tail, point) for point in moves_tail] )
    
end

function update_reverse_checking_space(move::EndMove, space::Space, oldpoint, newpoint)
    # handle the reverse checking space
    xnew, ynew = newpoint
    xold, yold = oldpoint
    assert(space.rspace[xnew, ynew, 1] == 0 && space.rspace[xold, yold, 1] != 0)
    space.rspace[xnew, ynew, 1] = space.rspace[xold, yold, 1]
    space.rspace[xnew, ynew, 2] = space.rspace[xold, yold, 2]
    space.rspace[xold, yold, 1], space.rspace[xold, yold, 2] = 0, 0
end

function move(move::EndMove, space::Space, polyid::Int, polytyp::ASCIIString)
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
    
    bond_filter(bpoint) = can_build_bond(space, newpoint, bpoint)
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
            bpoint = choice(bond_choice)
            create_bond(space, newpoint, bpoint)
        end
    end

end

end
