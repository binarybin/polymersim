include("Space.jl")

function weight(nbr_bond_inc, beta)
    exp(nbr_bond_inc*beta)
end

type SnakeMove
    beta :: Float64
    movetype :: ASCIIString
    SnakeMove() = new(0, "SnakeMove")
end

function avail_filter(space::Space, point::Tuple{Int64,Int64})
    space.space[point[1], point[2]] == 0
end

function get_possible_moves(move::SnakeMove, space::Space, poly::Polymer)
    # This functions is almost identical in different spaces except for the avail_filter
    endloc_head = 1
    endloc_tail = length(poly.locs)
    
    possible_moves_head = neighbor(space, poly.locs[endloc_head])
    possible_moves_tail = neighbor(space, poly.locs[endloc_tail])

    # this automatically chooses the filter for that geometry
    avail_f(point) = avail_filter(space, point) 
    
    moves_head = filter(avail_f, possible_moves_head)
    moves_tail = filter(avail_f, possible_moves_tail)
    
    append!([(endloc_head, endloc_tail, point) for point in moves_head] , # moves head, kills tail
    [(endloc_tail, endloc_head, point) for point in moves_tail] ) #moves tail, kills head
    
end

function update_reverse_checking_space(move::SnakeMove, space::Space, oldpoint, newpoint, poly)
    # handle the reverse checking space
    xnew, ynew = newpoint
    xold, yold = oldpoint
    assert(space.rspace[xnew, ynew, 1] == 0 && space.rspace[xold, yold, 1] != 0)
    
    polyid = space.rspace[xold, yold, 1]
    space.rspace[xold, yold, 1], space.rspace[xold, yold, 2] = 0, 0
    
    for (idx, (xtemp, ytemp)) in enumerate(poly.locs)
        space.rspace[xtemp, ytemp, 1] = polyid
        space.rspace[xtemp, ytemp, 2] = idx
    end
end

function move(move::SnakeMove, space::Space, polyid::Int, polytyp::ASCIIString)
    if polytyp == "sim"
        poly = space.Sims[polyid]
        sitevalue = 1 # code for free sim
    elseif polytyp == "sumo"
        poly = space.Sumos[polyid]
        sitevalue = -1 # code for free sumo
    else
        error("polymer type undefined")
    end
    
    nbr_bond_inc = 0
    
    possible_moves = get_possible_moves(move, space, poly)
    if isempty(possible_moves)
        return (false, nbr_bond_inc)
    end    
    (pointid, killpointid, newpoint) = possible_moves[rand(1:end)]
    oldpoint = poly.locs[killpointid]       
    
    
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
        
        # make changes on polymer.locs
        if pointid == 1 # move towards head
            pop!(poly.locs)
            unshift!(poly.locs, newpoint)
        elseif killpointid == 1 # move towards tail
            shift!(poly.locs)
            push!(poly.locs, newpoint)
        else
            error("end type undefined")
        end

        # handle the reverse checking space (This part needs to use specific geometry)
        update_reverse_checking_space(move, space, oldpoint, newpoint, poly)

        # try building new bonds
        if !isempty(bond_choice)
            bpoint = bond_choice[rand(1:end)]
            create_bond(space, newpoint, bpoint)
        end
        return (true, nbr_bond_inc)
    else
        return (false, nbr_bond_inc)
    end

end