module SpaceModule

export Polymer, Space
export initialize, neighbor, safe_remove, safe_create, create_bond, can_build_bond, exist_bond, in_a_bond, can_build_bond_tentative


using PyCall
@pyimport itertools

NOBOND = 0

type Polymer
    poly_id :: Int
    locs :: Array{Tuple{Int, Int}, 1}
    Polymer() = new(0, [])
end

type Space
    NSim :: Int
    LSim :: Int
    NSumo :: Int
    LSumo :: Int
    Lx :: Int
    Ly :: Int
    space :: Array{Int, 2}
    bond :: Array{Int, 3}
    rspace :: Array{Int, 3}
    Sims :: Array{Polymer, 1}
    Sumos :: Array{Polymer, 1}
    SimId :: Int
    SumoId :: Int
    
    Space(NSim :: Int,
            LSim :: Int,
            NSumo :: Int,
            LSumo :: Int,
            Lx :: Int,
            Ly :: Int) = 
    new(NSim, LSim, NSumo, LSumo, Lx, Ly,
    zeros(Lx, Ly),
    zeros(Lx, Ly, 2),
    zeros(Lx, Ly, 2),
    Polymer[],
    Polymer[],
    0,
    0)
end   


function place!(space::Space, typ::ASCIIString, poly_id::Int, locs::Array{Tuple{Int, Int}, 1})
    if typ == "sim"
        spacetype = 1
        polymers = space.Sims
    elseif typ == "sumo"
        spacetype = -1
        polymers = space.Sumos
    else
        error("Unrecognized polymer type.")
    end
    poly = Polymer()
    poly.poly_id = poly_id
    poly.locs = locs
    push!(polymers, poly)
    for pair in enumerate(locs)
        idx, (x, y) = pair
        space.space[x, y] = spacetype
        space.rspace[x, y, 1] = poly_id
        space.rspace[x, y, 2] = idx
    end
end

function Cprod(x::Array{Int, 1}, y::Array{Int, 1})
    convert(Array{Tuple{Int, Int}, 1}, collect(itertools.product(x, y)))
end

function dilute_init(space::Space, direction::Char)
    if direction == 'v'
        dl = div(space.Lx, (space.NSim + space.NSumo))
        for x in collect(1:space.NSim)
            place!(space, "sim", x, Cprod([(x-1)*dl+1], collect(1:space.LSim)))
        end
        for x in collect(1:space.NSumo)
            place!(space, "sumo", x, Cprod([space.Lx-(x-1)*dl], collect(1:space.LSumo)))
        end
    elseif direction == 'h'
        dl = div(space.Ly, (space.NSim + space.NSumo))
        for y in collect(1:space.NSim)
            palce!(space, "sim", y, Cprod([1:space.LSim], collect((y-1) * dl)+1))
        end
        for y in collect(1:space.NSumo)
            place!(space, "sumo", y, Cprod(collect(1:space.LSumo), [(space.Ly-(y-1))*dl]))
        end
    else
        error("initialization direction undefined")
    end
end


function initialize(space::Space)
    if max(space.LSim, space.LSumo) <= space.Ly && space.NSim+space.NSumo < space.Lx
        dilute_init(space, 'v')
    elseif max(space.LSim, space.LSumo) <= space.Lx && space.NSim+space.NSumo < space.Ly
        dilute_init(space, 'h')
    else
        error("Polymer too dense, cannot initiate.")
    end
end
    
function neighbor(space, point):
    """
    Get a set of neighbors of a point
    """
    x, y = point
    four_points = [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]
    [(mod(pp[1]-1, space.Lx)+1, mod(pp[2]-1, space.Ly)+1) for pp in four_points] 
end

function safe_remove(space, point):
    """
    Removes the monomer at point and
    safely register the change except on polymer.locs and rspace
    """
    x, y = point
    assert(space.space[x, y] != 0)
    space.space[x, y] = 0

    xb, yb = space.bond[x, y, 1], space.bond[x, y, 2]
    if xb != NOBOND && yb != NOBOND # we do have a bond to remove 
        assert(space.bond[xb, yb, 1] == x && space.bond[xb, yb, 2] == y) # assert bond reversibility
        assert(abs(space.space[xb, yb]) == 2) # assert the bond consistency
        space.bond[x, y, 1], space.bond[x, y, 2] = NOBOND, NOBOND # remove the bond from (x,y)
        space.bond[xb, yb, 1], space.bond[xb, yb, 2] = NOBOND, NOBOND # remove the bond to (x,y)
        space.space[xb, yb] = sign(space.space[xb, yb]) # register the space site to be no-bond
    end
end



function safe_create(space, point, polyvalue)
    """
    Create the monomer at point and
    safely register the change except on polymer.locs and rspace
    """
    x, y = point
    assert(space.space[x, y] == 0)
    space.space[x, y] = polyvalue 
end

function create_bond(space, point1, point2)
    """
    Create the bond between point1 = (x1, y1) and point2 = (x2, y2)
    Note that bond information is registered in self.bond
    """
    x1, y1 = point1
    x2, y2 = point2
    assert(can_build_bond(space, point1, point2)) # assert the bond consistency
    assert((space.bond[x1, y1, 1] == NOBOND) && (space.bond[x1, y1, 2] == NOBOND)) # assert the bond nonexistance
    assert((space.bond[x2, y2, 1] == NOBOND) && (space.bond[x2, y2, 2] == NOBOND)) # assert the bond nonexistance
    space.bond[x1, y1, 1], space.bond[x1, y1, 2] = x2, y2
    space.bond[x2, y2, 1], space.bond[x2, y2, 2] = x1, y1
    space.space[x1, y1] *= 2
    space.space[x2, y2] *= 2
end

function can_build_bond(space, point1, point2)
    x1, y1 = point1
    x2, y2 = point2
    space.space[x1, y1] * space.space[x2, y2] == -1
end

function exist_bond(space, point1, point2)
    x1, y1 = point1
    x2, y2 = point2
    space.space[x1, y1] * space.space[x2, y2] == -4
end

function in_a_bond(space, point)
    x, y = point
    abs(space.space[x, y]) != 1
end

function can_build_bond_tentative(space, sitevalue, bpoint)
    xb, yb = bpoint
    space.space[xb, yb] * sitevalue == -1
end
end