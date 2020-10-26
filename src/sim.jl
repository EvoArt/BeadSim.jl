using Pkg
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("Combinatorics")
Pkg.add("DataFrames")
using Distributions
using Random
using Combinatorics
using DataFrames

function Amat(s,c,alpha)
    """
    Genereate adjacency matrix of size s,
    with connectance c
    and interaction strength alpha
    """
A = zeros(s, s)

for i in 1:s
    for j in 1:s
    if rand(Uniform(0,1))< c/2
        A[i,j] = rand(Normal(0,alpha))
    end
end
end

return A
end

function Board(s)
    """
    Create s*s matrix for simulation
    to take place on
    """
    return zeros(Int, s, s)
end

function intc(x)
    """
    get coordinates for interactions
    """
    intCoords = []
    for i in -1:1
        for j in -1:1
            if i * j != i+j
                push!(intCoords, (i,j))
            end
        end
    end
    return intCoords
end


function fitnessModifier(board, A, intCoords,focal)
    """
    Calculate the change in fitness of focal individual
    due to local interactions
    """
    modifier = 0


    partners = [Int(board[p[1],p[2]])  for p in intCoords]

    for partner in partners
        if partner != 0
        modifier += A[focal,partner] #
    end
    end
    return modifier
end

function setBoard(board, density, nSpecies)
    """
    set out initial species distributions
    """

    for i in size(board)[1]
        for j in size(board)[1]
        if rand(Uniform(0,1))< density
            board[i,j] = rand(1:nSpecies)
        end
    end
    end
    return board
end
function checkSpace(board, intCoords)

    if sum([board[i[1],i[2]] == 0 for i in intCoords])!= 0
        return true
    else
        return false
    end
end

function Birth(board, intCoords,focal)

coords = [(i[1],i[2]) for i in intCoords if board[i[1],i[2]] == 0]
coord = coords[rand(1:length(coords))]
board[coord[1],coord[2]] = focal
return board

end

function cull(board, maxDens)
    """
    density regulation
    """
    Size = size(board)[1]

    density = sum([i != 0 for i in board])# current density
    while density > maxDens # if/while density is too high
        board[rand(1:Size),rand(1:Size)] = 0 #kill at random
        density = sum([i != 0 for i in board])# current density
    end
    return board
end

function disperse(board, sections)
    """
    Shuffle the peices on the board
    """
    Size = size(board)[1]
    deme = Size ÷ sections # length of individual demes
    demes = [i for i in 1:deme:Size] # list of starting points for demes
    for i in demes
        for j in demes
            # shuffle the individual demes
            #each deme is the value from the list + the deme length
            board[i: i+deme-1,j:j+deme-1] = shuffle(board[i: i+deme-1,j:j+deme-1])
        end
    end
    return board
end

function sim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens )
    # set up board and adjacency matrix
    board = Board(boardSize)
    board = setBoard(board, density, nSpecies)
    Amatrix = Amat(nSpecies, c, alpha)
    intCoords = intc(1)
    moves = combinations(hcat(1:boardSize,1:boardSize),2) |> collect

    # initialise counters
    births = 0
    richness = length(unique(board))
    gens = 0

    #run the simulation
    while births < bLim
        gens+=1
        #get randomly ordered vector of board indeces
        for (i,j) in shuffle(unique(moves))
                if board[i,j] != 0 #not empty
                    #evaluate fitness
                    focal = board[i,j] # species id

                    locality = [(i,j) .+ ic for ic in intCoords]

                    locality = locality[[l for l in 1:8 if checkbounds(Bool,board,locality[l][1],locality[l][2])]]

                    if rand(Uniform(-30.5,-30))< fitnessModifier(board, Amatrix, locality,focal)
                        # then give birth if space available
                        if checkSpace(board, locality) == true
                        board = Birth(board,locality,focal)
                        births += 1
                        end
                    end
                end

        end
        #dispersal
        if rand(Uniform(0,1)) < dispersal
            board = disperse(board, demes)
        end
        #density regulation
        board = cull(board,maxDens)

        # check for extinctions
        if length(unique(board)) < richness # if an extinction has occurred
            richness = length(unique(board)) #update richness
            births = 0 # reset counter
        end
    end
        # return data!
        return [density c nSpecies boardSize dispersal demes^2 alpha maxDens richness gens]
end

function threadSim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens)
# initialise dictionary for storing data
D = Dict()
print(Threads.nthreads())
Threads.@threads for i in 1:Threads.nthreads()
    D[Threads.threadid()] = sim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens)
end
return [vcat(V) for (K, V) in D ]
end


#x = sim(0.5, 100, 0.5, 100, 1, 4, 1,10,800 )
#x
y = threadSim(0.5, 100, 0.5, 100, 1, 4, 1,10,800 )
y

# run the project
dat = Array{Float64}(undef, 0, 10)
for c in LinRange(0.1,1,10)
    for nSpecies in [50,100,200]
        for density in [0.2,0.5,0.8]
            for boardSize in [100,1000]
                for dispersal in LinRange(0.1,1,10)
                    for demes in [boardSize ÷ d for d in 1:boardSize if boardSize/d == boardSize ÷ d]
                        for alpha in LinRange(0.1,2,20)
                            for bLim in boardSize^2 .* [0.1,0.5,1,2,3,]

                                dat = vcat(dat,threadSim(c, nSpecies, density, boardSize, dispersal, demes, alpha, bLim ,density))
                            end
                        end
                    end
                end
                df = DataFrame(dat)
                names!(df, [:c, :nSpecies, :density, :boardSize, :dispersal, :demes, :alpha,:bLim ,:maxDens, :gens])
                CSV.write("space.csv", df)

            end
        end
    end
end
