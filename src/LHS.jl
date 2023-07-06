module LHS

using Random 
using CairoMakie

export doLHSND

"""
functions for Latin hypercube sampling
"""

function dounscale01(x′,xₘᵢₙ,xₘₐₓ)
    """
    convert x′ ∈ [0,1] to x ∈ [xₘᵢₙ,xₘₐₓ]
    """
    Δx = xₘₐₓ - xₘᵢₙ
    x′ * Δx + xₘᵢₙ
end


function doLHSND(xₘᵢₙs,xₘₐₓs,n;seed=42)
    """
    do Latin hypercube sampling in N dimensions 
    inputs:
        * xₘᵢₙs [Array of dims: N]- minimum value for each dimension to sample 
        * xₘₐₓs [Array of dims: N] - maximum value for each dimension to sample 
        * n [Int] - number of Latin hypercube samples 
        * seed [Int] - seed for random number generator 
    output:
        * LHS [Array of dims: n × N] - Latin hypercube sample 
    """
    # set up random number generator 
    rng = MersenneTwister(seed)
    # number of dimensions to sample over 
    N = length(xₘᵢₙs)
    # throw error if incorrect inputs 
    if length(xₘₐₓs)!=N
        DomainError("xₘᵢₙs and xₘₐₓs do not have the same size! they must have the same size.")
    end
    # initial random sampling 
    # n samples over N dimensions 
    # here only between 0 and 1 
    LHS = rand(rng,n,N)
    # iterate over dimensions to sample 
    for iD ∈ 1:N
        # set up range of values for each dimension 
        xranges = LinRange(xₘᵢₙs[iD],xₘₐₓs[iD],n+1)
        # randomly assign position in xrange to LHS position 
        # (so samples aren't from the same relative portion of xrange across dimensions)
        for (iLHS,irange) ∈ enumerate(shuffle(rng, Vector(1:n)))
            # change scale from 0 to 1 to actual requested range 
            LHS[iLHS,iD] = dounscale01(LHS[iLHS,iD],xranges[irange],xranges[irange+1])
        end
    end
    LHS
end

function checkLHS()
    """
    visually check LHS works as it should
    """
    n = 5
    xₘᵢₙs = [3.,9.]
    xₘₐₓs = [5.,15.]
    LHS = doLHSND(xₘᵢₙs,xₘₐₓs,n)
    fig = Figure()
    ax = Axis(fig[1,1])
    for i ∈ 1:2
        xranges = LinRange(xₘᵢₙs[i],xₘₐₓs[i],n+1)
        if i==1
            vlines!(ax,xranges,color=:black,linestyle=:dash)
            
        else
            hlines!(ax,xranges,color=:black,linestyle=:dash)
        end
    end
    scatter!(ax,LHS[:,1],LHS[:,2],color=(:black,0.9))
    save("figs/checkLHS.pdf",fig)
end

end