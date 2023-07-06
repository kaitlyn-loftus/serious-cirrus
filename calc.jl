using Revise 
push!(LOAD_PATH,"./src/")
using ΜΦ
using LHS
using JLD2 

"""
perform microphysics calculations for Ding, Loftus, & Wordsworth (submitted)
[doi:TBD]

these calculations used for Figure 6 and Supplemental Figures 9-15


please contact Kaitlyn Loftus (kaitlyn.loftus@columbia.edu) with comments / questions / concerns
"""

println("doing microphysics calculations for Ding, Loftus, & Wordsworth (submitted):")
println("'Idealized 3D GCM and microphysical modeling of water ice clouds on early Mars:\nImplications for climate'")
println("\n\n")

# set up conditions to span 

# temperature 
nT = 100
Ts = LinRange(200.,270.,nT)
# pressure 
ps = [1e3,1e4,1e5]
np = length(ps)

# supersaturation 
𝒮s = [1e-4,1e-3,1e-2,1e-1]
n𝒮 = length(𝒮s)

# set up array to store results 
tstars = zeros(2,nT,np,n𝒮)


# do analytic spherical calculations 
# for Fig 6
println("starting calculations for spherical particles")
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=0.)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:ana_rt)[1]
            # do Earth calc 
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=0.)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:ana_rt)[1]
        end
    end
end
println("finished calculations for spherical particles")
# save 
jldsave("out/sphere.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)

# do analytic spherical calculations for Mars with Earth gravity  
# for SFig 9
println("starting calculations for spherical particles for Mars with Earth gravity")
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice(p,T,𝒮,:Earth,:CO2,isCc=true,rat=1.,issphereCS=true,r₀=0.)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:ana_rt)[1]
            # do Earth calc 
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=0.)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:ana_rt)[1]
        end
    end
end
println("finished calculations for spherical particles for Mars with Earth gravity")
# save 
jldsave("out/MarsEarthg.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)

# do analytic plate calculations 
# for SFig 10
println("starting calculations for plate-like particles")
rat_plate = 0.1
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=rat_plate,issphereCS=false,r₀=0.)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:ana_rt)[1]
            # do Earth calc
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=rat_plate,issphereCS=false,r₀=0.)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:ana_rt)[1]
        end
    end
end
println("finished calculations for plate-like particles")
# save 
jldsave("out/plate.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)


# do analytic column calculations 
# for SFig 11
println("starting calculations for column-like particles")
rat_col = 10.
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=rat_col,issphereCS=true,r₀=0.)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:ana_rt)[1]
            # do Earth calc 
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=rat_col,issphereCS=true,r₀=0.)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:ana_rt)[1]
        end
    end
end
println("finished calculations for column-like particles")
# save 
jldsave("out/column.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)

# do numerical spherical calculations 
# for SFig 13
println("starting calculations for spherical particles with numerical integration, αD=1")
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:num_rt)[1]
            # do Earth calc
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:num_rt)[1]
        end
    end
end
println("finished calculations for spherical particles with numerical integration, αD=1")
# save 
jldsave("out/noncontiuum.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)

# do numerical spherical calculations with αD=0.1
# for SFig 14
println("starting calculation for spherical particles with numerical integration, αD=0.1")
for (k,𝒮) ∈ enumerate(𝒮s)
    for (j,p) ∈ enumerate(ps)
        for (i,T) ∈ enumerate(Ts)
            # do Mars calc
            pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8,αD=0.1)
            tstars[1,i,j,k] = solvetHfall(pμϕ_Mars,:num_rt)[1]
            # do Earth calc
            pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8,αD=0.1)
            tstars[2,i,j,k] = solvetHfall(pμϕ_Earth,:num_rt)[1]
        end
    end
end
println("finished calculations for spherical particles with numerical integration, αD=0.1")
# save 
jldsave("out/noncontiuum_αD0p1.jld2",Ts=Ts,ps=ps,𝒮s=𝒮s,tstars=tstars,np=np,nT=nT,n𝒮=n𝒮)

# do LHS calculations

# number of p-T-𝒮 conditions 
nLHS = 100

logpTlog𝒮_LHS = doLHSND([3,Ts[1],-4],[5,Ts[end],-1],nLHS)
pT𝒮_LHS = deepcopy(logpTlog𝒮_LHS)
pT𝒮_LHS[:,1] .= 10. .^ logpTlog𝒮_LHS[:,1]
pT𝒮_LHS[:,3] .= 10. .^ logpTlog𝒮_LHS[:,3]

# axis ratios to span for each LHS p-T-𝒮 condition 
nrats = 100 
rats = 10 .^ LinRange(-2.,2.,nrats)
# mass accommodation coefficients to span for each LHS p-T-𝒮 condition 
nα = 100
αDs = 10 .^ LinRange(-3.,0.,nα)


# set up arrays to store LHS results 
tstar_ratsrats = zeros(nLHS,nrats)
tstar_ratsαDs = zeros(nLHS,nα)

# do LHS axis ratio calculations 
# for SFig 12
println("starting calculations for axis ratio with LHS")
for (j,rat) ∈ enumerate(rats)
    for i ∈ 1:nLHS
        p,T,𝒮 = pT𝒮_LHS[i,:]
        # do Mars calc
        pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=rat,issphereCS=true,r₀=0.)
        # do Earth calc
        pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=rat,issphereCS=true,r₀=0.)
        tstar_ratsrats[i,j] = solvetHfall(pμϕ_Mars,:ana_rt)[1] / solvetHfall(pμϕ_Earth,:ana_rt)[1]
    end
end
println("finished calculations for axis ratio with LHS")
# save 
jldsave("out/ratio_LHS.jld2",Ts=pT𝒮_LHS[:,2],ps=pT𝒮_LHS[:,1],𝒮s=pT𝒮_LHS[:,3],rats=rats,rat_tstar=tstar_ratsrats,nLHS=nLHS)

# do LHS αD calculations 
# for SFig 15
println("starting calculations for αD with LHS")
for (j,αD) ∈ enumerate(αDs)
    for i ∈ 1:nLHS
        p,T,𝒮 = pT𝒮_LHS[i,:]
        # do Mars calc
        pμϕ_Mars = setup_pμϕ_ice_Mars(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8,αD=αD)
        # do Earth calc
        pμϕ_Earth = setup_pμϕ_ice_Earth(p,T,𝒮,isCc=true,rat=1.,issphereCS=true,r₀=1e-8,αD=αD)
        tstar_ratsαDs[i,j] = solvetHfall(pμϕ_Mars,:num_rt)[1] / solvetHfall(pμϕ_Earth,:num_rt)[1]
    end
end
println("finished calculations for αD with LHS")
# save 
jldsave("out/αD_LHS.jld2",Ts=pT𝒮_LHS[:,2],ps=pT𝒮_LHS[:,1],𝒮s=pT𝒮_LHS[:,3],αDs=αDs,αD_tstar=tstar_ratsαDs,nLHS=nLHS)



