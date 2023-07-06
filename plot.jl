using Revise 
push!(LOAD_PATH,"./src/")
using ΜΦ
using JLD2
using CairoMakie 

"""
plot microphysics calculations for Ding, Loftus, & Wordsworth (submitted)
[doi:TBD]

this script makes Figure 6 and Supplemental Figures 9-15
using the results of calc.jl script 


please contact Kaitlyn Loftus (kaitlyn.loftus@columbia.edu) with comments / questions / concerns
"""

function plottstar(name_exp,fignum;cmap=:wong)
    """
    plot tstar for Mars conditions in days and in units of tstar for equivalent Earth conditions 
    """
    # load calcs 
    res = load("out/"*name_exp*".jld2")
    tstars = res["tstars"]
    Ts = res["Ts"]
    ps = res["ps"]
    𝒮s = res["𝒮s"]
    np = res["np"]
    nT = res["nT"]
    n𝒮 = res["n𝒮"]


    # plot results  
    pstyles = [:dot,:solid,:dash]
    lws = [3,5,3]
    if cmap==:wong
        cs = Makie.wong_colors()
    else
        cs = cgrad(cmap, 5, categorical = true)
    end

    𝒮strings = ["⁴","³","²","¹"]
    pstrs = ["³","⁴","⁵"]


    fig = Figure(fontsize=20)
    axt = Axis(fig[1,1],xlabel="𝘛 [K]", ylabel=rich(rich("t",font=:italic),superscript("∗"), " [days]"),yscale=log10,yminorticks=IntervalsBetween(9),yminorticksvisible=true,title="(a)",titlealign=:left)
    ylims!(axt,1e-1,1e2)
    xlims!(axt,Ts[1],Ts[end])
    hspan!(axt,15,1e2,color=(cs[n𝒮+1],0.35))
    hidexdecorations!(axt,grid=false,ticks=false)
    axrat = Axis(fig[2,1],xlabel="𝘛 [K]", ylabel=rich(rich("t",font=:italic),superscript("∗"),subscript("Mars")," [",rich("t",font=:italic),superscript("∗"),subscript("Earth"),"]"),title="(b)",titlealign=:left)
    if name_exp=="MarsEarthg"
        ylims!(axrat,0.9,1)
        axrat.yticks = [0.9,0.95,1.]
    else
        ylims!(axrat,2,3)
    end
    xlims!(axrat,Ts[1],Ts[end])
    for (k,𝒮) ∈ enumerate(𝒮s)
        for (j,p) ∈ enumerate(ps)
            lines!(axt,Ts,converts2Earthday.(tstars[1,:,j,k]),color=(cs[k],1),linestyle=pstyles[j],linewidth=lws[j],label="𝒮 = 10⁻"*𝒮strings[k]*", 𝘱 = "*pstrs[j]*" Pa") 
            lines!(axrat,Ts,tstars[1,:,j,k] ./ tstars[2,:,j,k],color=(cs[k],1),linestyle=pstyles[j],linewidth=lws[j]) 
        end
    end
    text!(axt,245.,20,text=rich("warm climate regime\ncloud lifetime timescales\n(",font=:bold,rich("τ",font=:bold_italic),subscript("c")," ≥ ",rich("τ",font=:bold_italic),subscript("tran"),")"),color=cs[n𝒮+1])
    rowsize!(fig.layout,2,Relative(0.25))

    axs = [axt,axrat]
    yspace = maximum(tight_yticklabel_spacing!, axs)
    for ax ∈ axs
        ax.yticklabelspace = yspace
    end
    # set up legend 
    fig[1:2, 2] = Legend(fig, [[LineElement(color=cs[i],linewidth=4) for i ∈ 1:n𝒮], 
    [LineElement(color=:black,linewidth=lws[i],linestyle=pstyles[i]) for i ∈ 1:np]],
    [["10⁻⁴","10⁻³","10⁻²","10⁻¹"], ["10³","10⁴","10⁵"]],
    ["𝓢 [ ]", "𝙥 [Pa]"], patchsize = (50, 35))

    save("figs/fig"*fignum*"_"*name_exp*".pdf",fig)
    nothing 
end

println("plotting microphysics calculations for Ding, Loftus, & Wordsworth (submitted):")
println("'Idealized 3D GCM and microphysical modeling of water ice clouds on early Mars:\nImplications for climate'")
println("\n\n")

# plot tstar as a function of T under different assumptions 
# Fig 6, Sfigs 9-11, 13-14
name_exps = ["sphere","MarsEarthg","plate","column","noncontiuum","noncontiuum_αD0p1"]
fignums = ["6","S9","S10","S11","S13","S14"]
for i ∈ 1:length(name_exps)
    plottstar(name_exps[i],fignums[i])
end


# plot tstar as a function of particle axis ratio and mass accommodation coefficient 
# with LHS environmental conditions 
LHS_rat = load("out/ratio_LHS.jld2")
rats = LHS_rat["rats"]
nLHS = LHS_rat["nLHS"]
rat_tstar = LHS_rat["rat_tstar"]

# Sfig 12 
α = 0.25
fig = Figure(fontsize=20)
ax = Axis(fig[1,1],xscale=log10,xminorticks=IntervalsBetween(9),xminorticksvisible=true)
for i ∈ 1:nLHS
    lines!(ax,rats,rat_tstar[i,:],color=(:black,α),linewidth=3)
end
xlims!(ax,1e-2,100)
ylims!(ax,2,3)
ax.xlabel = "axis ratio [ ]"
ax.ylabel = rich(rich("t",font=:italic),superscript("∗"),subscript("Mars")," [",rich("t",font=:italic),superscript("∗"),subscript("Earth"),"]")
save("figs/figS12_tstarratrat.pdf",fig)


# Sfig 15
LHS_rat = load("out/αD_LHS.jld2")
αDs = LHS_rat["αDs"]
nLHS = LHS_rat["nLHS"]
αD_tstar = LHS_rat["αD_tstar"]

α = 0.25
fig = Figure(fontsize=20)
ax = Axis(fig[1,1],xscale=log10,xminorticks=IntervalsBetween(9),xminorticksvisible=true)
for i ∈ 1:nLHS
    lines!(ax,αDs,αD_tstar[i,:],color=(:black,α),linewidth=3)
end
xlims!(ax,1e-3,1)
ylims!(ax,1.5,3)
ax.xlabel = "mass accommodation coefficient [ ]"
ax.ylabel = rich(rich("t",font=:italic),superscript("∗"),subscript("Mars")," [",rich("t",font=:italic),superscript("∗"),subscript("Earth"),"]")
save("figs/figS15_tstarratαD.pdf",fig)

println("completed making plots")
println("plots saved in figs/ directory")