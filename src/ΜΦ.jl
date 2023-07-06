module ΜΦ

using CairoMakie 
using Roots 
using OrdinaryDiffEq

export setup_pμϕ_ice,setup_pμϕ_ice_Mars,setup_pμϕ_ice_Earth,solvetHfall,converts2Marsday,converts2Earthday

"""
functions and structs for microphysics calculations 
"""

# constants 
# Boltzmann constant
const kB = 1.380649e-23 # [J K⁻¹] 
# molar mass of CO2
const μCO2 = 44.01e-3 # [kg mol⁻¹]
# molar mass of H2O
const μH2O = 18.01528e-3 # [kg mol⁻¹]
# molar mass of (dry) Earth air
const μEarth = 28.9647e-3 # [kg mol⁻¹]
# ideal gas constant 
const R = 8.3145 # [J K⁻¹ mol⁻¹]
# specific gas constant for water 
const Rₕ₂ₒ = R / μH2O  # [J kg⁻¹ K⁻¹]
# seconds in a Mars day 
const sinMarsday = 88620. # [s day⁻¹]
# seconds in an Earth day 
const sinEarthday = 86400. # [s day⁻¹]
# gravitational acceleration at Mars surface 
const gMars = 3.72 # [m s⁻¹]
# gravitational acceleration at Earth surface 
const gEarth = 9.81 # [m s⁻¹]


struct Paramμϕ
    """
    immutable structure for storing properties of
    cloud particle & local atmosphere
    needed for microphysics (μϕ) calculations
    """ 
    pair::Float64 # [Pa] local air pressure 
    Tair::Float64 # [K] local air temperature 
    𝒮air::Float64 # [ ] supersaturation (i.e., RH - 1 = pₕ₂ₒ(T)/pₕ₂ₒₛₐₜ(T) - 1)
    ρair::Float64 # [kg m⁻³] local air density 
    ηair::Float64 # [Pa s⁻¹] local air dynamic viscosity
    Rair::Float64 # [J K⁻¹ kg⁻¹] air specific gas constant (i.e., R/μair)
    cₚair::Float64 # [J kg⁻¹ K⁻¹] local air specific heat at constant pressure
    D::Float64 # [m² s⁻¹] diffusivity of H2O in air 
    K::Float64 # [W m⁻¹ K⁻¹] air thermal conductivity
    αD::Float64 # [ ] mass accommodation coefficient
    αT::Float64 # [ ] thermal accomodation coefficient
    mfp::Float64 # [m] local air mean free path 
    ρc::Float64 # [kg m⁻³] cloud particle density
    Fₖ::Float64 # [s m⁻¹ kg⁻¹] term for latent heat effects for condensational growth 
    Fd::Float64 # [s m⁻¹ kg⁻¹] term for diffusional supply of water vapor for condensational growth 
    rat::Float64 # [ ] cloud particle axis ratio (i.e., particle semi-axis along axis of rotation / other semi-axes)
    fshape::Float64 # [ ] shape correction factor for drag coefficient CD at low Re
    Cshape::Float64 # [ ] shape correction factor for drag coefficient CD at high Re
    issphereCS::Bool # [bool] boolean for if particle cross section is sphere-like (determines which CDstar calc to use)
    Cdivreq::Float64 # [ ] capacitance for particle shape divided by particle equivalent radius 
    Ashape::Float64 # [ ] shape correction factor for particle cross sectional area 
    r₀::Float64 # [m] initial cloud particle equivalent radius 
    g::Float64 # [m s⁻²] local graviational acceleration
    H::Float64 # [m] local atmospheric scale height 
    isCc::Bool # [bool] boolean for including Cunningham slip correction factor 
end

# make Paramμϕ structures broadcast like a single thing
Base.broadcastable(x::Paramμϕ) = Ref(x)

function calcmfp(p,T,d)
    """
    calculate air mean free path 
    inputs:
        * p [Pa] - air pressure 
        * T [K] - air temperature 
        * d [m] - (average) diameter of air molecule
    output:
        * mfp [m] - mean free path  
    """
    kB/√(2)/π*T/p/d/d
end

function calcmfpCO2(p,T)
    """
    calculate mean free path for CO2 gas 
    d for CO2 from Jacobson (2005) Table 16.1
    [doi:10.1017/cbo9781139165389] 
    inputs:
        * p [Pa] - air pressure 
        * T [K] - air temperature 
    output:
        * mfp [m] - mean free path  
    """
    d = 0.453e-9 # [m]
    calcmfp(p,T,d)
end

function calcmfpEarth(p,T)
    """
    calculate mean free path for Earth air gas 
    d for Earth air from Jacobson (2005) Table 16.1
    [doi:10.1017/cbo9781139165389] 
    inputs:
        * p [Pa] - air pressure 
        * T [K] - air temperature 
    output:
        * mfp [m] - mean free path  
    """
    d = 0.367e-9 #[m] 
    calcmfp(p,T,d)
end

function calcCc(mfp,r)
    """
    calculate the Cunningham slip correction factor 
    Knudsen number Kn ≡ air mean free path / particle radius
    correction factor from Allen & Raabe (1982)
    [doi:10.1016/0021-8502(82)90019-2]
    inputs:
        * mfp [m] - air mean free path 
        * r [m] - particle radius 
    outputs:
        * Cc [ ] - Cunningham slip correction factor 
    """
    α = 1.115
    β = 0.471
    γ = 0.596
    Kn = mfp/r
    1. + Kn*(α + β*exp(-γ/Kn))
end


function calcD′(r::Float64,pμϕ::Paramμϕ)::Float64
    """
    calculate the modified water diffusivity
    following Pruppacher & Klett (2010) eq (13.14)
    [doi:10.1007/978-0-306-48100-0]
    assuming Δᵥ = Cunningham factor × mfp
    inputs:
        * r [m] - particle equivalent radius 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    output:
        * D′ [m² s⁻¹] - modified diffusivity
    """
    pμϕ.D/(r / (r + calcCc(pμϕ.mfp,r)*pμϕ.mfp) + pμϕ.D/pμϕ.αD/r*(2*π/Rₕ₂ₒ/pμϕ.Tair)^0.5)
end

function calcK′(r::Float64,pμϕ::Paramμϕ)::Float64
    """
    calculate the modified thermal conductivity
    following Pruppacher & Klett (2010) eq (13.20)
    [doi:10.1007/978-0-306-48100-0]
    assuming Δₜ = Cunningham factor × mfp
    inputs:
        * r [m] - particle equivalent radius 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    output:
        * K′ [W m⁻¹ K⁻¹] - modified thermal conductivity

    """
    pμϕ.K/(r / (r + calcCc(pμϕ.mfp,r)*pμϕ.mfp) + pμϕ.K/pμϕ.αT/r/pμϕ.cₚair/pμϕ.ρair*(2*π/pμϕ.Rair/pμϕ.Tair)^0.5)
end

function calcpsatice(T)
    """
    calculate ice saturation pressure 
    valid for T ∈ [50,273.16] K
    from Wagner et al. (2011) eq (4) and table 3
    [doi:10.1063/1.3657937]
    input:
        * T [K] - local tempature 
    output:
        * p [Pa] - saturation pressure of ice
    """
    θ = T/273.16
    a = [-0.212144006e2,0.273203819e2,-0.610598130e1]
    b = [0.333333333e-2,0.120666667e1,0.170333333e1]
    lnπ = 0
    for i ∈ 1:3
        lnπ += a[i]*θ^(b[i])
    end
    lnπ /= θ
    exp(lnπ)*611.657 # [Pa]
end

function calcρice(T)
    """
    calculate density of ice 
    from Pruppacher & Klett (2010) eq 3.2 converted to SI units
    [doi:10.1007/978-0-306-48100-0]
    (factor of 1e3 converts g m⁻³ to kg m⁻³)
    valid for T ∈ [93.15,273.15] K
    input:
        * T [K] - tempature of ice
    output:
        * ρ [kg m⁻³] - density of ice
    """
    Tc = T - 273.15 # [C]
    (0.9167 - 1.75e-4*Tc - 5e-7*Tc^2)*1e3
end

function calcLice(T)
    """
    calculate latent heat of sublimation 
    from Feistel & Wagner (2007) eq (4.8)
    [doi:10.1016/j.gca.2006.08.034]
    input:
        * T [K] - local temperature 
    output:
        * L [J kg⁻¹] - latent heat of sublimation 
    """
    a = [2638742.45418107,
        400983.673912406,
        200812.111806393,
        -1486203.38485336,
        2290451.50230789,
        -1690159.93521118,
        479848.354373932]
    θ = T/273.16
    L = 0
    for i ∈ 0:6
        L += a[i+1]*θ^i
    end
    L
end

function calccpCO2(T)
    """
    calculate specific heat at constant pressure for CO2 gas 
    from Wordsworth & Pierrehumbert (2013) eq (2)
    [doi:10.1088/0004-637X/778/2/154]
    fit for T ∈ [175,600] K
    input:
        * T [K]
    output:
        * cp [J kg⁻¹ K⁻¹]
    """
    574.8 + 0.875*T
end

function calccpEarth(T)
    """
    calculate specific heat at constant pressure for Earth air 
    from Zografos et al. (1987) Table 1 
    [doi:10.1016/0045-7825(87)90003-X]
    fit for T ∈ [100-3000] K
    input:
        * T [K]
    output:
        * cp [J kg⁻¹ K⁻¹]
    """
    (1.3864e-13*T^4 - 6.4747e-10*T^3 + 1.0234e-6*T^2 - 4.3282e-4*T + 1.0613)*1e3
end

function calcKEarth(T)
    """
    calculate K from T for modern Earth composition (dry) air 
    following Beard & Pruppacher (1971)
    as given in Lohmann et al. (2016) eq (7.24)
    [doi:10.1017/CBO9781139087513]
    assumes no pressure dependence and pₕ₂ₒ dependence
    input:
        * T [K] - temperature
    output:
        * K [W m⁻¹ K⁻¹] - thermal conductivity 
    """
    4.1868e-3 * (5.69+0.017*(T-273.15))
end

function calcKCO2(T)
    """
    calculate K from T for CO2 gas 
    fit from CRC data (Section 6, Table 9) 
    [doi:]
    & 
    Hellmann (2014)
    [doi:10.1016/j.cplett.2014.08.057]
    input:
        * T [K] - temperature
    output:
        * K [W m⁻¹ K⁻¹] - thermal conductivity 
    """
    8.372826E-08*T^2 + 2.865834E-05*T + 4.993122E-04
end

function calcDCO2H2O(p,T)
    """
    calculate the diffusivity of H2O vapor in CO2 gas 
    from Hellmann+ (2019)
    [doi:10.1016/j.fluid.2018.11.033]
    inputs:
      * p [Pa] - pressure
      * T [K] - temperature
    output:
      * D [m² s⁻¹] - diffusivity of H2O in CO2
    """
    T_onethird = T^(1/3.)
    S = (-0.09647 + 4.8695 * T^(-1/6.) + 103.7 * T_onethird * exp(-T_onethird) - 4.04e4 * exp(-2. *T_onethird) + 2.1764e6 * exp(-3. * T_onethird)) # Hellmann (2019a) eq 12
    D = 8.31446261815324e-4 * T^1.5 / S / p # Hellmann (2019a) eq 11 modified
    D
end

function calcρCO2(p,T)
    """
    calculate density of CO2 gas 
    inputs: 
        * p [Pa] - pressure 
        * T [K] - temperature 
    outputs:
        * ρ [kg m⁻³] - density 
    """
    p/R/T*μCO2
end

function calcρEarth(p,T)
    """
    calculate density of Earth air 
    inputs: 
        * p [Pa] - pressure 
        * T [K] - temperature 
    outputs:
        * ρ [kg m⁻³] - density 
    """
    p/R/T*μEarth
end

function calcηCO2(T)
    """
    calculate dynamic viscosity of CO2 gas
    from Laesecke & Muzny (2017) eq (4)
    [doi:10.1063/1.4977429]
    fit for T ∈ [100,2000] K
    neglects effects of pressure
    factor of 1e-3 converts mPa s to Pa s
    input:
        * T [K] - temperature 
    output:
        * η [Pa s] - dynamic viscosity 
    """
    a₀ = 1749.354893188350
    a₁ = -369.069300007128
    a₂ = 5423856.34887691
    a₃ = -2.21283852168356
    a₄ = -269503.247933569
    a₅ = 73145.021531826
    a₆ = 5.34368649509278
    Tthird = T^(1/3.)
    Thalf = T^0.5
    Tsixth = T^(1/6.)
    1.0055*Thalf/(a₀ + a₁*Tsixth + a₂*exp(a₃*Tthird) + (a₄+a₅*Tthird)/(exp(Tthird)) + a₆*Thalf)*1e-3
end

function calcηEarth(T)
    """
    calculate dynamic viscosity of Earth air 
    from Zografos et al. (1987) Table 1 
    [doi:10.1016/0045-7825(87)90003-X]
    fit for T ∈ [100-3000] K
        neglects effects of pressure
    input:
        * T [K] - temperature 
    output:
        * η [Pa s] - dynamic viscosity
    """
    2.5914e-15*T^3 - 1.4346e-11*T^2 + 5.0523e-8*T + 4.1130e-6
end


function calcH(T,g,μair)
    """
    calculate atmospheric scale height 
    inputs:
        * T [K] - temperature
        * g [m s⁻²] - gravitational acceleration 
        * μair [kg mol⁻¹] - air molar mass 
    output:
        * H [m] - scale height 
    """
    R*T/μair/g
end

function calcDEarthH2O(p,T)
    """
    calculate water vapor in air diffusivity 
    for modern Earth composition (dry) air 
    following Hall & Pruppacher (1976)
    as given in Lohmann et al. (2016) eq (7.26)
    [doi:10.1017/CBO9781139087513]
    assumes no pₕ₂ₒ dependence
    inputs:
        * p [Pa] - pressure
        * T [K] - temperature
    output:
        * D [m² s⁻¹] - diffusivity of H2O in Earth air
    """
    2.11e-5 * (T/273.15)^1.94 * 101325. / p
end


function calcFₖ(T,L,K)
    """
    calculate term for latent heat effects for condensational growth 
    from Lohmann et al. (2016) eq (7.23)
    [doi:10.1017/CBO9781139087513]
    inputs:
        * T [K] - temperature 
        * L [J kg⁻¹] - latent heat of sublimation 
        * K [W m⁻¹ K⁻¹] - thermal conductivity
    output:
        * Fₖ [s m⁻¹ kg⁻¹] - term for latent heat effects for condensational growth 
    """
    (L/T/Rₕ₂ₒ - 1)*L/T/K
end

function calcFd(T,D,psat)
    """
    calculate term for diffusional supply of water vapor for condensational growth 
    from Lohmann et al. (2016) eq (7.25)
    [doi:10.1017/CBO9781139087513]
    inputs:
        * T [K] - temperature 
        * D [s m⁻²] - diffusivity of water vapor in air 
        * psat [Pa] - saturation pressure
    output:
        * Fd [s m⁻¹ kg⁻¹] - term for diffusional supply of water vapor for condensational growth
    """
    Rₕ₂ₒ*T / D / psat
end


function setup_pμϕ_ice(p,T,𝒮,gtype,airtype;rat=1.,r₀=0.,αD=1,αT=1,issphereCS=true,isCc=false)
    """
    setup Paramμϕ struct for ice calculation 
    inputs:
        * p [Pa] - pressure 
        * T [K] - temperature 
        * 𝒮 [ ] - supersaturation wrt ice 
        * gtype [symbol] - specify planetary gravitational acceleration (:Mars, :Earth)
        * airtype [symbol] - specify air type (:CO2, :Earth)
        * rat [ ] - particle axis ratio
        * r₀ [m] - initial equivalent radius at time 0 
        * αD [ ] - mass accommodation coefficient
        * αT [ ] - thermal accommodation coefficient
        * issphereCS [boolean] - boolean for if particle cross section is sphere-like (determines which CDstar calc to use)
        * isCc [boolean] - whether to correct fall velocity for noncontiuum effects 
    output:
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    """
    # determine properties dependent on air composition 
    if airtype==:CO2
        ρair = calcρCO2(p,T)
        ηair = calcηCO2(T)
        K = calcKCO2(T)
        D = calcDCO2H2O(p,T)
        mfp = calcmfpCO2(p,T)
        μair = μCO2 
        cₚair = calccpCO2(T)
    elseif airtype==:Earth
        ρair = calcρEarth(p,T)
        ηair = calcηEarth(T)
        K = calcKEarth(T)
        D = calcDEarthH2O(p,T)
        mfp = calcmfpEarth(p,T)
        μair = μEarth
        cₚair = calccpEarth(T)
    else
        DomainError(airtype,"airtype value not supported. currently only :CO2 and :Earth supported")
    end
    # set gravitational acceleration 
    if gtype==:Mars
        g = gMars 
    elseif gtype==:Earth 
        g = gEarth 
    else
        DomainError(gtype,"gtype value not supported. currently only :Mars and :Earth supported")
    end
    Rair = R/μair
    ρc = calcρice(T)
    L = calcLice(T)
    psat = calcpsatice(T)
    Fₖ = calcFₖ(T,L,K)
    Fd = calcFd(T,D,psat)
    H = calcH(T,g,μair)
    fshape = calcfshape(rat)
    Cshape = calcCshape(rat)
    Cdivreq = calcCdivreq(rat)
    Ashape = calcAshape(rat)
    Paramμϕ(p,T,𝒮,ρair,ηair,Rair,cₚair,D,K,αD,αT,mfp,ρc,Fₖ,Fd,rat,fshape,Cshape,issphereCS,Cdivreq,Ashape,r₀,g,H,isCc)
    
end

function setup_pμϕ_ice_Mars(p,T,𝒮;rat=1.,r₀=0.,αD=1,αT=1,isCc=false,issphereCS=true)
    """
    setup Paramμϕ struct for ice calculation assuming Mars surface gravity and dry air composition 
    inputs:
        * p [Pa] - pressure
        * T [K] - temperature 
        * 𝒮 [ ] - supersaturation of water vapor wrt ice 
        * rat [ ] - particle axis ratio (rat = 1 => spherical, rat > 1 => prolate spheroid, rat < 1 => oblate spheroid)
        * r₀ [m] - initial particle size at t = 0 s
        * αD [ ] - mass accommodation coefficient 
        * αT [ ] - thermal accommodation coefficient
        * isCc [Bool] - whether to correct particle velocity for non-continuum effects 
        * issphereCS [Bool] - whether particle has approximately spherical cross section (for CD calculation) 
    output:
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info  
    """
    setup_pμϕ_ice(p,T,𝒮,:Mars,:CO2;rat=rat,r₀=r₀,αD=αD,αT=αT,isCc=isCc,issphereCS=issphereCS)
end

function setup_pμϕ_ice_Earth(p,T,𝒮;rat=1.,r₀=0.,αD=1,αT=1,isCc=false,issphereCS=true)
    """
    setup Paramμϕ struct for ice calculation assuming Earth surface gravity and dry air composition 
    inputs:
        * p [Pa] - pressure
        * T [K] - temperature 
        * 𝒮 [ ] - supersaturation of water vapor wrt ice 
        * rat [ ] - particle axis ratio (rat = 1 => spherical, rat > 1 => prolate spheroid, rat < 1 => oblate spheroid)
        * r₀ [m] - initial particle size at t = 0 s
        * αD [ ] - mass accommodation coefficient 
        * αT [ ] - thermal accommodation coefficient
        * isCc [Bool] - whether to correct particle velocity for non-continuum effects 
        * issphereCS [Bool] - whether particle has approximately spherical cross section (for CD calculation) 
    output:
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 

    """
    setup_pμϕ_ice(p,T,𝒮,:Earth,:Earth;rat=rat,r₀=r₀,αD=αD,αT=αT,isCc=isCc,issphereCS=issphereCS)
end

function calcCDsphereCS(Re::Float64)::Float64
    """
    calculate CD for a sphere 
    (or a particle whose relative cross section is approximately circular)
    following Loth (2008) eq (29) 
    [doi:10.1016/j.powtec.2007.06.001])
    from Clift & Gauvin (1970) 
    input:
        * Re [ ] - Reynolds number or effective Reynolds number 
    ouput:
        * CD [ ] - drag coefficient 
    """
    24/Re * (1 + 0.15*Re^(229/333.)) + 0.42/(1 + 4.25e4*Re^(-1.16))
end

function calcCDnonsphereCS(Re::Float64)::Float64
    """
    calculate CD for particles with a non-spherical cross section 
    from Loth (2008) eq (31)
    [doi:10.1016/j.powtec.2007.06.001]
    input:
        * Re [ ] - effective Reynolds number 
    output:
        * CD [ ] - effective drag coefficient 
    """
    24/Re * (1 + 0.035*Re^0.74) + 0.42/(1 + 33/Re^0.5)
end

function calcCshape(rat)
    """
    calculate effect of shape on drag at large Re
    from Loth (2008)
    [doi:10.1016/j.powtec.2007.06.001]
    Astarsurf from eq (21)
    Cshape oblate from eq (25)
    Cshape prolate from eq (26)
    input:
        * rat [ ] - particle axis ratio 
    output:
        * Cshape [ ] - shape correction factor for drag coefficient CD at high Re
    """
    if rat < 1 # oblate spheriod 
        ε = (1 - rat^2)^0.5
        Astarsurf = 0.5*rat^(-2/3.) + 0.25*rat^(4/3.)/ε*log((1 + ε)/(1 - ε)) 
        Cshape = 1 + 1.5*(Astarsurf-1)^0.5 + 6.7*(Astarsurf-1) 
    elseif rat > 1 # prolate spheroid 
        ε2 = (1 - rat^(-2))^0.5
        Astarsurf = 0.5*rat^(-2/3.) + 0.5*rat^(1/3.)/ε2*asin(ε2)
        Cshape = 1. + 0.7*√(Astarsurf-1) + 2.4*(Astarsurf-1)
    else # spherical 
        Cshape = 1.
    end
    Cshape 
end

function calcfshape(rat)
    """
    calculate effect of shape on drag at small Re
    from Loth (2008) Table 1 
    (after comments following eq (28))
    assuming oblate spheoroid (rat<1), prolate spheroid (rat>1), or sphere (rat=1)
    from eq (17), fshape ≡ CD / CD_sphere for Re << 1 and sphere of same volume 
    [doi:10.1016/j.powtec.2007.06.001]
    inputs:
        * rat [ ] - particle axis ratio 
    output:
        * fshape [ ] - shape correction factor for drag coefficient CD at low Re 
    """
    if rat < 1 # oblate, follow f∥
        fshape = (4/3.)*rat^(-1/3.)*(1-rat^2)/(rat + (1-2*rat^2)*acos(rat)/((1 - rat^2)^0.5))
    elseif rat > 1 # prolate, follow f⟂
        fshape = (8/3. * rat^(-1/3.)*(rat^2-1))/(rat + (2*rat^2-3)*log(rat + (rat^2-1)^0.5)/((rat^2-1)^0.5))
    else # sphere, unity ratio
        fshape = 1.
    end
    fshape
end

function calcRestar(Re,Cshape,fshape)
    """
    calculate shape normalized Re 
    from Loth (2008) eq (28b)
    [doi:10.1016/j.powtec.2007.06.001]
    inputs:
        * Re [ ] - Reynolds number 
        * Cshape [ ] - shape correction factor for high Re
        * fshape [ ] - shape correction factor for low Re 
    output:
        Restar [ ] - shape normalized Re  
    """
    Cshape*Re/fshape
end

function calcCdivreq(rat)
    if rat < 1 # oblate spheroid 
        Cdivreq = calcCdivreq_obspheroid(rat)
    elseif rat > 1 # prolate spheroid  
        Cdivreq = calcCdivreq_prospheroid(rat)
    else # sphere
        Cdivreq = calcCdivreq_sphere()
    end
    Cdivreq

end

function calcCdivreq_sphere()
    """
    capacitance of sphere 
    from Pruppacher & Klett (2010) eq (13.77)
    C = r 
    therefore, 
    Cdivreq = 1

    input:
        * nothing
    output:
        * Cdivreq [ ] - capacitance of sphere normalized by equivalent radius (1)
    """
    1.
end

function calcCdivreq_obspheroid(rat)
    """
    capacitance of oblate spheriod 
    from Pruppacher & Klett (2010) eq (13.78)
    C = ε a (arcsin(ε))⁻¹ 
    where a is semi-major axis of spheriod and related to req by
    a = req * (rat)^(-1/3)
    from Loftus & Wordsworth (2021) eq (A1)
    [doi:10.1029/2020JE006653]
    therefore,
    Cdivreq = ε (arcsin(ε))⁻¹ × (rat)^(-1/3)

    inputs:
        * rat [ ] - particle axis ratio 
    output:
        * Cdivreq [ ] - capacitance of oblate spheriod normalized by equivalent radius
    
    """
    ε = (1 - rat^2)^0.5
    ε/(asin(ε)) * rat^(-1/3.)
end

function calcCdivreq_prospheroid(rat)
    """
    capacitance of prolate spheriod 
    from Pruppacher & Klett (2010) eq (13.79)
    (in their notation)
    C = √(a² - b²)/ln[(a + √(a² - b²))/b]
    then by geometry
    rat = a/b & 
    a = req rat^(2/3)
    let
    ε = √(1 - rat⁻²)
    simplifying
    => C = req × rat^(2/3) ε /ln[rat(1 + ε)]
    therefore,
    Cdivreq = rat^(2/3) ε /ln[rat(1 + ε)]

    inputs:
        * rat [ ] - particle axis ratio 
    output:
        * Cdivreq [ ] - capacitance of prolate spheriod normalized by equivalent radius 
    """
    ε = (1 - rat^(-2))^0.5
    rat^(2/3.)*ε/log(rat*(1+ε))
end

function calcAshape(rat)
    """
    calculate particle cross sectional area normalized by cross sectional area of sphere of same mass
    (i.e., Ashape = A(rat)/πreq²)
    for spheroid shaped particles 
    
    input:
        * rat [ ] - particle axis ratio 
    output: 
        * Ashape [ ] - particle cross sectional area normalized by cross sectional area of sphere of same mass
    """
    if rat==1
        Ashape = 1.
    elseif rat<1
        Ashape = rat^(-2/3.)
    elseif rat>1 
        Ashape = rat^(1/3.)
    end
    Ashape 
end

function findv0long(v::Float64,r::Float64,pμϕ::Paramμϕ)::Float64
    """
    calculate expression that evaluates to 0 when v = terminal velocity 
    see SI Text 2 for details 
    inputs:
        * v [m s⁻¹] - velocity 
        * r [m] - equivalent particle radius 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    outputs:
        * Δv [m s⁻¹] - expression that evaluates to 0 when v = terminal velocity 
    """
    # Reynolds number 
    Re = max(2*r*v*pμϕ.ρair/pμϕ.ηair,0.)
    # shape normalized Re 
    Restar = calcRestar(Re,pμϕ.Cshape,pμϕ.fshape)
    # Cunningham correction factor 
    Cc = 1.
    if pμϕ.isCc
        Cc = calcCc(pμϕ.mfp,r)
    end
    # shape normalized CD 
    # following approx shape of particle 
    # from Loth (2008) Fig 8 & 9 
    if pμϕ.issphereCS
        CDstar = calcCDsphereCS(Restar)
    else
        CDstar = calcCDnonsphereCS(Restar)
    end
    CD = CDstar*pμϕ.Cshape/Cc
    v -  (8/3. * (pμϕ.ρc - pμϕ.ρair) / pμϕ.ρair * pμϕ.g / CD * r / pμϕ.Ashape)^(0.5)
end 


function calcvT(r::Float64,pμϕ::Paramμϕ)::Float64
    """
    calculate terminal velocity 
    inputs:
        * r [m] - equivalent particle radius 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    output:
        * vT [m s⁻¹] - terminal velocity 
    """
    # anonymous function for root solver 
    findv0 = (v) -> findv0long(v,r,pμϕ)
    # return 0 if particle is very small (to remove divide by 0 error)
    if r <= 1e-7 
        vT = 0.
    else
        vT = -find_zero(findv0,(1e-9,50.))
    end
    vT
end


function converts2Marsday(t)
    """
    convert seconds into Mars days 
    input:
        t [s] - time in seconds 
    output:
        t [Mars days] - time in Mars days 
    """
    t / sinMarsday
end


function converts2Earthday(t)
    """
    convert seconds into Earth days 
    input:
        t [s] - time in seconds 
    output:
        t [Earth days] - time in Earth days 
    """
    t / sinEarthday
end

function validatepsatice()
    """
    validate implementation of calcpsatice
    """
    println("calculated psat(T = 230 K) = ",calcpsatice(230.), "Pa")
    println("expected psat(T = 230 K) = 8.94735 Pa")
    nothing
end



function calcrconds(t::Float64,pμϕ::Paramμϕ)::Float64
    """
    analytic solution for particle radius from condensational growth 
    assuming spheroid particle with fixed axis ratio 
    from Lohmann et al. (2016) eq (7.29)
    [doi:10.1017/CBO9781139087513]
    with corrections for non-spherical spheroids in SI Text 2
    inputs:
        * t [s] - time since r₀ 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
    output:
        * r [m] - particle radius
    """
    sqrt(pμϕ.r₀^2 + 2. * t * pμϕ.𝒮air * pμϕ.Cdivreq/(pμϕ.ρc*(pμϕ.Fₖ + pμϕ.Fd))) # [m]
end

function isfallH(u::Float64,t,integrator)::Float64
    """
    evaluates to 0 when particle has fallen distance H at time t
    """
    u + integrator.p.H
end

function isfallH(u::Array{Float64,1},t,integrator)::Float64
    """
    evaluates to 0 when particle has fallen distance H at time t
    """
    u[2] + integrator.p.H
end

function calcdzdt(z::Float64,pμϕ::Paramμϕ,t::Float64)::Float64
    """
    calculate dz/dt = vT 
    (call structure from OrdinaryDiffEq)
    inputs:
        * z [m] - height (z(t=0) = 0 m)
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
        * t [s] - time 
    output:
        * vT [m s⁻¹] - terminal velocity 

    """
    r = calcrconds(t,pμϕ)
    calcvT(r,pμϕ)
end

function calcdrzdt!(drzdt::Array{Float64,1},rz::Array{Float64,1},pμϕ::Paramμϕ,t::Float64)::Nothing
    """
    calculate drdt from condensational growth & dz/dt = vT 
    (call structure from OrdinaryDiffEq)
    inputs:
        * drzdt [[m s⁻¹,m s⁻¹]] - array of dr/dt and dz/dt 
        * rz [[m,m]]- array storing r and z values 
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
        * t [s] - time 
    output:
        * nothing (drzdt modified in place)

    """
    r = rz[1]
    Fk′ = pμϕ.Fₖ * pμϕ.K / calcK′(r,pμϕ)
    Fd′ = pμϕ.Fd * pμϕ.D / calcD′(r,pμϕ)
    drzdt[1] = pμϕ.𝒮air/(Fk′ + Fd′)/r/pμϕ.ρc*pμϕ.Cdivreq
    drzdt[2] = calcvT(r,pμϕ)
    nothing 
end

function solvetHfall(pμϕ::Paramμϕ,calctype::Symbol;tmax=3e7,dtmax=1e4,abstol=1e-9,reltol=1e-8,dt=0.1,solver=AutoTsit5(Rosenbrock23()))
    """
    solve differential equation for amount of time for cloud particle to fall 1 scale height 
    while growing via condensation   
    inputs:
        * pμϕ [Paramμϕ] - struct with necessary microphysics and environmnet info 
        * calctype [Symbol] - whether to solve for r(t) from condensational growth  analytically or numerically
            + :ana_rt - use analytical solution for condensational growth 
            + :num_rt - use numerical solution for condensational growth 
        * tmax [s] - maximum amount of time to integrate for 
        * dtmax [s] - maximum Δt allowed for integration
        * abstol [ ] - absolute tolerance for numerical solver (effective 0)
        * reltol [ ] - relative tolerance for numerical solver 
        * dt [s] - initial Δt
        * solver - solver to integrate with (see https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
    outputs:
        * tHfall [s] - time to fall scale height H 
        * sol - ODE solution struct
    """
    
    z₀ = 0. # set initial falling distance to 0 m
    # set up callback to terminate integration when cloud particle
    # has fallen 1 scale height 
    cb = ContinuousCallback(isfallH,terminate!)
    if calctype==:ana_rt
        prob = ODEProblem(calcdzdt, z₀, (0,tmax), pμϕ)
        sol = solve(prob,solver,dt=dt,dtmax=dtmax,callback=cb,abstol=abstol,reltol=reltol)
    elseif calctype==:num_rt
        prob = ODEProblem(calcdrzdt!, [pμϕ.r₀,z₀], (0,tmax), pμϕ)
        sol = solve(prob,solver,dt=dt,dtmax=dtmax,callback=cb,abstol=abstol,reltol=reltol)
    else
        # crash 
        DomainError(calctype,"calctype value is not supported. currently only :ana_rt and :num_rt supported.")
    end
    # only return value of tHfall if solver terminates because of isfallH callback
    tHfall = NaN
    if sol.retcode == ReturnCode.Terminated 
        tHfall = sol.t[end]
    end
    tHfall,sol
end

end