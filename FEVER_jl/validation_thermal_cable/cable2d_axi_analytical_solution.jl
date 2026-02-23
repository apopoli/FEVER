using Printf

# ---- Geometry (m) ----
const r_cu    = 20e-3
const r_sc1   = 21e-3
const r_xlpe  = 51e-3
const r_sc2   = 52e-3
const r_al    = 54e-3
const r_ext   = 58e-3

# ---- Materials (W/m/K) ----
const k_map = Dict(
    "cu"     => 3.8e2,
    "sc_in"  => 0.29,
    "xlpe"   => 0.29,
    "sc_out" => 0.29,
    "al"     => 2.37e2,
    "cover"  => 0.29
)

# ---- Electrical / source ----
const I = 1533          # A
const rho_e_cu = 1/5.8E7  # ohm*m (constant for benchmark)
const J = I / (pi * r_cu^2)
const qvol_cu = rho_e_cu * J^2           # W/m^3
const Wprime = qvol_cu * pi * r_cu^2     # W/m

# ---- Boundary / soil ----
const Tsoil = 293.15  # K
d = 1.3 # burial depth (m)
uu = d/r_ext
rho_thermal_soil = 1.3 # K*m/W (placeholder)
const RTg = rho_thermal_soil/(2*pi)*log(uu+sqrt(uu^2+1)) # IEC 60287

# If instead you have h_eff: RTg should satisfy RTg = 1/(2π r_ext h_eff)

# ---- Thermal resistances of cable layers (per unit length) ----
"""
Thermal resistance per unit length of a cylindrical layer:
R = (1/(2πk)) * ln(Rout/Rin)   [K*m/W]
"""
Rlayer(k, Rin, Rout) = (1.0 / (2pi*k)) * log(Rout/Rin)

function cable_layer_stack()
    # list from outside to inside (excluding copper bulk generation treatment)
    # We'll treat copper separately, but we still need its surface node at r_cu.
    return [
        ("cover",  r_al,  r_ext),
        ("al",     r_sc2, r_al),
        ("sc_out", r_xlpe, r_sc2),
        ("xlpe",   r_sc1, r_xlpe),
        ("sc_in",  r_cu,  r_sc1),
    ]
end

"""
Compute analytic interface temperatures at each radius, starting from outer boundary.
Returns Dict radius=>temperature.
"""
function analytic_interface_temperatures(; Tsoil=Tsoil, RTg=RTg, Wprime=Wprime, k_map=k_map)
    Tint = Dict{Float64,Float64}()

    # outer boundary temperature from soil resistance:
    T_rext = Tsoil + Wprime * RTg
    Tint[r_ext] = T_rext

    # march inward across each layer: T(Rin) = T(Rout) + W' * Rlayer
    for (name, Rin, Rout) in cable_layer_stack()
        k = k_map[name]
        ΔT = Wprime * Rlayer(k, Rin, Rout)
        Tint[Rin] = Tint[Rout] + ΔT
    end

    # now copper surface temperature is Tint[r_cu]
    # centerline in copper: add parabola contribution
    kcu = k_map["cu"]
    T0 = Tint[r_cu] + (qvol_cu / (4kcu)) * (r_cu^2)
    Tint[0.0] = T0

    return Tint
end

"""
Analytical temperature profile T(r) (K) for 0 <= r <= r_ext.
Piecewise: copper (parabolic) + layers (log).
"""
function T_analytic(r; Tsoil=Tsoil, RTg=RTg, Wprime=Wprime, k_map=k_map)
    @assert 0.0 <= r <= r_ext + 1e-12 "r out of range"

    Tint = analytic_interface_temperatures(; Tsoil, RTg, Wprime, k_map)

    if r <= r_cu
        # copper with generation
        kcu = k_map["cu"]
        return Tint[r_cu] + (qvol_cu/(4kcu))*(r_cu^2 - r^2)
    elseif r <= r_sc1
        k = k_map["sc_in"]; Rout = r_sc1; return Tint[Rout] + (Wprime/(2pi*k))*log(Rout/r)
    elseif r <= r_xlpe
        k = k_map["xlpe"];  Rout = r_xlpe; return Tint[Rout] + (Wprime/(2pi*k))*log(Rout/r)
    elseif r <= r_sc2
        k = k_map["sc_out"]; Rout = r_sc2; return Tint[Rout] + (Wprime/(2pi*k))*log(Rout/r)
    elseif r <= r_al
        k = k_map["al"];    Rout = r_al;  return Tint[Rout] + (Wprime/(2pi*k))*log(Rout/r)
    else
        k = k_map["cover"]; Rout = r_ext; return Tint[Rout] + (Wprime/(2pi*k))*log(Rout/r)
    end
end

# ---- Quick report ----
Tint = analytic_interface_temperatures()
@printf("qvol_cu = %.6g W/m^3\n", qvol_cu)
@printf("W'      = %.6g W/m\n", Wprime)
@printf("T(r_ext)= %.4f K\n", Tint[r_ext])
@printf("T(r_cu) = %.4f K\n", Tint[r_cu])
@printf("T(0)    = %.4f K\n", Tint[0.0])

# Example sampling:
# rs = range(0.0, r_ext; length=50)
# Ts = [T_analytic(r) for r in rs]


using DelimitedFiles
using Printf

# -----------------------------
# Leggi risultati FEM dal CSV (r,z,T) e confronta con analitica
# -----------------------------
femfile = "cable2d_axi_bottom_T.csv"

# legge tutto (header + dati)
raw = readdlm(femfile, ',', Any; header=true)
data = raw[1]  # matrice Any
# header = raw[2]  # se ti serve

# estrai colonne come Float64
r_fem = Float64.(data[:, 1])
z_fem = Float64.(data[:, 2])
T_fem = Float64.(data[:, 3])

plot(r_fem,T_analytic.(r_fem), labels="Analytic",linewidth=5)
xlabel!("r (m)")
xlabel!("r (m)")
scatter!(r_fem,T_fem, labels="FEM", color=:red,ms=2, ma=0.5)
