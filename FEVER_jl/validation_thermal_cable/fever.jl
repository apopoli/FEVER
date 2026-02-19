using Ferrite
using FerriteGmsh
using SparseArrays
using WriteVTK

# -----------------------------
# Input mesh Gmsh
# -----------------------------
const mshfile = "cable2d_axi.msh"

# Import grid + Physical Groups -> cellsets/facetsets
grid = FerriteGmsh.togrid(mshfile)

@info "Cellsets:  $(collect(keys(grid.cellsets)))"
@info "Facetsets: $(collect(keys(grid.facetsets)))"

# -----------------------------
# Parametri fisici
# -----------------------------
# Materias, see https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9863530
const k_map = Dict(
    "cu"     => 3.8E2,
    "sc_in"  => 0.29,
    "xlpe"   => 0.29,
    "sc_out" => 0.29,
    "al"     => 2.37E2,
    "cover"  => 0.29
)

# Joule heating nel rame
const I = 1533.097214951819
const rc = 20e-3
const rho_e_cu = 1/5.8E7                 # [Ω m] (placeholder, 20°C)
const J = I / (pi * rc^2)                # [A/m^2]
const qvol_cu = rho_e_cu * J^2           # [W/m^3]

qvol_for_set(name::String) = (name == "cu") ? qvol_cu : 0.0

# h_eff = 1/(2π r_ext RTg).  r_ext lo prendiamo dalla mesh (max r sui nodi dell'outer)
function estimate_r_ext(grid)
    # prende il r massimo tra i nodi
    rmax = -Inf
    for n in getnodes(grid)
        rmax = max(rmax, n.x[1])
    end
    return rmax
end
const r_ext = estimate_r_ext(grid)
const h_eff = 1.0 / (2pi * r_ext * RTg)

@info "Estimated r_ext = $r_ext m, h_eff = $h_eff W/m^2/K"

# Robin equivalente "interramento": -k dT/dn = h_eff (T - Tsoil)
const Tsoil = 293.15  # K (20°C)
# thermal resistance between the buried cable and the soil
d = 1.3 # burial depth (m)
uu = d/r_ext
rho_thermal_soil = 1.3 # K*m/W (placeholder)
const RTg = rho_thermal_soil/(2*pi)*log(uu+sqrt(uu^2+1)) # IEC 60287

# -----------------------------
# Spazio FEM
# -----------------------------
ip  = Lagrange{RefTriangle, 1}() # se la mesh è triangolare
qr  = QuadratureRule{RefTriangle}(2)
fqr = FacetQuadratureRule{RefTriangle}(2)

dh = DofHandler(grid)
add!(dh, :T, ip)
close!(dh)

cv = CellValues(qr, ip)
fv = FacetValues(fqr, ip)

# -----------------------------
# Assembly assialsimmetrico
# -----------------------------
function assemble_volume_for_set!(assembler, dh, cv; cellset, k::Float64, qvol::Float64)
    nbase = getnbasefunctions(cv)
    Ke = zeros(nbase, nbase)
    fe = zeros(nbase)

    for cell in CellIterator(dh, cellset)
        reinit!(cv, cell)
        fill!(Ke, 0.0)
        fill!(fe, 0.0)

        cell_coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, cell_coords)
            r = x[1]
            w = getdetJdV(cv, qp) * (2pi * r)

            for i in 1:nbase
                Ni  = shape_value(cv, qp, i)
                ∇Ni = shape_gradient(cv, qp, i)
                fe[i] += Ni * qvol * w

                for j in 1:nbase
                    ∇Nj = shape_gradient(cv, qp, j)
                    Ke[i, j] += k * (∇Ni ⋅ ∇Nj) * w
                end
            end
        end

        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return nothing
end

function assemble_robin_outer!(assembler, dh, fv; facetset, h::Float64, T∞::Float64)
    nbase = getnbasefunctions(fv)
    Ke = zeros(nbase, nbase)
    fe = zeros(nbase)

    for fc in FacetIterator(dh, facetset)
        reinit!(fv, fc)
        fill!(Ke, 0.0)
        fill!(fe, 0.0)

        facet_coords = getcoordinates(fc)

        for qp in 1:getnquadpoints(fv)
            x = spatial_coordinate(fv, qp, facet_coords)
            r = x[1]
            w = getdetJdV(fv, qp) * (2pi * r)

            for i in 1:nbase
                Ni = shape_value(fv, qp, i)
                fe[i] += (h * T∞) * Ni * w
                for j in 1:nbase
                    Nj = shape_value(fv, qp, j)
                    Ke[i, j] += h * Ni * Nj * w
                end
            end
        end

        assemble!(assembler, celldofs(fc), Ke, fe)
    end
    return nothing
end

# -----------------------------
# Build system
# -----------------------------
K = allocate_matrix(dh)
f = zeros(ndofs(dh))

assembler = start_assemble(K, f)   # <-- SOLO QUI

for name in ["cu", "sc_in", "xlpe", "sc_out", "al", "cover"]
    @assert haskey(grid.cellsets, name)
    assemble_volume_for_set!(assembler, dh, cv;
        cellset = getcellset(grid, name),
        k      = k_map[name],
        qvol   = qvol_for_set(name)
    )
end

outer = getfacetset(grid, "outer")
@info "outer facets = $(length(outer))"
assemble_robin_outer!(assembler, dh, fv; facetset=outer, h=h_eff, T∞=Tsoil)

# -----------------------------
# Constraints (nessun Dirichlet per ora)
# axis/top/bottom sono adiabatici naturali -> non si impone nulla
# -----------------------------
ch = ConstraintHandler(dh)
close!(ch)

apply!(K, f, ch)
T = K \ f
apply!(T, ch)

@info "Solved. Tmin=$(minimum(T)) K, Tmax=$(maximum(T)) K"

# -----------------------------
# Export VTK
# -----------------------------
outfile = mshfile[1:end-4]
VTKGridFile(outfile, dh) do vtk
    write_solution(vtk, dh, T)
end

@info "Wrote $outfile.vtu"
