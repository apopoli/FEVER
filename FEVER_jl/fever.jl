using Ferrite
using FerriteGmsh
using SparseArrays
using WriteVTK

include("src/Materials.jl")
include("src/CaseIO.jl")
include("src/ThermalBC.jl")
include("src/ElectricBC.jl")

using .Materials
using .CaseIO
using .ThermalBC
using .ElectricBC

# -----------------------------
# Input mesh from Gmsh
# -----------------------------

casepath = joinpath(@__DIR__, "case.toml")
case = CaseIO.load_case(casepath)

mshfile = case.meshfile
grid = FerriteGmsh.togrid(mshfile)

@info "Case regions:  $(collect(keys(grid.cellsets)))"
@info "Facetsets: $(collect(keys(grid.facetsets)))"

# -----------------------------
# Materias, see https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9863530
# -----------------------------
regmat = CaseIO.region_material_map(case)   # Dict(region_name => Material)

function estimate_r_ext(grid)
    # prende il r massimo tra i nodi
    rmax = -Inf
    for n in getnodes(grid)
        rmax = max(rmax, n.x[1])
    end
    return rmax
end
const r_ext = estimate_r_ext(grid)
@info "Estimated r_ext = $r_ext m"

# -----------------------------
# FEM space
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
# Assembly (axisymmetric)
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

function assemble_robin_outer!(assembler, dh, fv; facetset, hfun, T∞::Float64)
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
            w = getdetJdV(fv, qp) * (2*pi * r)

            h = hfun(x)

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

assembler = start_assemble(K, f)

function qvol_for_region(case, region, mat; Tref=293.15)
    srcs = get(get(get(case.raw, "physics", Dict()), "thermal", Dict()), "sources", Dict())
    haskey(srcs, region) || return 0.0

    spec = srcs[region]::Dict{String,Any}
    @info "Source spec for region=$region: $spec"

    typ = String(spec["type"])

    if typ == "constant"
        return Float64(spec["qvol"])

    elseif typ == "joule_J"
        J = Float64(spec["J"])
        sigma = mat.sigma(Tref, 0.0)
        return (1.0 / sigma) * J^2

    else
        error("Unknown source type '$typ' for region '$region'")
    end
end

Tref = 293.15
for (region, mat) in regmat
    @assert haskey(grid.cellsets, region) "Mesh missing cellset '$region'"
    qvol = qvol_for_region(case, region, mat; Tref=Tref)

    assemble_volume_for_set!(assembler, dh, cv;
        cellset = getcellset(grid, region),
        k      = mat.k(Tref),
        qvol   = qvol
    )
end

for (region, mat) in regmat
    @info "region=$region k=$(mat.k(293.15)) sigma=$(try mat.sigma(293.15,0.0) catch e e end)"
end

ThermalBC.apply_thermal_bcs!(assemble_robin_outer!, assembler, case, grid, dh, fv)
@info "thermal bc keys in case: " keys(case.raw["physics"]["thermal"]["bc"])
@info "outer facets = " length(getfacetset(grid, "outer"))

# -----------------------------
# Constraints
# axis/top/bottom adiabatic -> no forcing
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
# outer = getfacetset(grid, "outer") # utile per fare plot lungo una lines













# -----------------------------
# Electrostatics (axisymmetric Laplace) on Ω \ exclude_domain
# -----------------------------
# We solve: ∇·(ε ∇phi) = 0   (in axisym -> same 2πr weight as thermal)
# with Dirichlet on interface_cu and outer, symmetry on axis (natural).

# --- helper: union of included cellsets ---
function union_cellsets(grid, region_names::Vector{String})
    @assert !isempty(region_names) "No included regions for electric solve."
    sets = [getcellset(grid, r) for r in region_names]
    return reduce(union, sets)
end

# Regions included in electric
excluded = CaseIO.excluded_domains(case, "electric")
included_regions = [String(r) for r in keys(regmat) if !(String(r) in excluded)]

@info "Electric excluded regions: $excluded"
@info "Electric included regions: $included_regions"

cells_elec = union_cellsets(grid, included_regions)

# --- Dofs only on included regions (critical to avoid singular matrix) ---
dhφ = DofHandler(grid)
sdhφ = SubDofHandler(dhφ, cells_elec)
add!(sdhφ, :phi, ip)
close!(dhφ)

cvφ = CellValues(qr, ip)
fvφ = FacetValues(fqr, ip)  # not strictly needed for Dirichlet, but ok to keep

# --- Assembly function (like thermal, but no source term) ---
function assemble_laplace_for_set!(assembler, dh, cv; cellset, eps::Float64)
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
                ∇Ni = shape_gradient(cv, qp, i)
                for j in 1:nbase
                    ∇Nj = shape_gradient(cv, qp, j)
                    Ke[i, j] += eps * (∇Ni ⋅ ∇Nj) * w
                end
            end
        end

        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return nothing
end

Kφ = allocate_matrix(dhφ)
fφ = zeros(ndofs(dhφ))
assemblerφ = start_assemble(Kφ, fφ)

# simple permittivity: ε = ε0 * eps_r(Tref) per region (refine later)
const ε0 = 8.8541878128e-12

for region in included_regions
    mat = regmat[region]
    eps = ε0 * mat.eps_r(Tref)  # or use temperature-dependent choice later
    assemble_laplace_for_set!(assemblerφ, dhφ, cvφ;
        cellset = getcellset(grid, region),
        eps = eps
    )
end

# Dirichlet constraints from TOML
chφ = ElectricBC.build_electric_constraints!(case, grid, dhφ)

apply!(Kφ, fφ, chφ)
phi = Kφ \ fφ
apply!(phi, chφ)

@info "Solved electrostatics. phi_min=$(minimum(phi)) V, phi_max=$(maximum(phi)) V"

# Export
outfile_phi = outfile * "_phi"
VTKGridFile(outfile_phi, dhφ) do vtk
    write_solution(vtk, dhφ, phi)
end
@info "Wrote $outfile_phi.vtu"



