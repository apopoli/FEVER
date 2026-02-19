using Ferrite
using FerriteGmsh
using SparseArrays
using WriteVTK

include("src/Materials.jl")
include("src/CaseIO.jl")
include("src/ThermalBC.jl")

using .Materials
using .CaseIO
using .ThermalBC

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
regmat = region_material_map(case)   # Dict(region_name => Material)

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
outer = getfacetset(grid, "outer")
