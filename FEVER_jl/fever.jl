using Ferrite
using FerriteGmsh
using SparseArrays
using WriteVTK
using LinearAlgebra
using NonlinearSolve
using ADTypes

include("src/Materials.jl")
include("src/CaseIO.jl")
include("src/ThermalBC.jl")
include("src/ElectricBC.jl")
include("src/ElectricConduction.jl")

using .Materials
using .CaseIO
using .ThermalBC
using .ElectricBC
using .ElectricConduction

const ε0 = 8.8541878128e-12

function estimate_r_ext(grid)
    # prende il r massimo tra i nodi
    rmax = -Inf
    for n in getnodes(grid)
        rmax = max(rmax, n.x[1])
    end
    return rmax
end

# -----------------------------
# Assembly (axisymmetric)
# -----------------------------
function assemble_volume_for_set!(assembler, dh, cv; cellset, k::Float64, qvol::Float64)
    nbase = getnbasefunctions(cv)
    Ke = zeros(nbase, nbase)
    fe = zeros(nbase)

    for cell in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cell)
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
        Ferrite.reinit!(fv, fc)
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

"Return excluded + included region names for electric solve."
function electric_regions(case, regmat)
    excluded = CaseIO.excluded_domains(case, "electric")
    included = [String(r) for r in keys(regmat) if !(String(r) in excluded)]
    @assert !isempty(included) "No included regions for electric solve."
    return excluded, included
end

"Create dhφ restricted to included regions (SubDofHandler)."
function setup_electric_dofs(case, grid, regmat, ip)
    excluded, included_regions = electric_regions(case, regmat)

    @info "Electric excluded regions: $excluded"
    @info "Electric included regions: $included_regions"

    cells = reduce(union, (getcellset(grid, r) for r in included_regions))

    dhφ  = DofHandler(grid)
    sdhφ = SubDofHandler(dhφ, cells)
    add!(sdhφ, :phi, ip)
    close!(dhφ)

    return dhφ, included_regions
end

"Assembly for constant coefficient Laplace: ∫ a ∇Ni·∇Nj 2πr dΩ"
function assemble_laplace_constcoef_for_set!(assembler, dh, cv; cellset, a::Float64)
    nbase = getnbasefunctions(cv)
    Ke = zeros(nbase, nbase)
    fe = zeros(nbase)

    for cell in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cell)
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
                    Ke[i, j] += a * (∇Ni ⋅ ∇Nj) * w
                end
            end
        end

        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return nothing
end

function solve_electrostatic!(case, grid, regmat, ip, qr, outfile; Tref=293.15)
    dhφ, included_regions = setup_electric_dofs(case, grid, regmat, ip)
    cvφ = CellValues(qr, ip)

    Kφ = allocate_matrix(dhφ)
    fφ = zeros(ndofs(dhφ))
    assemblerφ = start_assemble(Kφ, fφ)

    for region in included_regions
        mat = regmat[region]
        eps = ε0 * mat.eps_r(Tref)  # simple ε(Tref)
        assemble_laplace_constcoef_for_set!(assemblerφ, dhφ, cvφ;
            cellset = getcellset(grid, region),
            a = eps
        )
    end

    chφ = ElectricBC.build_electric_constraints!(case, grid, dhφ)

    apply!(Kφ, fφ, chφ)
    phi = Kφ \ fφ
    apply!(phi, chφ)

    @info "Solved electrostatics. phi_min=$(minimum(phi)) V, phi_max=$(maximum(phi)) V"
    return dhφ, phi
end

function solve_conduction!(case, grid, regmat, ip, qr, dhT, T, outfile; Tref=293.15)
    elec = get(get(case.raw, "physics", Dict{String,Any}()), "electric", Dict{String,Any}())
    cond = get(elec, "conduction", Dict{String,Any}())
    nonlinear = Bool(get(cond, "nonlinear", false))

    if nonlinear
        return ElectricConduction.solve_conduction_nonlinear!(case, grid, regmat, ip, qr, dhT, T, outfile; Tref=Tref)
    else
        error("Linear conduction not implemented yet (set nonlinear=true or implement the linear branch).")
    end
end

function main()
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

    r_ext = estimate_r_ext(grid)

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
    # Build system
    # -----------------------------
    K = allocate_matrix(dh)
    f = zeros(ndofs(dh))

    assembler = start_assemble(K, f)

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
    # Export VTK - thermal
    # -----------------------------
    outfile = mshfile[1:end-4]
    VTKGridFile(outfile, dh) do vtk
        write_solution(vtk, dh, T)
    end

    @info "Wrote $outfile.vtu"
    # outer = getfacetset(grid, "outer") # utile per fare plot lungo una lines

    # ---

    # -----------------------------
    # ELECTROSTATICS / CONDUCTION
    # -----------------------------

    # --- Dispatch (MUST be after function definitions in a script) ---
    elec = get(get(case.raw, "physics", Dict{String,Any}()), "electric", Dict{String,Any}())
    etype = String(get(elec, "type", "electrostatic"))

    if etype == "electrostatic"
        dhφ, phi = solve_electrostatic!(case, grid, regmat, ip, qr, outfile; Tref=Tref)
    elseif etype == "conduction"
        dhφ, phi = solve_conduction!(case, grid, regmat, ip, qr, dh, T, outfile; Tref=Tref)
    else
        error("Unknown [physics.electric].type = $etype")
    end

    # --- Export once (choose suffix based on mode) ---
    suffix = (etype == "electrostatic") ? "_phi" : "_phi_cond"
    outfile_phi = outfile * suffix

    VTKGridFile(outfile_phi, dhφ) do vtk
        write_solution(vtk, dhφ, phi)
    end
    @info "Wrote $outfile_phi.vtu"
end

main()