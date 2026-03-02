module ElectricConduction

using Ferrite
using LinearAlgebra
using LinearSolve
using NonlinearSolve
using ADTypes
using SparseConnectivityTracer
using SparseMatrixColorings

import ..ElectricBC: build_electric_constraints!
import ..CaseIO: excluded_domains

export solve_conduction_nonlinear!

const EPS0 = 8.8541878128e-12

# --- XLPE conductivity law (E in V/m, beta in mm/kV) ---
E_Vm_to_kVmm(E) = E * 1e-6

function sigma_xlpe(T, E; sigma0, T0, alpha, beta_mm_per_kV)
    return sigma0 * exp(alpha * (T - T0) + beta_mm_per_kV * E_Vm_to_kVmm(E))
end

# --- choose σ(T,E): either from TOML law for XLPE or mat.sigma(T,E) ---
function sigma_for(mat, T, E, cond_cfg::Dict{String,Any})
    model = String(get(cond_cfg, "sigma_model", "material"))
    if model == "xlpe_exp" && mat.name == "XLPE"
        σ0 = Float64(get(cond_cfg, "sigma0", 1e-16))
        T0 = Float64(get(cond_cfg, "T0", 273.15))
        α  = Float64(get(cond_cfg, "alpha", 0.084))
        β  = Float64(get(cond_cfg, "beta_mm_per_kV", 0.0645))
        return sigma_xlpe(T, E; sigma0=σ0, T0=T0, alpha=α, beta_mm_per_kV=β)
    else
        return mat.sigma(T, E)
    end
end

# --- axisymmetric residual assembly for conduction: r_i = ∫ ∇N_i · (σ ∇φ) 2πr dΩ ---
function assemble_conduction_residual!(r, grid, dhφ, dhT, cvφ, uφ_full, Tvec, included_regions, regmat, cond_cfg)
    # ensure r is zeroed (respect element type)
    fill!(r, zero(eltype(r)))

    n = getnbasefunctions(cvφ)
    celldofsφ = zeros(Int, n)
    celldofsT = zeros(Int, n)

    # use element type matching the global r vector (so tracer types are allowed)
    T = eltype(r)
    re = Vector{T}(undef, n)    # local element residuals (generic element type)

    for region in included_regions
        mat = regmat[region]
        cset = getcellset(grid, region)

        for cell in CellIterator(dhφ, cset)
            cid = cellid(cell)
            celldofs!(celldofsφ, dhφ, cid)
            celldofs!(celldofsT, dhT, cid)

            φe = uφ_full[celldofsφ]
            Te = Tvec[celldofsT]

            Ferrite.reinit!(cvφ, cell)
            fill!(re, zero(T))

            cell_coords = getcoordinates(cell)
            for qp in 1:getnquadpoints(cvφ)
                x = spatial_coordinate(cvφ, qp, cell_coords)
                r_axi = x[1]
                w = getdetJdV(cvφ, qp) * (2pi * r_axi)

                ∇φ = function_gradient(cvφ, qp, φe)
                E = sqrt(sum(∇φ .* ∇φ))
                Tqp = function_value(cvφ, qp, Te)

                σ = sigma_for(mat, Tqp, E, cond_cfg)

                for i in 1:n
                    ∇Ni = shape_gradient(cvφ, qp, i)
                    re[i] += (∇Ni ⋅ (σ * ∇φ)) * w
                end
            end

            assemble!(r, celldofsφ, re)
        end
    end

    return nothing
end

# --- NonlinearSolve wrapper on free dofs only ---
struct ConductionData
    grid
    dhφ
    dhT
    cvφ
    Tvec::Vector{Float64}
    included_regions::Vector{String}
    regmat
    cond_cfg::Dict{String,Any}
    chφ
    free::Vector{Int}
    ufull::Vector      # allow different numeric/tracer types
    rfull::Vector
end

# Generic residual that works both for Float64 runs and tracer runs
function residual_free!(rf::AbstractVector, uf::AbstractVector, data::ConductionData)
    # number of global dofs
    nfull = ndofs(data.dhφ)

    # element type we must use for local full vectors (can be Float64 or a tracer type)
    T = eltype(uf)

    # local full solution / residual with element type matching uf
    ufull_local = fill(zero(T), nfull)
    rfull_local = fill(zero(T), nfull)

    # --- Build Dirichlet (prescribed) vector in Float64, then convert to T ---
    # We do this using the existing apply! which expects Float64 arrays in your current code.
    # This avoids calling apply! on tracer arrays directly.
    tmp_prescribed = zeros(Float64, nfull)
    apply!(tmp_prescribed, data.chφ)         # set prescribed values (Float64)
    @inbounds for i in 1:nfull
        # convert Float64 -> tracer type (or Float64->Float64 no-op)
        ufull_local[i] = convert(T, tmp_prescribed[i])
    end

    # --- Put free DOF values from uf into the local full vector ---
    @inbounds for (k, dof) in pairs(data.free)
        ufull_local[dof] = uf[k]
    end

    # --- assemble residual into rfull_local using your existing assembler ---
    # Note: assemble_conduction_residual! must accept the element type T for ufull_local/rfull_local,
    # i.e., it mustn't assume Float64 constants in places where tracer values flow through.
    assemble_conduction_residual!(
        rfull_local,
        data.grid, data.dhφ, data.dhT, data.cvφ,
        ufull_local, data.Tvec,
        data.included_regions, data.regmat,
        data.cond_cfg
    )

    # --- Extract entries corresponding to free DOFs into rf ---
    @inbounds for (k, dof) in pairs(data.free)
        rf[k] = rfull_local[dof]
    end

    return nothing
end


function solve_conduction_nonlinear!(case, grid, regmat, ip, qr, dhT, Tvec, outfile_prefix; Tref=293.15)
    elec = get(get(case.raw, "physics", Dict{String,Any}()), "electric", Dict{String,Any}())
    cond_cfg = get(elec, "conduction", Dict{String,Any}())

    excluded = excluded_domains(case, "electric")
    included_regions = [String(r) for r in keys(regmat) if !(String(r) in excluded)]
    @info "Conduction included regions: $included_regions"

    # union cellsets
    cells = reduce(union, (getcellset(grid, r) for r in included_regions))

    # dofs only on included regions
    dhφ = DofHandler(grid)
    sdhφ = SubDofHandler(dhφ, cells)
    add!(sdhφ, :phi, ip)
    close!(dhφ)

    cvφ = CellValues(qr, ip)

    # constraints
    chφ = build_electric_constraints!(case, grid, dhφ)
    free = Ferrite.free_dofs(chφ)

    # initial guess: solve linear conduction with E=0 (Picard seed)
    K = allocate_matrix(dhφ)
    f = zeros(ndofs(dhφ))
    assembler = start_assemble(K, f)

    # assemble linear K with σ(T,0)
    n = getnbasefunctions(cvφ)
    Ke = zeros(n, n)
    fe = zeros(n)
    celldofsφ = zeros(Int, n)
    celldofsT = zeros(Int, n)

    for region in included_regions
        mat = regmat[region]
        for cell in CellIterator(dhφ, getcellset(grid, region))
            cid = cellid(cell)
            celldofs!(celldofsφ, dhφ, cid)
            celldofs!(celldofsT, dhT, cid)
            Te = Tvec[celldofsT]

            Ferrite.reinit!(cvφ, cell)
            fill!(Ke, 0.0); fill!(fe, 0.0)

            cell_coords = getcoordinates(cell)
            for qp in 1:getnquadpoints(cvφ)
                x = spatial_coordinate(cvφ, qp, cell_coords)
                r_axi = x[1]
                w = getdetJdV(cvφ, qp) * (2pi * r_axi)

                Tqp = function_value(cvφ, qp, Te)
                σ = sigma_for(mat, Tqp, 0.0, cond_cfg)

                for i in 1:n
                    ∇Ni = shape_gradient(cvφ, qp, i)
                    for j in 1:n
                        ∇Nj = shape_gradient(cvφ, qp, j)
                        Ke[i, j] += σ * (∇Ni ⋅ ∇Nj) * w
                    end
                end
            end

            assemble!(assembler, celldofsφ, Ke, fe)
        end
    end

    apply!(K, f, chφ)
    phi0 = K \ f
    apply!(phi0, chφ)

    u0_free = phi0[free]

    # caches still Float64: OK because we will NOT use Duals (finite diff)
    nfull = ndofs(dhφ)
    data = ConductionData(
        grid, dhφ, dhT, cvφ, Tvec, included_regions, regmat, cond_cfg, chφ, free,
        Vector{Any}(undef, nfull), Vector{Any}(undef, nfull)
    )
    # initialize to 0.0 so later fill! works (and types start as Float64)
    fill!(data.ufull, 0.0)
    fill!(data.rfull, 0.0)

    # --------------------------------------------------------------------
    # Build sparse Jacobian *structure* via connectivity tracer + coloring
    # --------------------------------------------------------------------
    # closure for the residual depending only on the free dofs (data captured)
    f_jac = (rf, uf) -> residual_free!(rf, uf, data)

    # a zero-like vector of the same shape as u0_free used by the sparsity routine
    du0 = similar(u0_free)

    # compute jacobian sparsity pattern using connectivity tracer
    jac_sparsity = ADTypes.jacobian_sparsity(
    f_jac, du0, u0_free,
    SparseConnectivityTracer.TracerLocalSparsityDetector()
    )

    # optionally run coloring to reduce number of function evaluations during finite-diff
    # (ADTypes.jacobian_sparsity already returns a prototype with coloring info for many AD
    # backends; keep this here if you want to explicitly request a coloring object:)
    # coloring = SparseMatrixColorings.color(jac_sparsity)
    # NOTE: some workflows accept jac_prototype=coloring instead of plain sparse pattern.

    # build NonlinearFunction using the sparse jacobian prototype (gives solver the sparsity)
    f_nls = NonlinearFunction(
        (rf, uf, p) -> residual_free!(rf, uf, p);
        jac_prototype = jac_sparsity
    )

    prob = NonlinearProblem(f_nls, u0_free, data)

    # Choose Newton solver configuration:
    # - keep autodiff = AutoFiniteDiff() so NonlinearSolve uses finite-diff but respects sparsity
    # - optionally provide a sparse linsolver (requires LinearSolve / KLU / UMFPACK installed)
    #
    # Example with a sparse direct KLU factorization (if you have LinearSolve + KLU):
    # nl_solver = NewtonRaphson(linsolve = LinearSolve.KLUFactorization())
    #
    # If you don't have KLU, use the default NewtonRaphson and NonlinearSolve will use
    # the provided jac_prototype to speed up numerical differentiation.
    nl_solver = NewtonRaphson(; autodiff = AutoFiniteDiff())

    sol = solve(prob, nl_solver)

    # reconstruct full solution
    phi = copy(phi0)
    @inbounds for (k, dof) in pairs(free)
        phi[dof] = sol.u[k]
    end
    apply!(phi, chφ)

    @info "Nonlinear conduction solved. phi_min=$(minimum(phi)) V, phi_max=$(maximum(phi)) V"

    outfile = outfile_prefix * "_phi_cond"
    VTKGridFile(outfile, dhφ) do vtk
        write_solution(vtk, dhφ, phi)
    end
    @info "Wrote $outfile.vtu"

    return dhφ, phi
end

end # module