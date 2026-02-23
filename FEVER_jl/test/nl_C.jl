using Ferrite, SparseArrays, LinearAlgebra
using NonlinearSolve
using Plots

# -----------------------------
# Problem definition (manufactured solution)
# -----------------------------
α = 1.0

# Domain Ω = [-1, 1]×[-1, 1]
left  = Vec((-1.0, -1.0))
right = Vec(( 1.0,  1.0))

φ_exact(x::Vec{2}) = sin(pi * (x[1] + 1) / 2)

# ε_r(φ) and derivative
εr(φ)  = 1 + α * φ^2
dεr(φ) = 2α * φ

# Manufactured RHS f(x) for φ_exact on [-1,1]:
# For φ(x)=sin(π(x+1)/2), one gets
# f(x) = (π^2/4) * (1 - 2α) * φ + (3α π^2/4) * φ^3
function rhs_f(x::Vec{2})
    φ = φ_exact(x)
    return (pi^2 / 4) * (1 - 2α) * φ + (3α * pi^2 / 4) * φ^3
end

# -----------------------------
# Mesh / FE space
# -----------------------------
grid = generate_grid(Quadrilateral, (20, 20), left, right)

ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# -----------------------------
# Dirichlet boundary conditions: φ = φ_exact on ∂Ω
# -----------------------------
∂Ω = union(
    getfacetset(grid, "left"),
    getfacetset(grid, "right"),
    getfacetset(grid, "top"),
    getfacetset(grid, "bottom"),
)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ω, (x, t) -> φ_exact(x)))
close!(ch)
update!(ch, 0.0)

# Build a full vector containing only the prescribed Dirichlet values
u_bc = zeros(ndofs(dh))
apply!(u_bc, ch)  # set prescribed dofs in u_bc (common Ferrite pattern):contentReference[oaicite:2]{index=2}

# We solve only for free dofs (Dirichlet dofs eliminated)
free = Ferrite.free_dofs(ch)  # free dofs helper exists in Ferrite examples:contentReference[oaicite:3]{index=3}
nfree = length(free)

# Map global dof -> reduced dof index (0 if constrained)
g2r = zeros(Int, ndofs(dh))
for (k, I) in enumerate(free)
    g2r[I] = k
end

# Workspace to embed the reduced unknowns into a full dof vector
u_full = similar(u_bc)
Ftmp   = zeros(nfree)

# -----------------------------
# Assembly into reduced residual/Jacobian
# -----------------------------
function assemble_reduced!(F::Vector{Float64}, J::Union{Nothing,Matrix{Float64}},
                           u_free::Vector{Float64})
    fill!(F, 0.0)
    if J !== nothing
        fill!(J, 0.0)
    end

    # embed reduced unknowns into full vector with Dirichlet values
    u_full .= u_bc
    @inbounds for (k, I) in enumerate(free)
        u_full[I] = u_free[k]
    end

    nbf = getnbasefunctions(cv)

    for cell in CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        dofs = celldofs(cell)
        ue = @view u_full[dofs]  # element dofs (includes boundary values if present)

        cell_coords = getcoordinates(cell)  # coordinates for spatial_coordinate usage
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)

            # global coordinate of this quadrature point
            xqp = spatial_coordinate(cv, qp, cell_coords)  # same idea as in facet example:contentReference[oaicite:4]{index=4}

            # evaluate u and ∇u at qp (manual reconstruction)
            u_q = 0.0
            ∇u_q = zero(shape_gradient(cv, qp, 1))
            @inbounds for a in 1:nbf
                Na = shape_value(cv, qp, a)
                ∇Na = shape_gradient(cv, qp, a)
                u_q  += Na * ue[a]
                ∇u_q += ∇Na * ue[a]
            end

            ε  = εr(u_q)
            dε = dεr(u_q)
            f  = rhs_f(xqp)

            @inbounds for i in 1:nbf
                I = dofs[i]
                ri = g2r[I]
                ri == 0 && continue  # skip constrained test dofs

                Ni  = shape_value(cv, qp, i)
                ∇Ni = shape_gradient(cv, qp, i)
                gdot = (∇u_q ⋅ ∇Ni)

                # residual: ∫ ε(u) ∇u·∇δu - ∫ f δu
                F[ri] += (ε * gdot - f * Ni) * dΩ

                if J !== nothing
                    @inbounds for j in 1:nbf
                        Jg = dofs[j]
                        rj = g2r[Jg]
                        rj == 0 && continue  # skip constrained trial dofs

                        Nj  = shape_value(cv, qp, j)
                        ∇Nj = shape_gradient(cv, qp, j)

                        # tangent:
                        # ∂/∂u_j [ ε(u) ∇u·∇Ni ] =
                        #   ε(u) (∇Nj·∇Ni) + (dε/du) Nj (∇u·∇Ni)
                        J[ri, rj] += (ε * (∇Ni ⋅ ∇Nj) + dε * Nj * gdot) * dΩ
                    end
                end
            end
        end
    end

    return nothing
end

# -----------------------------
# NonlinearSolve interface
# -----------------------------
function f!(F, u_free, p)
    assemble_reduced!(F, nothing, u_free)
    return nothing
end

function jac!(J, u_free, p)
    assemble_reduced!(Ftmp, J, u_free)
    return nothing
end

Jproto = zeros(nfree, nfree)
nlf = NonlinearFunction(f!; jac = jac!, jac_prototype = Jproto)  # NonlinearFunction supports f!(du,u,p), jac(J,u,p) :contentReference[oaicite:5]{index=5}
u0_free = zeros(nfree)

prob = NonlinearProblem(nlf, u0_free, nothing)                   # NonlinearProblem(f,u0,p) :contentReference[oaicite:6]{index=6}
sol = solve(prob, NewtonRaphson(); abstol = 1e-12, reltol = 1e-12)  # NewtonRaphson solver :contentReference[oaicite:7]{index=7}

u_free = sol.u

# Reconstruct full solution vector (for export/evaluation)
u = copy(u_bc)
@inbounds for (k, I) in enumerate(free)
    u[I] = u_free[k]
end

# -----------------------------
# Compare FEM vs analytical on the top edge (y = 1)
# -----------------------------
points = [Vec((x, 1.0)) for x in range(-1.0, 1.0, length = 201)]
ph = PointEvalHandler(grid, points)
u_fem = evaluate_at_points(ph, dh, u, :u)  # point evaluation utility:contentReference[oaicite:8]{index=8}
u_an  = φ_exact.(points)

xs = getindex.(points, 1)
plot(xs, u_an, label = "analytical", xlabel = "x", ylabel = "ϕ")
plot!(xs, u_fem, label = "FEM", lw = 2)
