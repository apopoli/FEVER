using Ferrite
using NonlinearSolve
using SparseArrays
using LinearAlgebra
using Plots

# ----------------------------------------------------------------------
# 1. Grid e DofHandler
# ----------------------------------------------------------------------
grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
cellvalues = CellValues(qr, ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# ----------------------------------------------------------------------
# 2. Condizioni al Contorno e Soluzione Analitica
# ----------------------------------------------------------------------
u_exact(x::Vec{2}) = sin(π * x[1])

function source_term(x::Vec{2})
    s = sin(π * x[1])
    c = cos(π * x[1])
    return π^2 * s * (1.0 + 3.0 * π^2 * c^2)
end

∂Ω = union(
    getfacetset(grid, "left"), 
    getfacetset(grid, "right"), 
    getfacetset(grid, "top"), 
    getfacetset(grid, "bottom")
)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ω, (x, t) -> u_exact(x)))
close!(ch)
update!(ch, 0.0)

# ----------------------------------------------------------------------
# 3. Funzione Residuo
# ----------------------------------------------------------------------
function residual!(R::Vector, u::Vector, p::NamedTuple)
    dh, cellvalues, ch = p
    fill!(R, 0.0)
    
    n_basefuncs = getnbasefunctions(cellvalues)
    Re = zeros(n_basefuncs)
    
    for cell in CellIterator(dh)
        Ferrite.reinit!(cellvalues, cell)
        fill!(Re, 0.0)
        
        ue = u[celldofs(cell)]
        coords = getcoordinates(cell)
        
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            x = spatial_coordinate(cellvalues, q_point, coords)
            
            grad_u = zero(Vec{2})
            for j in 1:n_basefuncs
                grad_u += ue[j] * shape_gradient(cellvalues, q_point, j)
            end
            
            # Proprietà non lineare: ε = 1 + |∇u|²
            epsilon = 1.0 + norm(grad_u)^2
            ρ = source_term(x)
            
            for i in 1:n_basefuncs
                δu = shape_value(cellvalues, q_point, i)
                ∇δu = shape_gradient(cellvalues, q_point, i)
                
                Re[i] += (∇δu ⋅ (epsilon * grad_u)) * dΩ - (δu * ρ) * dΩ
            end
        end
        
        for (i, dof) in enumerate(celldofs(cell))
            R[dof] += Re[i]
        end
    end
    
    # Gestione BC - VERSIONE CORRETTA PER FERRITE RECENTE
    # Per i DOF vincolati, impostiamo R[i] = 0.0
    # (u0 già soddisfa le BC, e Newton mantiene questa proprietà)
    for i in ch.dofs
        R[i] = 0.0
    end
    
    return nothing
end

# ----------------------------------------------------------------------
# 4. Setup e Solve
# ----------------------------------------------------------------------
u0 = zeros(ndofs(dh))
apply!(u0, ch)

params = (dh=dh, cellvalues=cellvalues, ch=ch)

prob = NonlinearProblem(residual!, u0, params)
sol = solve(prob, NewtonRaphson(); abstol = 1e-12, reltol = 1e-12)

u_fem = sol.u

# ----------------------------------------------------------------------
# 5. Plot e Confronto
# ----------------------------------------------------------------------
x_plot = []
u_fem_plot = []
u_exact_plot = []

node_to_dof = Dict{Int, Int}()
for cell in CellIterator(dh)
    for (i, dof) in enumerate(celldofs(cell))
        node_to_dof[cell.nodes[i]] = dof
    end
end

for (node_idx, node) in enumerate(grid.nodes)
    if haskey(node_to_dof, node_idx)
        if abs(node[2] - 0.5) < 1e-6
            dof = node_to_dof[node_idx]
            push!(x_plot, node[1])
            push!(u_fem_plot, u_fem[dof])
            push!(u_exact_plot, u_exact(node))
        end
    end
end

p = sortperm(x_plot)
plot(x_plot[p], u_exact_plot[p], label="Analitica", linewidth=2, linestyle=:dash)
plot!(x_plot[p], u_fem_plot[p], label="FEM (Nonlineare)", seriestype=:scatter, markersize=3)
xlabel!("x")
ylabel!("u")
title!("Confronto Soluzione (y=0.5)")
display(plot!)

VTKGridFile("nonlinear_poisson", dh) do vtk
    write_solution(vtk, dh, u_fem)
end