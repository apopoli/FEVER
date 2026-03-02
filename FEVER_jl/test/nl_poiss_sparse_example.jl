import NonlinearSolve as NLS
import LinearAlgebra
import SparseArrays
import LinearSolve as LS
import SparseConnectivityTracer, SparseMatrixColorings
import ADTypes
import Plots

# --- problem parameters -----------------------------------------------------
L = 0.02                # domain size (m) = 2 cm
N = 50                  # grid points per dimension (tweak for accuracy/cost)
dx = L / (N - 1)
dy = dx

V_left  = 50e3          # 1000 V on left boundary
V_right = 0.0           # 0 V on right boundary

# Temperature: user said "uniform T = 50 deg" -> treat as 50 °C
T_C = 50
T_K = T_C + 273.15

# Conductivity function (as provided by user)
sigma = (T, E) -> begin
    σ0 = 1e-16
    T0 = 273.15
    α  = 0.084
    βmmkV = 0.0645
    # convert E [V/m] -> kV/mm: 1 V/m = 1e-6 kV/mm
    E_kV_per_mm = E * 1e-6
    σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
end

# --- helpers / geometry ----------------------------------------------------
# boolean masks for Dirichlet/Neumann boundaries
is_dirichlet = falses(N, N)
dirichlet_value = zeros(N, N)

# left boundary (i=1) Dirichlet = 1000 V
for j in 1:N
    is_dirichlet[1, j] = true
    dirichlet_value[1, j] = V_left
end
# right boundary (i=N) Dirichlet = 0 V
for j in 1:N
    is_dirichlet[N, j] = true
    dirichlet_value[N, j] = V_right
end
# top and bottom (j=1 and j=N) are Neumann -> leave is_dirichlet false

# initial guess: linear ramp left->right
function init_guess()
    u = zeros(N, N)
    for j in 1:N, i in 1:N
        x = (i-1)*dx
        u[i,j] = V_left*(1.0 - x/L) + V_right*(x/L)  # linear interpolation
        if is_dirichlet[i,j]
            u[i,j] = dirichlet_value[i,j]
        end
    end
    u
end

# clamp helper for indices (not periodic)
inbounds(i) = clamp(i, 1, N)

# --- residual function for NonlinearSolve ----------------------------------
function resid!(F, u, p)
    T = p
    Tu = eltype(u)

    Ex = zeros(Tu, N, N)
    Ey = zeros(Tu, N, N)
    fx = zeros(Tu, N+1, N)
    fy = zeros(Tu, N, N+1)

    for j in 1:N, i in 1:N
        im = max(i-1, 1); ip = min(i+1, N)
        jm = max(j-1, 1); jp = min(j+1, N)
        Ex[i,j] = (u[ip, j] - u[im, j]) / (2*dx)
        Ey[i,j] = (u[i, jp] - u[i, jm]) / (2*dy)
    end

    Ecenter = sqrt.(Ex.^2 .+ Ey.^2)

    σc = sigma.(T, Ecenter)

    for j in 1:N, i in 1:(N-1)
        σface = 0.5*(σc[i,j] + σc[i+1,j])
        fx[i+1, j] = -σface * (u[i+1, j] - u[i, j]) / dx
    end
    for i in 1:N, j in 1:(N-1)
        σface = 0.5*(σc[i,j] + σc[i,j+1])
        fy[i, j+1] = -σface * (u[i, j+1] - u[i, j]) / dy
    end

    for j in 1:N, i in 1:N
        if is_dirichlet[i,j]
            F[i,j] = u[i,j] - dirichlet_value[i,j]   # ok: Float64 promotes to Dual
        else
            div = (fx[i+1, j] - fx[i, j]) / dx + (fy[i, j+1] - fy[i, j]) / dy
            F[i,j] = div
        end
    end
    return nothing
end

# --- build problem & solve -------------------------------------------------
u0  = init_guess()
du0 = similar(u0)

# freeze p = T_K for the sparsity detection call:
f! = (F, u) -> resid!(F, u, T_K)

jac_sparsity = ADTypes.jacobian_sparsity(
    f!, du0, u0,
    SparseConnectivityTracer.TracerSparsityDetector()
)

ff = NLS.NonlinearFunction(resid!; jac_prototype = jac_sparsity)
prob_sp = NLS.NonlinearProblem(ff, u0, T_K; abstol=1e-10, reltol=1e-10)

sol = NLS.solve(
    prob_sp,
    NLS.NewtonRaphson(linsolve = LS.KLUFactorization());
    verbose = true
)

# solution in sol.u
U = sol.u # N×N array of potentials

sol.stats

x = range(0, L, length=N)
y = range(0, L, length=N)

Plots.heatmap(x, y, sol.u'; xlabel="x (m)", ylabel="y (m)", title="Potential (V)", aspect_ratio=1)

Plots.plot(x,sol.u[:,50])