module ElectricBC

using Ferrite
export build_electric_constraints!

tofloat(x) = x isa Real ? Float64(x) : error("Expected a number, got $(typeof(x))")

"""
Reads case.raw["physics"]["electric"]["bc"] and builds a ConstraintHandler for :phi.

Supported:
- type="dirichlet"  with parameter phi
  Optional parameter facetset="name" (otherwise the TOML key is used as facetset)
- type="symmetry"   (natural Neumann 0 -> do nothing)
"""
function build_electric_constraints!(case, grid, dh)
    ch = ConstraintHandler(dh)

    # Pull BC dict safely
    haskey(case.raw, "physics") || (close!(ch); return ch)
    haskey(case.raw["physics"], "electric") || (close!(ch); return ch)
    bcs = get(case.raw["physics"]["electric"], "bc", Dict{String,Any}())

    for (keyname, spec_any) in bcs
        spec = spec_any::Dict{String,Any}
        bctype = String(spec["type"])

        # allow bc keys != actual facetset name
        fsname = String(get(spec, "facetset", keyname))

        if bctype == "dirichlet"
            phi = tofloat(spec["phi"])
            haskey(grid.facetsets, fsname) || error(
                "Facetset '$fsname' not found. Available: $(collect(keys(grid.facetsets)))"
            )
            add!(ch, Dirichlet(:phi, getfacetset(grid, fsname), (x, t) -> phi))

        elseif bctype == "symmetry"
            # natural BC: nothing
            nothing
        else
            error("Unknown electric BC type '$bctype' (key '$keyname')")
        end
    end

    close!(ch)
    update!(ch, 0.0)  # safe even for time-independent BCs
    return ch
end

end # module