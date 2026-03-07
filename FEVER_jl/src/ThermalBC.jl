module ThermalBC

using Ferrite

export apply_thermal_bcs!

tofloat(x) = x isa Real ? Float64(x) : error("Expected a number, got $(typeof(x))")

# --- IEC 60287: h(x) using local r  ---
function h_buried_iec60287(; burial_depth::Float64, rho_thermal_soil::Float64)
    return function (x)
        r = x[1]
        @assert r > 0.0 "buried_iec60287: r must be > 0 on the boundary"
        u = burial_depth / r
        RTg = rho_thermal_soil/(2*pi) * log(u + sqrt(u^2 + 1))
        return 1.0 / (2*pi * r * RTg)
    end
end

"""
Legge case.raw["physics"]["thermal"]["bc"] e applica le BC termiche.
Supporta:
- type="convection"         (Robin con h costante)
- type="robin"              (alias di convection)
- type="buried_iec60287"    (Robin con h(x))
- type="dirichlet"          (vincolo essenziale su T)
- type="symmetry"/"adiabatic" (naturale: non fa nulla)

Formato atteso del TOML:
[physics.thermal.bc]
facetset_name = { type = "...", ... }
"""
function apply_thermal_bcs!(assemble_robin_outer!, assembler, case, grid, dh, fv, ch)
    raw = case.raw
    haskey(raw, "physics") || return nothing
    haskey(raw["physics"], "thermal") || return nothing
    haskey(raw["physics"]["thermal"], "bc") || return nothing

    bcs = raw["physics"]["thermal"]["bc"]  # Dict facetset_name => Dict params

    for (fsname, spec_any) in bcs
        spec = spec_any::Dict{String,Any}
        bctype = String(spec["type"])

        haskey(grid.facetsets, fsname) || error(
            "Facetset '$fsname' not found in mesh. Available: $(collect(keys(grid.facetsets)))"
        )

        facetset = getfacetset(grid, fsname)

        if bctype == "convection" || bctype == "robin"
            h = tofloat(spec["h"])
            Tinf = tofloat(spec["Tinf"])
            hfun = x -> h

            assemble_robin_outer!(assembler, dh, fv;
                facetset = facetset,
                hfun = hfun,
                T∞ = Tinf
            )

        elseif bctype == "buried_iec60287"
            Tinf = tofloat(spec["Tinf"])
            d    = tofloat(spec["burial_depth"])
            rho  = tofloat(spec["rho_thermal_soil"])

            hfun = h_buried_iec60287(; burial_depth = d, rho_thermal_soil = rho)

            assemble_robin_outer!(assembler, dh, fv;
                facetset = facetset,
                hfun = hfun,
                T∞ = Tinf
            )

        elseif bctype == "dirichlet"
            haskey(spec, "T") || error("Thermal Dirichlet on '$fsname' requires parameter 'T'")
            Tval = tofloat(spec["T"])

            add!(ch, Dirichlet(:T, facetset, (x, t) -> Tval))

        elseif bctype == "symmetry" || bctype == "adiabatic"
            nothing

        else
            error("Unknown thermal BC type '$bctype' on facetset '$fsname'")
        end
    end

    return nothing
end

end # module