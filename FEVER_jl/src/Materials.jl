module Materials

export Material, material_db, canonical_material_key

"""
Material properties are callables so they can be constant or state-dependent.
- k(T)                  thermal conductivity [W/m/K]
- sigma(T, E)           electrical conductivity [S/m]
- eps_r(T)              relative permittivity [-] (optional, for electrostatics)
"""
Base.@kwdef struct Material
    name::String
    k::Function                 = T -> throw(ArgumentError("k(T) not defined for $name"))
    sigma::Function             = (T, E) -> throw(ArgumentError("sigma(T,E) not defined for $name"))
    eps_r::Function             = T -> 1.0
    rho::Function               = T -> NaN       # optional
    cp::Function                = T -> NaN       # optional
    aliases::Vector{String}     = String[]
end

"Return database: canonical_key => Material"
function material_db()
    Dict{String,Material}(
        "copper" => Material(
            name="Copper",
            aliases=["cu", "Cu", "copper"],
            k = T -> 3.8e2,
            sigma = (T, E) -> 5.8e7,
            eps_r = T -> 1.0,
        ),
        "aluminum" => Material(
            name="Aluminum",
            aliases=["al", "Al", "aluminum"],
            k = T -> 2.37e2,
            sigma = (T, E) -> 3.5e7,
            eps_r = T -> 1.0,
        ),
        "xlpe" => Material(
            name="XLPE",
            aliases=["xlpe", "XLPE"],
            k = T -> 0.29,
            sigma = (T, E) -> 1e-14,  # placeholder
            eps_r = T -> 2.3,         # placeholder
        ),
        "semicon" => Material(
            name="Semiconductive layer",
            aliases=["sc_in","sc_out","semicon","sc"],
            k = T -> 0.29,
            sigma = (T, E) -> 1.0,    # placeholder
            eps_r = T -> 10.0,        # placeholder
        ),
        "cover" => Material(
            name="Cover/Jacket",
            aliases=["cover", "jacket"],
            k = T -> 0.29,
            sigma = (T, E) -> 1e-15,  # placeholder
            eps_r = T -> 2.5,         # placeholder
        ),
    )
end

"Normalize a user/mesh string to a canonical material key."
function canonical_material_key(key::AbstractString, db=material_db())
    s = lowercase(strip(String(key)))
    # direct hit
    if haskey(db, s)
        return s
    end
    # search aliases
    for (ck, mat) in db
        if any(a -> lowercase(a) == s, mat.aliases)
            return ck
        end
    end
    throw(ArgumentError("Unknown material key '$key'. Known: $(collect(keys(db)))"))
end

end # module
