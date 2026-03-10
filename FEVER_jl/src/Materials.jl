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
            # E is assumed in V/m
            sigma = (T, E) -> begin
                σ0 = 1e-16
                T0 = 273.15
                α  = 0.084
                βmmkV = 0.0645
                E_kV_per_mm = E * 1e-6          # V/m -> kV/mm
                σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
            end,
            eps_r = T -> 2.3,
        ),
        "xlpe_LIMES" => Material(
            name="XLPE_LIMES",
            aliases=["XLPE_LIMES"],
            k = T -> 0.29,
            # E is assumed in V/m
            sigma = (T, E) -> begin
                σ0 = 1e-17
                T0 = 273.15
                α  = 0.074
                βmmkV = 0.10
                E_kV_per_mm = E * 1e-6          # V/m -> kV/mm
                σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
            end,
            eps_r = T -> 2.3,
        ),
        "LSR" => Material(
            name="LSR",
            aliases=["lsr", "LSR"],
            k = T -> 0.25,
            # E is assumed in V/m
            sigma = (T, E) -> begin
                σ0 = 7.6872e-16
                T0 = 273.15
                α  = 0.0411
                βmmkV = 0.094
                E_kV_per_mm = E * 1e-6          # V/m -> kV/mm
                σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
            end,
            eps_r = T -> 4.1,
        ),
        "LDPE" => Material(
            name="LDPE",
            aliases=["ldpe", "LDPE"],
            k = T -> 0.29,
            # E is assumed in V/m
            sigma = (T, E) -> begin
                σ0 = 3.4797e-17
                T0 = 273.15
                α  = 0.0706
                βmmkV = 0.0966
                E_kV_per_mm = E * 1e-6          # V/m -> kV/mm
                σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
            end,
            eps_r = T -> 2.2,
        ),
        "EVA" => Material(
            name="EVA",
            aliases=["eva", "EVA"],
            k = T -> 0.38,
            # E is assumed in V/m
            sigma = (T, E) -> begin
                σ0 = 2.5714e-17
                T0 = 273.15
                α  = 0.0831
                βmmkV = 0.1291
                E_kV_per_mm = E * 1e-6          # V/m -> kV/mm
                σ0 * exp(α*(T - T0) + βmmkV*E_kV_per_mm)
            end,
            eps_r = T -> 2.2,         
        ),
        "semicon" => Material(
            name="Semiconductive layer",
            aliases=["sc_in","sc_out","semicon","sc"],
            k = T -> 0.34,
            sigma = (T, E) -> 6E3,    
            eps_r = T -> 2.3,        
        ),
        "semicon_AP" => Material(
            name="Semiconductive layer",
            aliases=["semicon_AP"],
            k = T -> 0.29,
            sigma = (T, E) -> 6E-6,  
            eps_r = T -> 2.3,
        ),
        "cover" => Material(
            name="Cover/Jacket",
            aliases=["cover", "jacket"],
            k = T -> 0.29,
            sigma = (T, E) -> 1e-15,  
            eps_r = T -> 2.5,         
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
