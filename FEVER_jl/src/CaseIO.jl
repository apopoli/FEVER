module CaseIO
using TOML
import ..Materials: material_db, canonical_material_key

export load_case, region_material_map

struct Case
    meshfile::String
    region_to_material_key::Dict{String,String}
    raw::Dict{String,Any}
end

function load_case(path::AbstractString)
    raw = TOML.parsefile(path)

    meshfile = raw["mesh"]["file"]
    region_to_material_key = Dict{String,String}(raw["materials"])

    return Case(meshfile, region_to_material_key, raw)
end

"Return region_name => Material object"
function region_material_map(case::Case)
    db = material_db()
    out = Dict{String,Any}()
    for (region, key) in case.region_to_material_key
        ck = canonical_material_key(key, db)
        out[region] = db[ck]
    end
    return out
end

end # module
