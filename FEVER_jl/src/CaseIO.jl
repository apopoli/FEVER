module CaseIO
using TOML
import ..Materials: material_db, canonical_material_key

export load_case, region_material_map, physics_cfg, excluded_domains, physics_bc

struct Case
    meshfile::String
    region_to_material_key::Dict{String,String}
    raw::Dict{String,Any}
end

function load_case(path::AbstractString)
    raw = TOML.parsefile(path)
    casedir = dirname(abspath(path))

    # Mesh path: relativo al case.toml
    meshfile_rel = raw["mesh"]["file"]
    meshfile = isabspath(meshfile_rel) ? meshfile_rel : joinpath(casedir, meshfile_rel)

    # Materials mapping (region => material key)
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

"Return the dictionary under [physics.<phys>], or an empty Dict if missing."
function physics_cfg(case::Case, phys::AbstractString)
    phys_all = get(case.raw, "physics", Dict{String,Any}())
    return get(phys_all, String(phys), Dict{String,Any}())
end

"Return exclude_domain list under [physics.<phys>], defaulting to []."
function excluded_domains(case::Case, phys::AbstractString)
    phys_all = get(case.raw, "physics", Dict{String,Any}())
    cfg = get(phys_all, String(phys), Dict{String,Any}())
    ex = get(cfg, "exclude_domain", Any[])
    return [String(x) for x in ex]
end

"Return boundary condition dict under [physics.<phys>.bc], defaulting to empty."
function physics_bc(case::Case, phys::AbstractString)
    cfg = physics_cfg(case, phys)
    return get(cfg, "bc", Dict{String,Any}())
end

end # module
