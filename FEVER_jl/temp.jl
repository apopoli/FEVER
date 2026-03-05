using Ferrite, FerriteGmsh, StaticArrays

function rescale_and_shift_grid(grid::Grid; scale = 1/1000.0, axis = 2)
    # info sul tipo Node
    NodeType = typeof(grid.nodes[1])
    fnames = fieldnames(NodeType)          # campi della struct Node
    # cerco il campo coords (vector) e l'id (Integer)
    coords_field = nothing
    id_field = nothing
    for f in fnames
        v = getfield(grid.nodes[1], f)
        if isa(v, AbstractVector) || isa(v, SVector)
            coords_field = f
        elseif isa(v, Integer)
            id_field = f
        end
    end
    if coords_field === nothing
        error("Non ho trovato un campo coordinate nei nodi (campi: $(fnames)).")
    end
    # estraggo tutte le coordinate scalate e calcolo ymin
    coords = [getfield(n, coords_field) for n in grid.nodes]
    # converti a vettori Float64 (compatibilità SVector/Vector)
    coordsf = [collect(Float64.(c)) for c in coords]
    for i in eachindex(coordsf)
        coordsf[i] .= coordsf[i] .* scale
    end
    ymin = minimum(c[axis] for c in coordsf)
    # applica traslazione
    for c in coordsf
        c[axis] -= ymin
    end

    # costruiamo nuovi nodi: proviamo diverse firme di costruttore
    newnodes = Vector{NodeType}(undef, length(grid.nodes))
    success = false
    # possibile firma 1: NodeType(id, coords)
    try
        for (i, n) in enumerate(grid.nodes)
            idv = id_field === nothing ? i : getfield(n, id_field)
            newcoords = SVector{length(coordsf[i])}(coordsf[i]...)  # prova SVector
            newnodes[i] = NodeType(idv, newcoords)
        end
        success = true
    catch err
        # fallisce; continueremo a provare altre firme
        @warn "Costruttore NodeType(id, coords) fallito: $err"
    end

    if !success
        # possibile firma 2: NodeType(coords, id) (inversa)
        try
            for (i, n) in enumerate(grid.nodes)
                idv = id_field === nothing ? i : getfield(n, id_field)
                newcoords = SVector{length(coordsf[i])}(coordsf[i]...)
                newnodes[i] = NodeType(newcoords, idv)
            end
            success = true
        catch err
            @warn "Costruttore NodeType(coords, id) fallito: $err"
        end
    end

    if !success
        # possibile firma 3: NodeType(...) con keyword args: proviamo a creare via reinterpret/copy (ultimo tentativo)
        @warn "Non sono riuscito a costruire nuovi Node in memoria: farò fallback a riscrittura MSH su disco."
        return nothing  # segnaliamo fallback necessario
    end

    # se qui, abbiamo nuovi nodi; costruiamo una nuova Grid (copia del grid originale ma con nodes nuovi)
    newgrid = deepcopy(grid)   # copia struttura base
    # tentiamo di sostituire il campo `nodes` se è mutabile; altrimenti usiamo setfield!
    try
        setfield!(newgrid, :nodes, newnodes)
    catch err
        # se setfield! non è permesso, proviamo assign diretto (potrebbe non funzionare)
        try
            newgrid.nodes = newnodes
        catch err2
            @warn "Impossibile assegnare nodes al nuovo Grid in memoria: $err2"
            return nothing
        end
    end

    return newgrid
end

# Uso:
mshfile = "joint_LIMES/joint_LIMES.msh"
grid_orig = FerriteGmsh.togrid(mshfile)
grid_fixed = rescale_and_shift_grid(grid_orig; scale=1/1000.0, axis=2)

if grid_fixed === nothing
    # fallback: riscrivo .msh su disco (MSH v2 ASCII) e ricarico
    @info "Eseguo fallback: riscrivo il file .msh scalato/traslato e ricarico"
    mshfile2 = replace(mshfile, r"\.msh$" => "_scaled.msh")
    scale_and_shift_msh(mshfile, mshfile2; scale=1/1000.0)  # usa la funzione che ti ho dato prima
    grid_fixed = FerriteGmsh.togrid(mshfile2)
end

# ora grid_fixed è pronto per l'uso