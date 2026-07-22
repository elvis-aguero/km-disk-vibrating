"""
Shared bootstrap included at the top of every script in this directory.

Activates this directory's own environment (`scripts/Project.toml`, which pulls in
`MAT.jl`/`JLD2.jl`/`TOML.jl` — dependencies only the migration/generation scripts need,
kept out of the main package's `Project.toml` so installing `FaradayDisk` itself doesn't
require a MATLAB-file-reading dependency), then `Pkg.develop`s the parent package
in-place if it isn't already resolved. `Pkg.develop(path=...)` is used deliberately
instead of a `[sources]` entry in `Project.toml` — the latter needs a newer Pkg than this
project's `julia = "1.10"` compat floor guarantees; `Pkg.develop` has worked the same way
since early Julia 1.x.
"""
import Pkg

Pkg.activate(@__DIR__)

let manifest = Pkg.project().dependencies
    if !haskey(manifest, "FaradayDisk")
        Pkg.develop(path = joinpath(@__DIR__, ".."))
    end
end

Pkg.instantiate()

using FaradayDisk
