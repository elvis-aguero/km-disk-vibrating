function A = generate_dtn(nr, D, NameValueArgs)
% GENERATE_DTN  Dense DtN matrix of the domain (nr radial nodes, bath diameter D in disk radii). Mirrors Julia generate_dtn.
%   A = generate_dtn(nr, D)                          single-layer, order 6 (default, O(nr^2))
%   A = generate_dtn(nr, D, 'method', 'polar')       the paper's scheme (DTNVectorized, O(nr^3), slow)
%   A = generate_dtn(nr, D, 'order', 4)              finite-difference order (single_layer only)
% The old script that built and saved the D20Quant20 cache is now generate_dtn_cache_D20.m.
arguments
    nr (1, 1) double
    D (1, 1) double
    NameValueArgs.method (1, 1) string {mustBeMember(NameValueArgs.method, ["single_layer", "polar"])} = "single_layer"
    NameValueArgs.order (1, 1) double = 6
end
if NameValueArgs.method == "single_layer"
    A = generate_dtn_single_layer(nr, D, 'order', NameValueArgs.order);
else
    % DTNVectorized saves a .mat in the current folder as a side effect: run it in a throw-away folder.
    work = tempname();  mkdir(work);  here = pwd;
    cleanup = onCleanup(@() cleanupDir(here, work)); %#ok<NASGU>
    cd(work);
    A = DTNVectorized(nr, D);
end
end

function cleanupDir(here, work)
cd(here);
if isfolder(work), rmdir(work, 's'); end
end
