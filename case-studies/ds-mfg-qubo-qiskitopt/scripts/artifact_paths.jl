function public_artifact_path(path::AbstractString, root::AbstractString)
    normalized = normpath(abspath(path))
    root_normalized = normpath(abspath(root))
    relative = relpath(normalized, root_normalized)
    if relative == ".." || startswith(relative, "../") || startswith(relative, "..\\")
        return basename(normalized)
    end
    return relative
end
