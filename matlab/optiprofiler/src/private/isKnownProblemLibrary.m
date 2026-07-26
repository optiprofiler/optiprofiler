function known = isKnownProblemLibrary(name)
%ISKNOWNPROBLEMLIBRARY reports whether NAME has an explicit definition.

    definitions = getProblemLibraryDefinitions();
    known = any(strcmp({definitions.name}, name));
end
