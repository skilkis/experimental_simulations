function [relativePath] = relativepath(path)
%RELATIVEPATH Transforms path into a relative path within the pwd
    if contains(path, pwd)
        relativePath = fullfile('.', erase(path,pwd));
    else
        error('%s not in %s', path, pwd)
    end
end

