function hash = hashdir(directory, varargin)
%HASHDIR Wrapper around util.datahash to allow hashing of directories
%   See documentation of datahash for optional arguments

if exist(directory, 'dir')
    x = dir(directory);
    hash = util.datahash(x, varargin{:});
else
    error('"%s" is not a valid directory', directory)
end
