function [hash] = hashstring(string)
%HASHSTRING Hashes a string using SHA-1
%
%AUTHOR: 
%   Oliver Woodford
%
%USAGE:
%
%   hash = string2hash(string)
%
%IN:
%   string - a string!
%
%OUT:
%   hash - a 64 character string, encoding the SHA-1 hash of string
c = convertStringsToChars(string); % Must be char to convert to uint8
persistent md
if isempty(md)
    md = java.security.MessageDigest.getInstance('SHA-1');
end
hash = upper(...
    sprintf('%2.2x', typecast(md.digest(uint8(c)), 'uint8')')...
);
end
