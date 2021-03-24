function [hash] = hashstring(string)
%HASHSTRING Hashes a string using SHA-1
    c = convertStringsToChars(string); % Must be char to convert to uint8
    hasher = System.Security.Cryptography.SHA1Managed; % .NET SHA-1 hasher
    sha1 = dec2hex(uint8(hasher.ComputeHash(uint8(c))));
    hash = convertCharsToStrings(sha1');
end
