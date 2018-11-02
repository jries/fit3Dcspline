function out=contains(str,pattern)
%legacy implementation. Remove for Matlab 2017a or newer
out=any(strfind(str,pattern));
end