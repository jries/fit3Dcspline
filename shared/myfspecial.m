function filt = myfspecial(type,varargin)
% Generates normalized filters

switch type
    case 'average'
        S=varargin{1};
        filt=ones(S);
    case 'gaussian'
        S=varargin{1};
        sig=varargin{2};
        S=S-1;
        x = [-S/2:1:S/2];
        x = exp(-x.^2/(2*sig^2));
        filt = x'*x;
end
filt = filt./sum(filt(:));