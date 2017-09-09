function y = bla(x1,varargin)
switch nargin
    case 1
        x2 = 2;
        x3 = 3;
    case 2
        x2 = varargin{1};
        x3 = 3;
    case 3
        x2 = varargin{1};
        x3 = varargin{2};
end
y = x1*x2*x3;
end