function [K,sup] = generateKernels(t,sig,varargin)

%%% INPUT ARGUMENTS %%%

%   t: time vector for kernels (must be evenly spaced about zero, e.g. t = linspace(-1000,1000,10000)) (ms)
%   sig: kernel width (ms)

%%% OUTPUT ARGUMENTS %%%

%   K: struct variable containing y-values for each kernel, according to inputs t and sig
%   sup: struct variable defining domain over which kernel functions assume non-zero values

if nargin > 2
    switch varargin{1}
        case 'boxcar'
            K = 1/(2*sig*sqrt(3));
            sup = [-sig*sqrt(3), sig*sqrt(3)];
        case 'triangle'
            K = (1/(6*sig^2))*((sig*sqrt(6)) - abs(t));
            sup = [-sig*sqrt(6), sig*sqrt(6)];
        case 'epan'
            K = (3/(4*sig*sqrt(5)))*(1 - (t.^2/(5*sig^2)));
            sup = [-sig*sqrt(5), sig*sqrt(5)];
        case 'gauss'
            K = (1/(sig*sqrt(2*pi)))*exp(-t.^2/(2*sig^2));
            sup = [-inf, inf];
        otherwise
            error('Unexpected kernel shape')
    end
else

    % Define some standard kernel types K(t,sig)
    K.boxcar = 1/(2*sig*sqrt(3));
    K.triangle = (1/(6*sig^2))*((sig*sqrt(6)) - abs(t));
    K.epan = (3/(4*sig*sqrt(5)))*(1 - (t.^2/(5*sig^2)));
    K.gauss = (1/(sig*sqrt(2*pi)))*exp(-t.^2/(2*sig^2));

    % Define support vectors for kernel types
    % i.e. domain over which kernel functions assume non-zero values
    sup.boxcar = [-sig*sqrt(3), sig*sqrt(3)];
    sup.triangle = [-sig*sqrt(6), sig*sqrt(6)];
    sup.epan = [-sig*sqrt(5), sig*sqrt(5)];
    sup.gauss = [-inf, inf];
end
end