function setup_paths()
%SETUP_PATHS  Put the library folder lib/ on the MATLAB path.
%
%   Called automatically by params.m, so every script that starts with
%       P = params();
%   can find the library functions. Safe to call any number of times.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, 'lib'));
end
