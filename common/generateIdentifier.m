function [foldername, filename] = generateIdentifier(S, root)
%%GENERATEIDENTIFIER creates two strings characterising the solution structure
%array S
%   [foldername, filename] = GENERATEIDENTIFIER(S) creates names for the folder
%   and the file
%
%   [foldername, filename] = GENERATEIDENTIFIER(S, root) creates the folder name
%   inside the given root directory
%

% Copyright 2020 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

    if nargin < 2
        root = '';
    else
        if numel(root) > 0 && ~strcmp(root(end), '/')
            root = [root, '/'];
        end
    end

    % Create string characterising the refinement method
    if S.theta == 1
        refinement = 'unif';
    elseif strcmp(S.method, 'Approx')
        refinement = ['rho', num2str(100*S.rho)];
    elseif strcmp(S.refinementIndicator, 'sep')
        refinement = ['adap', num2str(100*S.theta), '_',...
                      S.refinementIndicator, kappa2str(S.kappa)];
    else
        refinement = ['adap', num2str(100*S.theta), '_',...
                      S.refinementIndicator];
    end

    % Create folder and file name
    if strcmp(S.method, 'Approx')
        foldername = [root, 'Approx_', S.name];
    else
        foldername = [root, S.problem, '_', S.name, '_', S.method];
    end
    filename = [refinement, '_', num2str(S.minNdof), 'dofs_', ...
                S.identifier];

end


function s = kappa2str(kappa)
%%KAPPA2STR converts kappa parameter to a very short string
    exponent = floor(log10(kappa));
    coefficient = kappa / 10^exponent;
    if exponent == 0
        s = [num2str(coefficient)];
    else
        s = [num2str(coefficient), 'e', num2str(exponent)];
    end
end
