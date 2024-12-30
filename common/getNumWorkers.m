function nWorkers = getNumWorkers()
%%GETNUMWORKERS returns the number of workers for parallel computing (1 in
%case of serial computing)
%   nWorkers = GETNUMWORKERS()

% Copyright 2021 Philipp Bringmann
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


    if isOctave()
        nWorkers = 1;
    elseif isempty(ver('parallel'))
        nWorkers = 1;
    else
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            nWorkers = poolobj.NumWorkers;
        else
            nWorkers = 1;
        end
    end
end
