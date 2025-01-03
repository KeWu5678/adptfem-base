function hostname = getHostname()
%%GETHOSTNAME returns the hostname on UNIX systems
%   hostname = GETHOSTNAME()

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


    %% READ HOSTNAME
    if isunix()
        [~, hostname] = system('cat /etc/hostname');
        hostname = hostname(1:end-1);
    elseif ispc()
        hostname = 'win';
    elseif ismac()
        hostname = 'mac';
    else
        hostname = 'unknown';
    end        

    %% CHECK RESULT
    if isempty(hostname)
        hostname = 'no_hostname';
    end

end
