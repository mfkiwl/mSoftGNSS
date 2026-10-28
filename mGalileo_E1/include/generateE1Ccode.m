function E1Ccode = generateE1Ccode(PRN,codeType)
% This function generates Galileo E1C CBOC(-) or BOC(1,1) code for PRNs
% 1 to 50. Both code types use 12 samples per primary-code chip.
% The primary memory codes are stored in the file E1c.dat.
%
% E1Ccode = generateE1Ccode(PRN, codeType)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%       codeType    - 'CBOC' for full-band CBOC(-), or 'BOC11' for
%                     narrowband BOC(1,1).
%
%   Outputs:
%       E1Ccode     - row vector containing the E1C waveform, with
%                     4092 * 12 samples per primary-code period.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: generateE1Bcode.m,v 1.1.2.5 2017/04/24 22:00:00 dpl Exp $

%--- Generate the primary code --------------------------------------------
% primary code length
PRIMARY_CODE_LENGTH = 4092;

% Number of Galileo E1C PRN codes
NUM_GAL_PRNS = 50;

% File containing the primary memory codes
fid = fopen('E1c.dat');

% Read all primary code sequences
allprncode = fscanf(fid,'%d', PRIMARY_CODE_LENGTH*NUM_GAL_PRNS);

% Select the primary code sequence for the specified PRN
e1b_primary_code= allprncode(...
    (PRN-1)*PRIMARY_CODE_LENGTH+1:PRN*PRIMARY_CODE_LENGTH);

% Change to bipolar format (-1, +1)
e1cCodeRaw = 1-2*e1b_primary_code;

% Close the fid
fclose(fid);
%--- Select the E1C pilot subcarrier waveform -----------------------------
% 'CBOC' uses CBOC(-); 'BOC11' uses BOC(1,1) on the same 12-point grid.
if strcmp(codeType, 'CBOC')
    plus0 = sqrt(10/11) + sqrt(1/11);
    plus1 = sqrt(10/11) - sqrt(1/11);
    minus0 = -sqrt(10/11) + sqrt(1/11);
    minus1 = -sqrt(10/11) - sqrt(1/11);
elseif strcmp(codeType, 'BOC11')
    plus0 = 1;
    plus1 = 1;
    minus0 = -1;
    minus1 = -1;
end

% Apply the selected 12-point subcarrier template to each primary chip.
boc61Wave = [plus1, plus0, plus1, plus0,plus1, plus0,...
    minus1, minus0, minus1, minus0, minus1, minus0];
cbocWave = repmat(boc61Wave,1,PRIMARY_CODE_LENGTH);
E1Ccode = repelem(e1cCodeRaw,12)' .* cbocWave;
