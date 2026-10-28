function [eph, firstSubFrame,TOW] = NAVdecoding(I_P)

%Find a CRC-valid Galileo F/NAV page boundary and decode ephemeris data.
%Preamble spacing and convolutional decoding are checked before CRC-24Q.
%
%[eph, firstSubFrame, TOW] = NAVdecoding(I_P)
%
%   Inputs:
%       I_P             - prompt-I tracking outputs at 1 ms intervals.
%
%   Outputs:
%       eph             - decoded satellite ephemeris and message flags.
%       firstSubFrame   - one-based page-start index in I_P, or inf when no
%                         valid candidate is found.
%       TOW             - decoded GST time-of-week reference [s], or inf
%                         when no valid candidate is found.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Jakob Almqvist, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------
%
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
%--- Initialize ephemeris structure  --------------------------------------
% Initialize ephemeris fields before decoding the available navigation pages.
eph = eph_structure_init();

% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 0;

% Starting positions of the first message in the input bit stream
firstSubFrame = inf;

% TOW of the first message
TOW = inf;

%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Rate-1/2, constraint-length-7 Galileo F/NAV convolutional-code trellis.
% Note that the difference from GPS is that the second branch G2 is
% inverted at the end (see ICD)
trellis = poly2trellis(7,[171 133]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%--- Generate the sync pattern --------------------------------------------
% Secondary code is "842E9"
secondCode = [-1 1 1 1   1 -1 1 1   1 1 -1 1   -1 -1 -1 1   -1 1 1 -1];

% The preamble is [1 0 1 1 0 1 1 1 0 0 0 0], and
% the antipodal form of the preamble pattern is:
preamble_bits = [-1 1 -1 -1 1 -1 -1 -1 1 1 1 1];
sync_bits = preamble_bits <0;

% Expand each preamble symbol by the 20-chip secondary code to search
% at the 1 ms tracking-epoch resolution.
preamble_ms = kron(preamble_bits, secondCode);

% Use the prompt correlator as the symbol stream, skip start of record if
%   set in initSettings to avoid tracking loop transients
bits = I_P(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

% Correlate tracking output with the preamble
tlmXcorrResult = xcorr(bits, preamble_ms);

% Find all starting points of all preamble-like patterns -----------------
clear index
xcorrLength = (length(tlmXcorrResult) +  1) /2;

% Find at what index/ms the preambles start
index = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 239.99)';

% Analyze detected preamble-like patterns ================================
for ind = 1:numel(index) % For each occurrence
    
    % Check distance to all other possible sync patterns for current start
    index2 = index - index(ind);
    
    % Require another page 10 s later and at least five 500-symbol F/NAV pages.
    if (~isempty(find(index2 == 500*20,1))) && ...
            ( length(bits(index(ind):end))> 2500*20 )
        
        %=== Read bit values for CRC-24Q check and ephemeris decoding ======
        % Search every possible preamble pattern.
        temp_bits = bits(index(ind):index(ind)+ 2500*20 -1);
        
        % Group every 10 I_P to a row for corresponding bit
        I_P_group = reshape(temp_bits,20,[])';
        
        % Wipe off the 2nd code and form data bits
        navBits = sum(I_P_group .* repmat(secondCode,2500,1),2)';
        navBits = navBits <0;
        
        %--- Correct polarity of the all data bits according to preamble bits
        if(~isequal(navBits(1:12),sync_bits))
            navBits = not(navBits);
        end   
        
        % Decoding the even page part -------------------------------------
        % Pull out implied pages from the detected preamble
        pageSymInt = navBits(13:500);
        
        % De-interleave symbols
        symMat = reshape(pageSymInt,61,8)';
        pageSym = reshape(symMat,1,[])';
        
        % Restore the inverted G2 convolutional-code branch.
        pageSym(2:2:end) = ~pageSym(2:2:end);
        % Remove convolutional encoding from implied pages
        decBits = vitdec(pageSym,trellis,tblen,'trunc','hard');
        
        % Remove tail bits
        decBits = decBits(1:238);
        
        pageBits = decBits';
        
        % Detect errors in input data using CRC ---------------------------
        [~,frmError] = step(crcDet,pageBits');
        
        % CRC-24Q check was OK. Then to decode ephemeris message by message.
        if (~frmError)
            %--- Just save for first message ------------------------------
            % firstSubFrame is the starting positions of the first message
            % in the input bit stream dataBits
            firstSubFrame = index(ind) + searchStartOffset;
            
            %--- Ephemeris decoding ---------------------------------------
            eph = ephemeris(navBits',eph);
            TOW = eph.TOW;
            break
        end % if CRC is OK ...
        
    end % Page spacing and available-data checks
end % for ind = 1:size(index)
