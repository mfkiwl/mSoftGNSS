function [eph, firstSubFrame,TOW] = BCNAV3decoding(I_P_InputBits, Q_P_InputBits,PRNin)
%BCNAV3DECODING Find BDS B2b frames and decode B-CNAV3 messages.
%Frame polarity is resolved from the preamble. The 972 coded symbols are
%LDPC-decoded, checked with CRC-24Q, and passed to the ephemeris decoder.
%
%[eph, firstSubFrame,TOW] = BCNAV3decoding(I_P_InputBits,Q_P_InputBits,PRNin)
%
%   Inputs:
%       I_P_InputBits   - prompt I output, one value per 1 ms code period.
%       Q_P_InputBits   - prompt Q output used for noise estimation.
%       PRNin           - satellite PRN used to check the frame header.
%   Outputs:
%       firstSubFrame   - index of the first CRC-valid frame in the 1 ms
%                       tracking stream; inf if no valid frame is found.
%       TOW             - BDS second of week of that frame, in seconds.
%       eph             - decoded satellite ephemeris.
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
%--------------------------------------------------------------------------
%--- Initialize ephemeris structute  --------------------------------------
eph = eph_structure_init_B2b();
% Starting positions of the first message in the input stream trackResults.I_P
firstSubFrame = inf;
% TOW of the first message
TOW = inf;
%% Bit and frame synchronization ====================================
searchStartOffset = 0;
% [1 1 1 0 1 0 1 1 1 0 0 1 0 0 0 0], and the antipodal 
% form of the preamble pattern is:%B2b  EB90
preamble_bits = [-1 -1 -1 1 -1 1 -1 -1 -1 1 1 -1 1 1 1 1];
preamble_ms = preamble_bits;

% ==================== [Soft and Hard Information Extraction Initialization] ====================
% 1. Extract raw soft bits from the I-branch (containing signal and noise)
raw_soft_bits_I = I_P_InputBits(1 + searchStartOffset : end);
% 2. Extract raw soft bits from the Q-branch (containing complex noise or degraded real noise)
raw_soft_bits_Q = Q_P_InputBits(1 + searchStartOffset : end);
% 3. Generate hard decision bits specifically for preamble cross-correlation
bits_hard = raw_soft_bits_I;
bits_hard(bits_hard > 0)  =  1;
bits_hard(bits_hard <= 0) = -1;

% Correlate tracking output with the preamble (Find synchronization header using hard decisions)
tlmXcorrResult = xcorr(bits_hard, preamble_ms);
% Find all starting points off all preamble like patterns -----------------
clear index
xcorrLength = (length(tlmXcorrResult) +  1) /2;
index = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 15.9)';

% Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Path Configuration for LDPC decoder
addpath ('include/ldpcDecoder');

%% B-CNAV3 decoding =================================================
for i = 1:length(index) % For each occurrence
    % Ensure the search for i has the number of a whole message(1000)
    if ((length(bits_hard) - index(i) + 1) >= 1000 )    
         
        % Synchronized Extraction of Soft and Hard Sequences ==============
        % Extract continuous level data for the I and Q branches of the 
        % current frame (1000 bits)
        navBits_soft_I = raw_soft_bits_I(index(i) : index(i)+1000-1);
        navBits_soft_Q = raw_soft_bits_Q(index(i) : index(i)+1000-1);
        
        % Extract hard bits for polarity determination and PRN parsing
        navBits_hard = bits_hard(index(i) : index(i)+1000-1);
        
        % Polarity reversal logic: Due to BPSK phase ambiguity, if the 
        % preamble is inverted, the signal in the complex plane is entirely
        %  rotated by 180 degrees
        if(~isequal(navBits_hard(1:16), preamble_bits))
            navBits_hard   = -navBits_hard;
            navBits_soft_I = -navBits_soft_I; 
            navBits_soft_Q = -navBits_soft_Q;
        end
        
        % Parse PRN (based on hard bits)
        temp_PRN = dec2bin(navBits_hard(17:22) < 0.5);
        PRN = bin2dec(temp_PRN');
        if (PRN < 6) || (PRN > 58 || (PRN ~= PRNin))
            continue; 
        end
        eph.PRN = PRN;
        
        % ------------- B2b LDPC decoder ----------------------------------
        decodedNavBits = B2bLDPCDecoder(navBits_soft_I, navBits_soft_Q);
      
        %--- To do CRC-24Q check for current frame ------------------------
        NavDataLogic = (decodedNavBits == 1);
        [~,frmError] = step(crcDet, NavDataLogic);
        
        % Decode ephemeris message by message  ----------------------------
        if (~frmError)
            navBitsChar = char(NavDataLogic + '0');
            eph = ephemeris_B2b(navBitsChar', eph);
            if isinf(firstSubFrame)
                firstSubFrame = index(i) + searchStartOffset;
                TOW = eph.SOW;
            end
        end % if CRC is OK ...
    end
end   
end
