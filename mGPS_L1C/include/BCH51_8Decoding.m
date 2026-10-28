function [flag, decodedBits]= BCH51_8Decoding(checkBits)

% This function decodes the 52-symbol first subframe of CNAV-2 using
% 2^8 BCH(51,8) hypotheses. The best hypothesis and its polarity recover
% the 9-bit Time of Interval (TOI).
%
% [flag,decodedBits] = BCH51_8Decoding(checkBits)
%
%   Inputs:
%       checkBits    - Row vector of 52 message symbols in bipolar (-1, +1).
%   Outputs:
%       flag         - Indicator for  decoding of the BCH(51,8) data 
%                      symbols, 1:success, 0: fail.
%       decodedBits  - Nine decoded TOI bits, or [] if decoding fails.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Written by Yafeng Li
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
%$Id: BCH51_8Decoding.m,v 1.1.2.5 2018/07/24 22:00:00 dpl Exp $

% A correlation threshold of 40 permits up to six errors in 52 symbols.
threshold = 40; 

% Ensure input is a row vector
checkBits = checkBits(:).';

% Check input length
if length(checkBits) ~= 52
    error('Input checkBits must be length 52 for L1C TOI decoding.');
end

% Total 2^8 = 256 hypotheses for the 8 LSBs
numHypotheses = 2^8;
hypoBits8 = zeros(numHypotheses, 8);
correValue = zeros(1, numHypotheses);

%% Generate BCH hypotheses for the eight least-significant TOI bits =======
% Reference: IS-GPS-800J Figure 3.2-4 / Text Section 3.2.3.2


for i = 0 : numHypotheses - 1
    idx = i + 1;
    
    % Candidate TOI bits in order b8,...,b1.
    current8Bits = bitget(i, 8:-1:1); 
    hypoBits8(idx, :) = current8Bits;
    
    % Load b1,...,b8 into the BCH shift register.
    reg = fliplr(current8Bits); 
    
    generatedSeq = zeros(1, 51);
    
    for k = 1:51
        % The last register stage (8th) supplies the encoded output.
        outBit = reg(8);
        generatedSeq(k) = outBit;
        
        % Feedback taps for 1+x+x^4+x^5+x^6+x^7+x^8 (763 octal).
        feedback = xor(reg(1), reg(4));
        feedback = xor(feedback, reg(5));
        feedback = xor(feedback, reg(6));
        feedback = xor(feedback, reg(7));
        feedback = xor(feedback, reg(8));
        
        % Shift toward stage 8 and insert feedback at stage 1.
        reg = [feedback, reg(1:7)];
    end
    
    % Form [b9, BCH xor b9] with b9=0; correlation polarity determines b9.
    
    hypoSequenceBin = [0, generatedSeq];
    
    % Convert to Bipolar for Correlation (+1 for 0, -1 for 1)
    hypoSyms = zeros(1, 52);
    hypoSyms(hypoSequenceBin == 0) = 1;
    hypoSyms(hypoSequenceBin == 1) = -1;
    
    % Correlate the received symbols with the b9=0 hypothesis.
    correValue(idx) = sum(hypoSyms .* checkBits);
end

%% Select the hypothesis and recover the TOI most-significant bit =========
[maxAbsVal, pos] = max(abs(correValue));

% Retrieve the raw correlation value (signed) to determine b9
rawCorr = correValue(pos);

% Check threshold
if maxAbsVal >= threshold
    flag = 1;
    
    % Extract the 8 LSBs corresponding to the best match
    best8Bits = hypoBits8(pos, :);
    
    % Positive correlation indicates b9=0; negative correlation indicates b9=1.
    if rawCorr > 0
        b9 = 0;
    else
        b9 = 1;
    end
    
    % Combine to form full 9-bit TOI message
    decodedBits = [b9, best8Bits];
else
    flag = 0;
    decodedBits = [];
end

end
 
