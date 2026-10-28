function [eph, firstSubFrame, TOW] = BCNAV1decoding(trackResults, channelNr, settings)
% BCNAV1decoding decodes BeiDou B1C B-CNAV1 navigation message.
%
% This function follows the B1C BCH/LDPC decoding structure:
%
%   decodedNavBits = B1CLDPCDecoder(navBits_soft_I, navBits_soft_Q);
%
% The main function is responsible for:
%   1) frame synchronization using pilot secondary code;
%   2) extracting one complete 1800-symbol B-CNAV1 frame;
%   3) calling B1CLDPCDecoder;
%   4) CRC-24Q checking for SF2 and SF3;
%   5) ephemeris decoding.

% -------------------------------------------------------------------------
% Initialize ephemeris structure
% -------------------------------------------------------------------------
eph = eph_structure_init();

firstSubFrame = inf;
TOW = inf;

%% Frame synchronization using pilot secondary code ========================
if settings.pilotACQflag == 1
    syncBits_soft = trackResults(channelNr).Pilot_I_P;
else
    % Use the Q prompt for frame synchronization when pilotACQflag is zero.
    syncBits_soft = trackResults(channelNr).Pilot_Q_P;
end

% Hard decision only for synchronization
syncBits_hard = syncBits_soft;
syncBits_hard(syncBits_hard > 0)  =  1;
syncBits_hard(syncBits_hard <= 0) = -1;

% Generate B1C pilot secondary code
PRNin = trackResults(channelNr).PRN;
Secondary = generate2ndCode(PRNin);

% Correlate tracking output with pilot secondary code
XcorrResult = xcorr(syncBits_hard, Secondary);

xcorrLength = (length(XcorrResult) + 1) / 2;
XcorrResult = XcorrResult(xcorrLength : xcorrLength * 2 - 1);

% Each B-CNAV1 frame has 1800 symbols
index = find(abs(XcorrResult) >= 1799.5)';

% CRC-24Q detector
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Decoder path
addpath('include/ldpcDecoder');

%% B1C data decoding and ephemeris extraction ==============================
raw_soft_bits_I = trackResults(channelNr).I_P;

% Q branch is used as noise/reference information for LLR generation.
% Use the Q prompt for noise estimation; substitute zeros if it is absent.
if isfield(trackResults, 'Q_P')
    raw_soft_bits_Q = trackResults(channelNr).Q_P;
else
    raw_soft_bits_Q = zeros(size(raw_soft_bits_I));
end

% Make sure I/Q branches have the same length
minLen = min(length(raw_soft_bits_I), length(raw_soft_bits_Q));
raw_soft_bits_I = raw_soft_bits_I(1:minLen);
raw_soft_bits_Q = raw_soft_bits_Q(1:minLen);

for i = 1:length(index)

    % Ensure one complete B-CNAV1 frame is available
    if (length(raw_soft_bits_I) - index(i) + 1) < 1800
        continue;
    end

    % ---------------------------------------------------------------------
    % Extract one complete B-CNAV1 frame
    % ---------------------------------------------------------------------
    navBits_soft_I = raw_soft_bits_I(index(i) : index(i) + 1800 - 1);
    navBits_soft_Q = raw_soft_bits_Q(index(i) : index(i) + 1800 - 1);

    % ---------------------------------------------------------------------
    % B1C BCH + LDPC / EMS decoder
    % ---------------------------------------------------------------------
    decodedNavBits = B1CLDPCDecoder(navBits_soft_I, navBits_soft_Q);

    % If BCH decoding fails, decodedNavBits is empty
    if isempty(decodedNavBits)
        continue;
    end

    % ---------------------------------------------------------------------
    % CRC-24Q check for SF2 and SF3
    % ---------------------------------------------------------------------
    decodedNav_SF2 = decodedNavBits(15:614);
    decodedNav_SF3 = decodedNavBits(615:878);

    checkBits2 = (decodedNav_SF2 == 1);
    [~, frmError1] = step(crcDet, checkBits2(:));

    checkBits3 = (decodedNav_SF3 == 1);
    [~, frmError2] = step(crcDet, checkBits3(:));

    % ---------------------------------------------------------------------
    % Ephemeris decoding
    % ---------------------------------------------------------------------
    if (~frmError1) && (~frmError2)

        navBitsChar = char(decodedNavBits(:).' + '0');

        eph = ephemeris(navBitsChar, eph);

        if isinf(TOW)
            TOW = eph.TOW;
            firstSubFrame = index(i);
        end
    end
end

end
