function [eph, firstSubFrame,TOW] = CNAV2decoding(trackResults,channelNr,settings)
%CNAV2DECODING Decode GPS L1C CNAV-2 navigation data.
%   Uses pilot overlay-code correlation for frame synchronization and
%   decodes the CNAV-2 LDPC subframes with Min-Sum and the H2/H3
%   parity-check matrices. CRC-24Q is checked before
%   ephemeris decoding.
%
%[eph, firstSubFrame, TOW] = CNAV2decoding(trackResults, channelNr, settings)
%
%   Inputs:
%       trackResults  - Tracking results with data and pilot prompt I/Q.
%       channelNr     - Index of the tracking channel to decode.
%       settings      - Receiver settings and LDPC decoder options.
%   Outputs:
%       eph           - Decoded CNAV-2 ephemeris and auxiliary parameters.
%       firstSubFrame - One-based start epoch of the first accepted frame
%                       in the 10 ms tracking stream; inf if unavailable.
%       TOW           - GPS seconds of week at that frame's start;
%                       inf if no ephemeris is accepted.

%--- Initialize ephemeris structure --------------------------------------
eph = eph_structure_init();
firstSubFrame = inf;
TOW = inf;

%% Frame sync using pilot channel overlay code ============================
% Joint tracking uses pilot I; data-only loop tracking uses pilot Q.
if isfield(settings, 'pilotTRKflag') && settings.pilotTRKflag == 1 && ...
        isfield(trackResults, 'Pilot_I_P')
    syncSoft = trackResults(channelNr).Pilot_I_P;
else
    syncSoft = trackResults(channelNr).Pilot_Q_P;
end

syncBits = syncSoft;
syncBits(syncBits > 0)  =  1;
syncBits(syncBits <= 0) = -1;

secondary = generate2ndCode(trackResults(channelNr).PRN);

xcorrResult = xcorr(syncBits, secondary);
xcorrLength = (length(xcorrResult) + 1) / 2;
xcorrResult = xcorrResult(xcorrLength : xcorrLength * 2 - 1);
index = find(abs(xcorrResult) >= 1750).';
% LDPC decoder path
thisDir = fileparts(mfilename('fullpath'));
ldpcDir = fullfile(thisDir, 'ldpcDecoder');
addpath(ldpcDir);

persistent H2 H3
if isempty(H2) || isempty(H3)
    h2File = load(fullfile(ldpcDir, 'H2.mat'), 'H2');
    h3File = load(fullfile(ldpcDir, 'H3.mat'), 'H3');
    H2 = h2File.H2;
    H3 = h3File.H3;
end

dataSoftI = trackResults(channelNr).I_P;
if isfield(trackResults, 'Q_P')
    dataSoftQ = trackResults(channelNr).Q_P;
else
    dataSoftQ = zeros(size(dataSoftI));
end

frameLength = 1800;
toiLength = 52;
decodedNavLength = 883;

for i = 1:numel(index)
    frameStart = index(i);
    frameEnd = frameStart + frameLength - 1;

    if frameStart < 1 || frameEnd > length(dataSoftI) || frameEnd > length(dataSoftQ)
        continue;
    end

    frameSoftI = dataSoftI(frameStart:frameEnd);
    frameSoftQ = dataSoftQ(frameStart:frameEnd);

    hardBits = frameSoftI;
    hardBits(hardBits > 0)  = 1;
    hardBits(hardBits <= 0) = 0;

    checkBits = 1 - 2 * hardBits(1:toiLength);
    [flag, decodedToi] = BCH51_8Decoding(checkBits);

    if flag == 0
        hardBits = 1 - hardBits;
        checkBits = 1 - 2 * hardBits(1:toiLength);
        [flag, decodedToi] = BCH51_8Decoding(checkBits);

        if flag == 0
            continue;
        end

        frameSoftI = -frameSoftI;
        frameSoftQ = -frameSoftQ;
    end

    decodedNavBits = zeros(1, decodedNavLength);
    decodedNavBits(1:9) = decodedToi;

    deinterleavedSoftI = cnav2Deinterleave(frameSoftI(53:frameLength));
    deinterleavedSoftQ = cnav2Deinterleave(frameSoftQ(53:frameLength));

    sf2SoftI = deinterleavedSoftI(1:1200).';
    sf2SoftQ = deinterleavedSoftQ(1:1200).';
    sf3SoftI = deinterleavedSoftI(1201:end).';
    sf3SoftQ = deinterleavedSoftQ(1201:end).';

    sf2Settings = cnav2LdpcSettings(settings, 'L1C_2');
    sf3Settings = cnav2LdpcSettings(settings, 'L1C_3');

    [sf2Decoded, sf2Iter, sf2HasError] = ...
        cnav2DecodeLdpc(sf2SoftI, sf2SoftQ, H2, sf2Settings);

    if sf2HasError
        [sf2DecodedNeg, sf2IterNeg, sf2HasErrorNeg] = ...
            cnav2DecodeLdpc(-sf2SoftI, -sf2SoftQ, H2, sf2Settings);

        if ~sf2HasErrorNeg
            decodedNavBits(1) = 1 - decodedNavBits(1);
            sf2Decoded = sf2DecodedNeg;
            sf2Iter = sf2IterNeg;
            sf2HasError = false;
            sf3SoftI = -sf3SoftI;
            sf3SoftQ = -sf3SoftQ;
        end
    end

    [sf3Decoded, sf3Iter, sf3HasError] = ...
        cnav2DecodeLdpc(sf3SoftI, sf3SoftQ, H3, sf3Settings);

    if ~sf2HasError && ~sf3HasError
        decodedNavBits(10:609) = sf2Decoded(1:600).';
        decodedNavBits(610:883) = sf3Decoded(1:274).';

        %--- Check CRC-24Q of Subframes 2 and 3 ---------------------------
        sf2HasError = cnav2CheckCrc(decodedNavBits(10:609));
        sf3HasError = cnav2CheckCrc(decodedNavBits(610:883));

        if sf2HasError || sf3HasError
            continue;
        end

        navBitsChar = char(decodedNavBits + '0');
        eph = ephemeris(navBitsChar, eph);

        if isinf(TOW) && ~isempty(eph.flag) && eph.flag == 1
            TOW = eph.TOW - 18;
            firstSubFrame = frameStart;
        end

        if isfield(settings, 'plotNavigation') && settings.plotNavigation
            fprintf(' PRN %02d CNAV2 LDPC decoded: SF2 iter=%d, SF3 iter=%d\n', ...
                trackResults(channelNr).PRN, sf2Iter, sf3Iter);
        end
    else
        fprintf([' PRN %02d LDPC matrix check failed: residual bit errors remain. ', ...
            'SF2_Err=%d, SF3_Err=%d\n'], ...
            trackResults(channelNr).PRN, sf2HasError, sf3HasError);
    end
end
end

function deinterleaved = cnav2Deinterleave(interleaved)
% Deinterleave the 1748 encoded symbols of CNAV-2 Subframes 2 and 3.
matDeint = reshape(interleaved, 38, 46).';
deinterleaved = matDeint(:).';
end

function [decodedCodeword, realIter, hasError] = ...
    cnav2DecodeLdpc(softI, softQ, hMatrix, ldpcSettings)
% Form LLRs, run Min-Sum decoding, and check the parity-check syndrome.
[normIp, normSigma2] = noiseEst(softI(:), softQ(:));

% In this receiver's CNAV-2 prompt stream, positive prompt-I samples map to
% logical 1 after the TOI polarity check. initLlr uses the conventional
% 0 -> positive BPSK sign, so invert the normalized samples here.
llr = initLlr(-normIp(:), normSigma2);

[decodedCodeword, realIter] = msDecode(llr, hMatrix, ldpcSettings);
decodedCodeword = decodedCodeword(:);
hasError = any(mod(sparse(logical(hMatrix)) * decodedCodeword, 2));
end

function ldpcSettings = cnav2LdpcSettings(settings, messageType)
% Select the CNAV-2 subframe dimensions and Min-Sum decoder parameters.
ldpcSettings.messageType = messageType;
ldpcSettings.maxIterations = getOptionalSetting(settings, ...
    'ldpcMaxIterations', getOptionalSetting(settings, 'maxIterations', 50));
ldpcSettings.alpha = getOptionalSetting(settings, ...
    'ldpcMsAlpha', getOptionalSetting(settings, 'alpha', 1));
ldpcSettings.offset = getOptionalSetting(settings, ...
    'ldpcMsOffset', getOptionalSetting(settings, 'offset', 0));

switch messageType
    case 'L1C_2'
        ldpcSettings.ldpcN = 1200;
        ldpcSettings.ldpcK = 600;
    case 'L1C_3'
        ldpcSettings.ldpcN = 548;
        ldpcSettings.ldpcK = 274;
end
end

function value = getOptionalSetting(settings, fieldName, defaultValue)
% Read an optional decoder setting or use its default value.
if isstruct(settings) && isfield(settings, fieldName)
    value = settings.(fieldName);
else
    value = defaultValue;
end
end
