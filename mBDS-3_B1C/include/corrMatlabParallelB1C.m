function correValues = corrMatlabParallelB1C(settings, rawSignal, B1CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
    startIdx, chSampSize, isDataRead)
%CORRMATLABPARALLELB1C MATLAB implementation of the channel-parallel BDS-3 B1C tracking correlator.
%   This function keeps the same calling interface as the SIMD / GPU
%   correlators so it can be used directly in trkChannelsParallel.m when
%   settings.correlatorType == 0.
%
%   Each output column contains six E/P/L I/Q values per branch:
%       rows 1:6   - data BOC(1,1)
%       rows 7:12  - pilot BOC(1,1)
%       rows 13:18 - pilot BOC(6,1), zero when fullBandEn == 0.

persistent rawSignalI rawSignalQ rawSignalT

% For complex data, split the interleaved buffer [I0 Q0 I1 Q1 ...] into I
% and Q row vectors once for each newly read data block. If the current
% call reuses the same rawSignal block, keep the cached I/Q buffers.
if (settings.fileType == 1) && (isDataRead || isempty(rawSignalT))
    rawSignalT = rawSignal';
elseif (settings.fileType == 2) && (isDataRead || isempty(rawSignalI) || isempty(rawSignalQ))
    rawSignalI = rawSignal(1:2:end)';
    rawSignalQ = rawSignal(2:2:end)';
end

% Number of active tracking channels in this coherent integration.
channelCnt = numel(startIdx);

% Output order per branch: [IE QE IP QP IL QL]^T.
correValues = zeros(18, channelCnt);

% Early/late offset in 12-point-per-primary-chip code-table units.
earlyLateSpc = settings.dllCorrelatorSpacing;

for channelNr = 1:channelCnt
    % Find the size of the current code period in whole samples for the
    % current channel, and the start index of this channel block within the
    % 1-second rawSignal buffer.
    blksize = chSampSize(channelNr);
    codeStartIdx = startIdx(channelNr);

    % Get a vector with the local code for the current channel. The table
    % already contains the guard chips required by ceil(tcode) + 1.
    B1CCode = B1CCodeTable(channelNr, :);
    codeLen = length(B1CCode)/3;
    B1CCodeD = B1CCode(1:codeLen);
    B1CCodeP = B1CCode(1+codeLen:codeLen*2);
    pilotBOC61 = B1CCode(1+codeLen*2:end);


    % Extract the current signal block to be processed by this channel.
    if settings.fileType == 1
        rawSignalBlock = rawSignalT(codeStartIdx : codeStartIdx + blksize - 1);
    else
        rawSignalBlockI = rawSignalI(codeStartIdx : codeStartIdx + blksize - 1);
        rawSignalBlockQ = rawSignalQ(codeStartIdx : codeStartIdx + blksize - 1);
    end

    %% Set up all the code phase tracking information ---------------------
    sampleIndex = (0:blksize-1);
    codePhase = remCodePhase(channelNr) + codePhaseStep(channelNr) .* sampleIndex;

    % Define index into early code vector
    tcode = ceil(codePhase - earlyLateSpc) + 1;
    earlyCodeD = B1CCodeD(tcode);
    earlyCodeP = B1CCodeP(tcode);
    p61_earlyCode = pilotBOC61(tcode);


    % Define index into late code vector
    tcode = ceil(codePhase + earlyLateSpc) + 1;
    lateCodeD = B1CCodeD(tcode);
    lateCodeP = B1CCodeP(tcode);
    p61_lateCode = pilotBOC61(tcode);

    % Define index into prompt code vector
    tcode = ceil(codePhase) + 1;
    promptCodeD = B1CCodeD(tcode);
    promptCodeP = B1CCodeP(tcode);
    p61_promptCode  = pilotBOC61(tcode);

    %% Generate the carrier frequency to mix the signal to baseband -------
    % carrPhaseStep is already the carrier phase step in radians per
    % sample, so the local carrier phase is formed directly with sample
    % indices instead of a time vector in seconds.
    trigarg = carrPhaseStep(channelNr) * sampleIndex + remCarrPhase(channelNr);
    carrCos = cos(trigarg);
    carrSin = sin(trigarg);

    %% Accumulate E/P/L I/Q correlations for each signal branch ------------
    % First mix to baseband
    if settings.fileType == 1
        iBasebandSignal = rawSignalBlock .* carrCos;
        qBasebandSignal = -rawSignalBlock .* carrSin;
    else
        iBasebandSignal = rawSignalBlockI .* carrCos + rawSignalBlockQ .* carrSin;
        qBasebandSignal = rawSignalBlockQ .* carrCos - rawSignalBlockI .* carrSin;
    end

    % Now get early, late, and prompt values for each
    I_E = sum(earlyCodeD  .* iBasebandSignal);
    Q_E = sum(earlyCodeD  .* qBasebandSignal);
    I_P = sum(promptCodeD .* iBasebandSignal);
    Q_P = sum(promptCodeD .* qBasebandSignal);
    I_L = sum(lateCodeD   .* iBasebandSignal);
    Q_L = sum(lateCodeD   .* qBasebandSignal);

    % Output order matches the SIMD / GPU correlator interface used by
    % trkChannelsParallel.m.
    correValues(1:6, channelNr) = [I_E; Q_E; I_P; Q_P; I_L; Q_L];

    % Now get early, late, and prompt values for each
    I_E_P = sum(earlyCodeP  .* iBasebandSignal);
    Q_E_P = sum(earlyCodeP  .* qBasebandSignal);
    I_P_P = sum(promptCodeP .* iBasebandSignal);
    Q_P_P = sum(promptCodeP .* qBasebandSignal);
    I_L_P = sum(lateCodeP   .* iBasebandSignal);
    Q_L_P = sum(lateCodeP   .* qBasebandSignal);
    correValues(7:12, channelNr) = [I_E_P; Q_E_P; I_P_P; Q_P_P; I_L_P; Q_L_P];

    if (settings.fullBandEn == 1)
        % Correlation values for pilot BOC(6,1) spreading waveform
        p61_I_E = sum(p61_earlyCode  .* iBasebandSignal);
        p61_Q_E = sum(p61_earlyCode  .* qBasebandSignal);
        p61_I_P = sum(p61_promptCode .* iBasebandSignal);
        p61_Q_P = sum(p61_promptCode .* qBasebandSignal);
        p61_I_L = sum(p61_lateCode   .* iBasebandSignal);
        p61_Q_L = sum(p61_lateCode   .* qBasebandSignal);
        correValues(13:18, channelNr) = [p61_I_E, p61_Q_E, p61_I_P, p61_Q_P, p61_I_L, p61_Q_L];
    end

end
