function correValues = corrMatlabParallelE1(settings, rawSignal, E1CodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
    startIdx, chSampSize, isDataRead)
%Correlate all active Galileo E1B/E1C channels with a shared IF buffer.
%
%   Output layout:
%       correValues = [I_E; Q_E; I_P; Q_P; I_L; Q_L;
%                      pilot_I_E; pilot_Q_E; pilot_I_P;
%                      pilot_Q_P; pilot_I_L; pilot_Q_L]
%   where each column corresponds to one tracking channel.

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

% Number of active tracking channels in this integration epoch.
channelCnt = numel(startIdx);

% Six data and six pilot correlation values per channel.
correValues = zeros(12, channelCnt);

% Early-late offset on the 12-point-per-primary-chip local-code grid.
earlyLateSpc = settings.dllCorrelatorSpacing;

for channelNr = 1:channelCnt
    % Find the size of the current code period in whole samples for the
    % current channel, and the start index of this channel block within the
    % 1-second rawSignal buffer.
    blksize = chSampSize(channelNr);
    codeStartIdx = startIdx(channelNr);

    % Each data/pilot branch contains 12 periodic guard points at each end.
    % remCodePhase includes the corresponding 11-point index offset.
    E1Code = E1CodeTable(channelNr, :);
    codeLen = length(E1Code)/2;
    E1BCodeD = E1Code(1:codeLen);
    if (settings.pilotTRKflag == 1)
        E1CCodeP = E1Code(1+codeLen:end);
    end

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
    earlyCodeD = E1BCodeD(tcode);
    if (settings.pilotTRKflag == 1)
        earlyCodeP = E1CCodeP(tcode);
    end

    % Define index into late code vector
    tcode = ceil(codePhase + earlyLateSpc) + 1;
    lateCodeD = E1BCodeD(tcode);
    if (settings.pilotTRKflag == 1)
        lateCodeP = E1CCodeP(tcode);
    end

    % Define index into prompt code vector
    tcode = ceil(codePhase) + 1;
    promptCodeD = E1BCodeD(tcode);
    if (settings.pilotTRKflag == 1)
        promptCodeP = E1CCodeP(tcode);
    end

    %% Generate the carrier frequency to mix the signal to baseband -------
    % carrPhaseStep is already the carrier phase step in radians per
    % sample, so the local carrier phase is formed directly with sample
    % indices instead of a time vector in seconds.
    trigarg = carrPhaseStep(channelNr) * sampleIndex + remCarrPhase(channelNr);
    carrCos = cos(trigarg);
    carrSin = sin(trigarg);

    %% Accumulate Early/Prompt/Late I/Q values for each code branch --------
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
    if (settings.pilotTRKflag == 1)
        % Now get early, late, and prompt values for each
        I_E_P = sum(earlyCodeP  .* iBasebandSignal);
        Q_E_P = sum(earlyCodeP  .* qBasebandSignal);
        I_P_P = sum(promptCodeP .* iBasebandSignal);
        Q_P_P = sum(promptCodeP .* qBasebandSignal);
        I_L_P = sum(lateCodeP   .* iBasebandSignal);
        Q_L_P = sum(lateCodeP   .* qBasebandSignal);
        correValues(7:12, channelNr) = [I_E_P; Q_E_P; I_P_P; Q_P_P; I_L_P; Q_L_P];
    end
end
