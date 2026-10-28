function correValues = corrMatlabParallelL1C(settings, rawSignal, L1CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
    startIdx, chSampSize, isDataRead)
%Correlate all active GPS L1C data/pilot channels using MATLAB.
%
%correValues = corrMatlabParallelL1C(settings, rawSignal, L1CCodeTable, ...
%    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
%    startIdx, chSampSize, isDataRead)
%
%   Inputs:
%       settings      - Receiver settings; dllCorrelatorSpacing is in
%                       1/12-primary-chip code-table points.
%       rawSignal     - Shared buffer of real samples or interleaved I/Q.
%       L1CCodeTable  - One [data, pilot] row per channel; each branch has
%                       12 points per primary chip and 12 guards at each end.
%       remCarrPhase  - Initial carrier phases [rad].
%       carrPhaseStep - Carrier phase increments [rad/IF sample].
%       remCodePhase  - Initial phases on the 12-point grid, including the
%                       caller's +11 guard-point offset.
%       codePhaseStep - Code phase increments [1/12 primary chip/IF sample].
%       startIdx      - One-based start index of each channel in the buffer.
%       chSampSize    - Number of logical IF samples for each channel.
%       isDataRead    - 1 for a newly read buffer; 0 to reuse cached samples.
%   Outputs:
%       correValues   - 12 x channelCnt matrix; rows 1:6 are data
%                       [I_E Q_E I_P Q_P I_L Q_L], rows 7:12 are pilot.
%                       Pilot prompt is also used for CNAV-2 frame sync.
%   startIdx counts logical samples: one real value or one complex I/Q pair.

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

% Output order per branch: [IE QE IP QP IL QL]^T.
correValues = zeros(12, channelCnt);

% Define early-late offset on the 12-point-per-primary-chip grid.
earlyLateSpc = settings.dllCorrelatorSpacing;

for channelNr = 1:channelCnt
    % Find the size of the current code period in whole samples for the
    % current channel, and the start index of this channel block within the
    % 1-second rawSignal buffer.
    blksize = chSampSize(channelNr);
    codeStartIdx = startIdx(channelNr);

    % Get a vector with the local code for the current channel. The table
    % contains 12 guard points at each end of each data/pilot branch.
    L1CCode = L1CCodeTable(channelNr, :);
    % Length of one guarded branch, not the combined data/pilot table.
    codeLen = length(L1CCode)/2;
    L1CCodeD = L1CCode(1:codeLen);
    L1CCodeP = L1CCode(1+codeLen:end);

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
    earlyCodeD = L1CCodeD(tcode);
    earlyCodeP = L1CCodeP(tcode);

    % Define index into late code vector
    tcode = ceil(codePhase + earlyLateSpc) + 1;
    lateCodeD = L1CCodeD(tcode);
    lateCodeP = L1CCodeP(tcode);

    % Define index into prompt code vector
    tcode = ceil(codePhase) + 1;
    promptCodeD = L1CCodeD(tcode);
    promptCodeP = L1CCodeP(tcode);

    %% Generate the carrier frequency to mix the signal to baseband -------
    % Local carrier phase in radians.
    trigarg = carrPhaseStep(channelNr) * sampleIndex + remCarrPhase(channelNr);
    carrCos = cos(trigarg);
    carrSin = sin(trigarg);

    %% Generate six accumulated I/Q values per branch ---------------------
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

    % Data correlations in E/P/L order.
    correValues(1:6, channelNr) = [I_E; Q_E; I_P; Q_P; I_L; Q_L];
    % Pilot correlations for tracking and CNAV-2 frame synchronization.
    I_E_P = sum(earlyCodeP  .* iBasebandSignal);
    Q_E_P = sum(earlyCodeP  .* qBasebandSignal);
    I_P_P = sum(promptCodeP .* iBasebandSignal);
    Q_P_P = sum(promptCodeP .* qBasebandSignal);
    I_L_P = sum(lateCodeP   .* iBasebandSignal);
    Q_L_P = sum(lateCodeP   .* qBasebandSignal);
    correValues(7:12, channelNr) = [I_E_P; Q_E_P; I_P_P; Q_P_P; I_L_P; Q_L_P];
end
