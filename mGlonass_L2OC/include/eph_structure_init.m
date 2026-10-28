function eph = eph_structure_init()
% Initialize GLONASS L2OC ephemeris/status structure.
%
% Note:
%   The ICD "GLONASS CDMA L2, Edition 1.0, 2016" defines
%   L2OCp signal structure, PRN generator, OC2 and L2CSI/L2OCp TDM.
%   It does not define the L2 CSI navigation message format or ephemeris
%   bit fields. This structure records tracking status without decoded
%   ephemeris data.

%--------------------------------------------------------------------------
% GLONASS L2OC signal metadata
%--------------------------------------------------------------------------
eph.system          = 'GLONASS';
eph.signal          = 'L2OCp';
eph.component       = 'pilot';
eph.navComponent    = 'L2CSI';
eph.navDecoded      = 0;

% Nominal L2 CDMA carrier frequency from ICD: 1248.06 MHz
eph.carrierFreq     = 1248.06e6;

% L2OCp PRN parameters from ICD
eph.L2OCpCodeLength = 10230;      % PRN_L2OCp chips
eph.L2OCpCodeRate   = 0.5115e6;   % chips/s
eph.L2OCpPeriod     = 0.020;      % 20 ms

% Insert zero-valued L2CSI slots to form a 20460-point TDM code table.
eph.equivCodeLength = 20460;
eph.equivCodeRate   = 1.023e6;
eph.equivCodePeriod = 0.020;

%--------------------------------------------------------------------------
% Undecoded navigation/ephemeris fields
%--------------------------------------------------------------------------
% No L2 CSI navigation message is decoded; ephemeris fields remain empty.
eph.SVID        = [];
eph.PRN         = [];
eph.slot        = [];
eph.freqChannel = [];

% Time-related fields
eph.tow         = [];
eph.tk          = [];
eph.tb          = [];
eph.N4          = [];
eph.NT          = [];

% Clock correction fields
eph.tau_n       = [];
eph.gamma_n     = [];
eph.tau_c       = [];
eph.tau_gps     = [];

% Orbit state-vector fields
eph.x           = [];
eph.y           = [];
eph.z           = [];
eph.xDot        = [];
eph.yDot        = [];
eph.zDot        = [];
eph.xDotDot     = [];
eph.yDotDot     = [];
eph.zDotDot     = [];

% Health/status fields
eph.health      = [];
eph.P           = [];
eph.P1          = [];
eph.P2          = [];
eph.P3          = [];
eph.P4          = [];

% Raw data/debug fields
eph.rawNavBits      = [];
eph.rawNavBitsLength = 0;
eph.decodeStatus    = 'L2OCp pilot tracking only; L2 CSI navigation message is not decoded.';

% Decoded flag
% 0: ephemeris not available
% 1: ephemeris decoded
eph.flag        = 0;

end
