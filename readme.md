An Open-Source MATLAB GNSS SDR Toolbox with SIMD/GPU Acceleration: mSoftGNSS 
===============================================================================



Overview
-------------------------------------------------------------------------------
mSoftGNSS is an open-source MATLAB-based toolbox for post-processing recorded 
GNSS intermediate-frequency (IF) signals. Building on the SoftGNSS receiver 
architecture, it supports GPS, Galileo, GLONASS, and BeiDou signals within a 
configurable processing framework. Signal-dependent capabilities include 
CPU/GPU acquisition, data/pilot tracking, navigation-message decoding, 
pseudorange generation, and positioning. The toolbox provides channel-serial 
and channel-parallel tracking modes, each supporting MATLAB-based, 
SIMD-accelerated, and GPU-accelerated correlators. SIMD implementations include
x86 AVX2 and ARM64 NEON, selected automatically for the running platform.
Channel-parallel tracking shares IF data buffers across channels and batches 
correlation tasks to reduce repeated data access. Computationally intensive 
correlators are implemented in C++ and CUDA C++ and integrated through MEX 
interfaces, while receiver scheduling, tracking-loop control, navigation 
decoding, and positioning remain in MATLAB. This separation combines the 
accessibility of high-level algorithm development with efficient native 
computation. The toolbox also supports modernized GNSS signals and LDPC 
decoding for algorithm research, receiver prototyping, and evaluation using 
recorded IF data.



Authors
-------------------------------------------------------------------------------
* Yafeng Li
    * E-Mail: <lyf8118@126.com>
    * Wechat: lyf8118521
    * QQ Group for Technical Discussions on GNSS Software Receivers: 147304049


* Dennis Akos  
    * E-Mail: <dma@colorado.edu>
    * HP: <http://www.colorado.edu/aerospace/dennis-akos>



Features
-------------------------------------------------------------------------------
* GNSS signal processing functions written in MATLAB
    * Local code generation
    * CPU and GPU acquisition
    * Channel-serial and channel-parallel tracking
    * MATLAB, SIMD MEX, and GPU/CUDA MEX correlators
    * Data/pilot tracking for modernized signals
    * Carrier-to-noise-density ratio (C/N0) estimation and tracking plots
    * Navigation-message decoding (including LDPC decoder)
    * Pseudorange generation
    * Receiver position/clock-bias calculation and result plotting
* SIMD/GPU acceleration and correlator options
    * settings.gpuACQflag = 0/1 selects CPU/GPU acquisition.
    * settings.trkMode = 0/1 selects channel-serial/channel-parallel tracking.
    * settings.correlatorType = 0 uses MATLAB-based correlators.
    * settings.correlatorType = 1 uses CPU SIMD MEX correlators.
    * settings.correlatorType = 2 uses CUDA GPU MEX correlators.
    * MATLAB correlator variants include BPSK, QPSK, TMBPSK, and QMBOC
      implementations inside each receiver's include folder.
    * Shared MEX correlator backends include BPSK, QPSK, and QMBOC
      SIMD and CUDA variants.
* Supported signals
    * GPS L1 C/A (data)
    * GPS L1C (data + pilot)
    * GPS L2C (data + pilot)
    * GPS L5 (data + pilot)
    * Galileo E1 (data + pilot)
    * Galileo E5a (data + pilot)
    * Galileo E5b (data + pilot)
    * Galileo E6 (data + pilot; HAS page decoding)
    * GLONASS L1OF (data)
    * GLONASS L2OF (data)
    * GLONASS L1OC (data + pilot)
    * GLONASS L2OCp (pilot)
    * GLONASS L3OC (data + pilot)
    * BeiDou B1I (data)
    * BeiDou B3I (data)
    * BDS-3 B1C (data + pilot)
    * BDS-3 B2a (data + pilot)
    * BeiDou B2b (data)
* IF binary-file post-processing
    * settings.dataType selects int8 or int16. Both types support real and
      complex I/Q samples.
    * settings.fileType = 1 selects real samples; 2 selects interleaved
      I/Q samples (I0,Q0,I1,Q1,...). dataType applies to each I/Q component.
    * LimeSDR complex-sample and NUT4NT real-sample examples are provided below.
      Set the data format, IF, and sampling rate to match the recording.
    * The receivers have been tested on Windows with MATLAB R2025a/R2025b.
	


Directory and Files
-------------------------------------------------------------------------------
Root folders

    ./Doc                     Documentation, ICD material, papers, and receiver
                              summary documents.
    ./IF_Data_Set             Optional location for IF data files and metadata.
    ./native_Correlators      Shared SIMD and CUDA MEX correlator source files
                              and locally compiled MEX binaries.

Receiver folders

    Each receiver folder is self-contained except for the shared ./native_Correlators,
    ./IF_Data_Set, and ./Doc folders. The current receiver folders are:

    ./mGPS_L1CA               GPS L1 C/A SDR receiver
    ./mGPS_L1C                GPS L1C SDR receiver
    ./mGPS_L2C                GPS L2C SDR receiver
    ./mGPS_L5                 GPS L5 SDR receiver
    ./mGalileo_E1             Galileo E1 SDR receiver
    ./mGalileo_E5a            Galileo E5a SDR receiver
    ./mGalileo_E5b            Galileo E5b SDR receiver
    ./mGalileo_E6             Galileo E6B/E6C SDR receiver (HAS reception)
    ./mGlonass_L1_L2          GLONASS L1OF/L2OF SDR receiver
    ./mGlonass_L1OC           GLONASS L1OC SDR receiver
    ./mGlonass_L2OC           GLONASS L2OCp SDR receiver
    ./mGlonass_L3OC           GLONASS L3OC SDR receiver
    ./mBDS_B1I                BeiDou B1I SDR receiver
    ./mBDS_B3I                BeiDou B3I SDR receiver
    ./mBDS-3_B1C              BDS-3 B1C SDR receiver
    ./mBDS-3_B2a              BDS-3 B2a SDR receiver
    ./mBDS_B2b                BeiDou B2b SDR receiver

Common layout inside each receiver folder

    ./init.m                  Receiver startup script. It sets paths, loads
                              settings, displays IF waveforms and spectra, and
                              starts processing after user confirmation.
    ./initSettings.m          Receiver-specific parameter configuration.
    ./postProcessing.m        Top-level processing script for acquisition,
                              channel initialization, tracking, navigation,
                              result saving, and plotting.
    ./acquisitionCPU.m        MATLAB/CPU acquisition implementation.
    ./acquisitionGPU.m        MATLAB GPU-array acquisition implementation.
    ./trkChannelsSerial.m     Channel-serial tracking implementation.
    ./trkChannelsParallel.m   Channel-parallel tracking implementation.
    ./postNavigation.m        Navigation-message decoding and position computation.
                              For ./mGlonass_L3OC this file is located under
                              ./include/postNavigation.m.
    ./include                 Signal-specific helper functions, local-code
                              generators, MATLAB correlators, navigation-message
                              decoders, plotting helpers, and channel setup.
                              Keep E1b.dat/E1c.dat in mGalileo_E1/include and
                              the HAS Reed-Solomon matrix in mGalileo_E6/include.
    ./include/ldpcDecoder     LDPC decoder support files, present only in
                              receivers that need LDPC-based navigation decoding
                              such as mBDS-3_B1C, mBDS-3_B2a, mBDS_B2b, and mGPS_L1C.
    ./Common                  Receiver-local common math and positioning utilities such as
                              pseudorange calculation, coordinate transforms,
                              least-squares positioning, and loop coefficients.

Shared correlator backends in ./native_Correlators

    corrX86SIMDSerial*.cpp       SIMD serial tracking correlator source files.
    corrX86SIMDParallel*.cpp     SIMD channel-parallel tracking correlator source files.
    corrNeonSIMD*.cpp            ARM64 NEON serial/channel-parallel correlator sources.
    corrGPUSerial*.cu            CUDA serial tracking correlator source files.
    corrGPUParallel*.cu          CUDA channel-parallel tracking correlator source files.
    corrGPUParallelFused*.cu     CUDA fused channel-parallel tracking correlator
                                 source files where available.
    *.mexw64 / *.mexa64          Windows x64 / Linux x86-64 MEX binaries.
    *.mexmaca64                  Native Apple silicon MEX binaries (NEON).

Correlator selection

    settings.correlatorType = 0
                              MATLAB-based correlators in the receiver's
                              ./include folder. These are easiest to inspect
                              and debug.
    settings.correlatorType = 1
                              CPU SIMD MEX correlators from ./native_Correlators.
                              The code automatically selects x86 or ARM64 NEON
                              according to the running MATLAB architecture.
                              These accelerate tracking without requiring CUDA.
    settings.correlatorType = 2
                              CUDA GPU MEX correlators from ./native_Correlators.
                              These are intended for high-throughput serial or
                              channel-parallel tracking.



Software Dependencies
-------------------------------------------------------------------------------
* MATLAB is required. Compile MEX binaries for the operating system and
  architecture of the running MATLAB.
* MATLAB Signal Processing Toolbox is used by IF spectrum plotting (pwelch).
* MATLAB Communications Toolbox is needed by signal-dependent navigation decoding:
  -- Create Galois field array: gf()
  -- BCH decoder: bchdec()
  -- Detect errors in input data using CRC: comm.CRCDetector()
  -- Convert convolutional code polynomials to trellis description: poly2trellis()
  -- Convolutionally decode binary data using Viterbi algorithm: vitdec()
  -- Detect errors in input data using CRC: step()
* MATLAB Parallel Computing Toolbox and a CUDA-capable GPU are needed for
  GPU-array acquisition. Direct CUDA MEX execution requires a compatible
  NVIDIA GPU, driver, and CUDA runtime.
* x86 SIMD MEX correlators require AVX2 and FMA; ARM64 uses the separate NEON
  sources. A supported C++ compiler is needed to build them with mex.
  A supported CUDA toolkit/compiler setup is needed to rebuild CUDA MEX files
  with mexcuda (Parallel Computing Toolbox).
* CPU-only operation is possible by selecting MATLAB or SIMD correlators and
  disabling GPU acquisition in each receiver's initSettings.m.
  
  
 
Compile correlators
-------------------------------------------------------------------------------
This section describes how to compile MATLAB MEX correlators for mSoftGNSS.
Choose the instructions for your operating system and CPU architecture. Run
the commands in the MATLAB Command Window, starting in the mSoftGNSS root.
Enter this project's native_Correlators folder to compile the SIMD and CUDA
sources.

The cudaArch value '89' below is an example for the RTX 4080 Laptop GPU.
For another GPU, use its compute capability and
ensure that the selected CUDA Toolkit supports the target. Choose the section
for your platform.

| Platform | SIMD requirements | Output files |
| --- | --- | --- |
| Windows x86-64 | x86-64 CPU with AVX2 and FMA support | .mexw64 |
| Linux x86-64 | x86-64 CPU with AVX2 and FMA support | .mexa64 |

### Windows x86-64

The SIMD example uses MinGW-w64/G++. Select a MATLAB-supported MinGW64
Compiler (C++) with mex -setup C++ and compile the six x86 SIMD correlators.

```matlab
cd('native_Correlators');
mex -setup C++
sources = dir('corrX86SIMD*.cpp');
for k = 1:numel(sources)
    mex('-R2017b', '-O', ...
        'CXXFLAGS=$CXXFLAGS -std=c++14 -mavx2 -mfma', ...
        sources(k).name);
end
```

Then switch to Microsoft Visual C++ supported by both MATLAB and the CUDA
Toolkit to compile the CUDA correlators.

```matlab
mex -setup C++
cudaArch = '89';
cudaFlags = sprintf('NVCCFLAGS=-gencode=arch=compute_%s,code=sm_%s', ...
    cudaArch, cudaArch);
sources = dir('corrGPU*.cu');
for k = 1:numel(sources)
    mexcuda('-R2017b', '-O', cudaFlags, sources(k).name);
end
```

After compiling both the SIMD and CUDA correlators, the .mexw64 files are
generated in this directory.

### Linux x86-64

Use a MATLAB-supported GCC/G++. Add it to PATH before starting MATLAB, then
configure it with mex -setup C++. First compile the x86 SIMD correlators:

```matlab
cd('native_Correlators');
mex -setup C++
sources = dir('corrX86SIMD*.cpp');
for k = 1:numel(sources)
    mex('-R2017b', '-O', ...
        'CXXFLAGS=$CXXFLAGS -std=c++14 -mavx2 -mfma', sources(k).name);
end
```

Then compile the CUDA correlators in the same directory. The host GCC/G++
compiler must be supported by both MATLAB and the CUDA Toolkit.

```matlab
cudaArch = '89';
cudaFlags = sprintf('NVCCFLAGS=-gencode=arch=compute_%s,code=sm_%s', ...
    cudaArch, cudaArch);
sources = dir('corrGPU*.cu');
for k = 1:numel(sources)
    mexcuda('-R2017b', '-O', cudaFlags, sources(k).name);
end
```

After compiling both the SIMD and CUDA correlators, the .mexa64 files are
generated in this directory.

Keep the source basenames for all MEX files and place them directly in this
project's native_Correlators folder. Check the outputs in the MATLAB Command
Window:

```matlab
mexext
dir(['*.' mexext])
```

How to use
-------------------------------------------------------------------------------
* Step 1: Enter the receiver folder for the target signal, for example
          "mBDS-3_B1C" or "mGPS_L1CA".
* Step 2: Put the IF data file under "IF_Data_Set", or set an absolute file
          path in settings.fileName.
* Step 3: Configure processing, IF-data, acquisition, tracking, correlator, and
          navigation parameters in "initSettings.m". skipNumberOfSamples counts
          logical samples: one real value or one I/Q pair counts as one sample.
* Step 4: Select CPU/GPU acquisition with settings.gpuACQflag, channel-serial or
          channel-parallel tracking with settings.trkMode, and script-based/SIMD/GPU
          correlation with settings.correlatorType.
* Step 5: Run "init.m" from the receiver folder. Inspect the IF-data plots,
          then enter 1 to start processing or 0 to exit.
* Step 6: After processing, review the acquisition, tracking, navigation, and
          positioning results and their plots.



Implementation details
-------------------------------------------------------------------------------
The receiver set follows the original SoftGNSS-style MATLAB SDR organization
and extends it with modern GNSS signals, CPU/GPU acquisition, channel-serial and
channel-parallel tracking, and shared SIMD/CUDA MEX correlator backends.

Channel-parallel tracking uses an outer loop over tracking epochs, shared IF
buffering, and batched correlator calls to process all active channels before
updating their individual tracking loops. This architecture preserves independent
channel states and loop feedback, reduces redundant IF reads and interface-call
overhead, and supports multichannel batch processing with SIMD/GPU acceleration.

See the reference for the BDS-3 B1C/B2a receiver design:
Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and implementation of an open-source 
BDS-3 B1C/B2a SDR receiver. GPS Solut 23, 60 (2019). 
https://doi.org/10.1007/s10291-019-0853-z


	   
Test signal (collected by LimeSDR) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C.bin
       -- for GPS L1 C/A, Galileo E1 and BDS B1C 
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1NdZeVrEdEjDnJiV91Cg0qg?pwd=bvh5 password: bvh5
    * GPS L1C.iq
       -- for GPS L1C 
       -- dataType: int8
       -- fileType: 2 (complex I/Q samples; each component is int8)
       -- IF: 0 Hz
       -- samplingFreq: 25e6 Hz
       -- link: https://pan.baidu.com/s/1hZ4S8R4VWB5OtPr9VyPCYQ?pwd=tv8b password: tv8b
    * GPS_L2.bin
       -- for GPS L2C
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1KqdyWyuVPAgx0EoUdoLwag?pwd=5tya password: 5tya
    * L5_E5a_B2a.bin
       -- for GPS L5, Galileo E5a and BDS B2a
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1rX2lguJgMmgjAdZIJK1FNg?pwd=bnzm password: bnzm
    
    * E5b_B2b.bin   
       -- for Galileo E5b, BDS B2b
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/128dwHIKUG9iM0MigS_7vJw?pwd=x4uq password: x4uq
    * GAL_E6.bin  
       -- for Galileo E6
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1OrHtpNG-dQxSmmVc9rX6tA?pwd=qfx8 password: qfx8
    * GLO_L1OC_L1.bin 
       -- for GLONASS L1OF  
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 1.005e6 Hz
       -- samplingFreq: 30e6 Hz
       -- GLONASS L1OC
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1w-fYtmM3oW_pxR4ZH_HGyw?pwd=js93 password: js93
    * GLO_L2OC_L2.bin  
       -- for GLONASS L2OF 
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: -2.06e6 Hz
       -- samplingFreq: 30e6 Hz
    
       -- for GLONASS L2OC 
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1jefXiCkeI9qHSqwZ_RXr2g?pwd=33w6 password: 33w6
    * GLO_L3OC.bin  
       -- for GLONASS L3OC	
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1ySNHn3W66Q6n9fxERR7nOQ?pwd=289m password: 289m
    * BDS_B1I.bin
       -- for BDS B1I 
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1qV-TCtROGHL8iJZlqw5mZA?pwd=39vd password: 39vd
    * BDS_B3I.bin       
       -- for BDS B3I    
       -- dataType: int16
       -- fileType: 2 (complex I/Q samples; each component is int16)
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1T_g1to1e8OW5OFtHlQsZeg?pwd=6w1u password: 6w1u 



Test signal (collected by NUT4NT) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C_real.bin
       -- for GPS L1 C/A, Galileo E1 and BDS B1C
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -14.58e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/1cFGCBwcWA0vHBeiL7WU_2Q?pwd=r8au password: r8au
    * B1I_real.bin
       -- for BDS B1I
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -28.902e6 Hz
       -- samplingFreq: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1y4LfFZSO-Qm3Tvshvg_Gzg?pwd=86ii password: 86ii
    * GLO_L1OF_L1OC_real.bin
       -- for GLONASS L1OF
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +12e6 Hz
       -- samplingFreq: 53e6 Hz
    
       -- for GLONASS L1OC
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +10.995e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/1-fxXe6VW36izl2ve3dMQsw?pwd=j7p9 password: j7p9
    * GLO_L2OF_L2OC_real.bin
       -- for GLONASS L2OF
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -14e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for GLONASS L2OC
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -11.94e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/1P-xWe-aAUu3esymnbsPd-w?pwd=pnr8 password: pnr8
    * L2C_real.bin
       -- for GPS L2C (L2CM and L2CL)
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -7.4e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/16OZPphh6SAObSKf1t2ykPA?pwd=ss7k password: ss7k
    * E5b_B2b_GLO_L3OC_real.bin
       -- for Galileo E5b, BDS B2b
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +17.14e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for GLONASS L3OC
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +12.025e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/1Z-LK5UZ_DQgQk654cY0LgA?pwd=puzw password: puzw
    * L5_E5a_B2a_real.bin
       -- for GPS L5, Galileo E5a and BDS B2a
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: -13.55e6 Hz
       -- samplingFreq: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1xlftENWf76QwEVzzmxSUDw?pwd=m74w password: m74w
    * B3I_E6_real.bin
       -- for BDS B3I
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +8.52e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for Galileo E6B
       -- dataType: int8
       -- fileType: 1 (real samples)
       -- IF: +18.75e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/14q3K1MKNnAtxz64mUu9DKw?pwd=unn8 password: unn8
