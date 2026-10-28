An Open-Source MATLAB GNSS SDR Toolbox with SIMD/GPU Acceleration: mSoftGNSS 
===============================================================================



Overview
-------------------------------------------------------------------------------
mSoftGNSS 是一个基于MATLAB的开源工具箱，用于对录制的GNSS中频（IF）信号进行后处理。
该工具箱基于SoftGNSS接收机架构，在可配置的处理框架下支持GPS、Galileo、GLONASS和北斗信号。
支持CPU/GPU捕获、数据与导频跟踪、导航电文解码、伪距计算和定位，具体功能因信号而异。
工具箱提供通道串行和通道并行两种跟踪模式，每种模式均支持MATLAB相关器、SIMD加速相关器和
GPU加速相关器。SIMD包含x86 AVX2和ARM64 NEON两种实现，由程序根据运行平台自动选择。
通道并行跟踪使各通道共用同一中频数据缓冲区，批量处理相关任务，减少重复读取。
计算密集型相关器由C++和CUDA C++实现，并通过MEX接口集成；接收机调度、跟踪环路控制、
导航电文解码和定位仍在MATLAB中完成，兼顾算法开发的便利性与计算效率。
工具箱还支持现代化GNSS信号和LDPC解码，可用于算法研究、接收机原型开发和基于录制中频数据的评估。



Authors
-------------------------------------------------------------------------------
* Yafeng Li
    * E-Mail: <lyf8118@126.com>
    * Wechat: lyf8118521
    * GNSS软件接收机技术讨论QQ群：147304049


* Dennis Akos  
    * E-Mail: <dma@colorado.edu>
    * HP: <http://www.colorado.edu/aerospace/dennis-akos>




Features
-------------------------------------------------------------------------------
* 基于 MATLAB 实现的 GNSS 信号处理功能
    * 本地码生成
    * CPU 和 GPU 捕获
    * 通道串行和通道并行跟踪
    * MATLAB、SIMD MEX 和 CUDA MEX 相关器
    * 现代化信号的数据/导频跟踪
    * 载噪比（C/N0）估计和跟踪结果绘图
    * 导航电文解码（包括 LDPC 解码器）
    * 伪距生成
    * 接收机位置与钟差解算及结果绘图
* SIMD/GPU 加速和相关器选项
    * settings.gpuACQflag = 0/1 分别选择 CPU/GPU 捕获。
    * settings.trkMode = 0/1 分别选择通道串行/通道并行跟踪。
    * settings.correlatorType = 0 表示使用 MATLAB 脚本相关器。
    * settings.correlatorType = 1 表示使用 CPU SIMD MEX 相关器。
    * settings.correlatorType = 2 表示使用 CUDA GPU MEX 相关器。
    * MATLAB 相关器包含 BPSK、QPSK、TMBPSK、QMBOC等变体，
      位于各接收机的 include 目录中。
    * 共用 MEX 相关器包含 BPSK、QPSK 和 QMBOC 的
      SIMD 和 CUDA 实现。
* 支持的信号
    * GPS L1 C/A （数据）
    * GPS L1C （数据与导频）
    * GPS L2C （数据与导频）
    * GPS L5 （数据与导频）
    * Galileo E1 （数据与导频）
    * Galileo E5a （数据与导频）
    * Galileo E5b （数据与导频）
    * Galileo E6 （数据与导频；HAS 页解码）
    * GLONASS L1OF （数据）
    * GLONASS L2OF （数据）
    * GLONASS L1OC （数据与导频）
    * GLONASS L2OCp （导频）
    * GLONASS L3OC （数据与导频）
    * 北斗 B1I （数据）
    * 北斗 B3I （数据）
    * BDS-3 B1C （数据与导频）
    * BDS-3 B2a （数据与导频）
    * 北斗 B2b （数据）
* 支持以二进制文件存储的中频采样数据
    * settings.dataType 可设为 int8 或 int16，两种类型均支持实数和复数采样。
    * settings.fileType = 1 表示实数采样，2 表示 I/Q 交织复数采样
      （I0,Q0,I1,Q1,...）；dataType 指定每个 I/Q 分量的数据类型。
    * 下文提供 LimeSDR 复采样和 NUT4NT 实采样数据及配置示例。
      使用其他数据时，应根据录制参数设置数据格式、中频和采样率。
    * 接收机已在 Windows 下的 MATLAB R2025a/R2025b 中进行测试。



Directory and Files
-------------------------------------------------------------------------------
根目录

    ./Doc                     文档、ICD、论文和各接收机说明材料。
    ./IF_Data_Set             可放置待处理 IF 数据及对应的元数据文件。
    ./native_Correlators      所有接收机共用的 SIMD/CUDA MEX 相关器源码和
                              在本机编译的 MEX 文件。

接收机目录

    每个接收机目录都是相对独立的软件接收机；除 ./native_Correlators、
    ./IF_Data_Set 和 ./Doc 外，各接收机代码互不依赖。当前接收机目录如下：

    ./mGPS_L1CA               GPS L1 C/A 软件接收机
    ./mGPS_L1C                GPS L1C 软件接收机
    ./mGPS_L2C                GPS L2C 软件接收机
    ./mGPS_L5                 GPS L5 软件接收机
    ./mGalileo_E1             Galileo E1 软件接收机
    ./mGalileo_E5a            Galileo E5a 软件接收机
    ./mGalileo_E5b            Galileo E5b 软件接收机
    ./mGalileo_E6             Galileo E6B/E6C 软件接收机（HAS 接收）
    ./mGlonass_L1_L2          GLONASS L1OF/L2OF 软件接收机
    ./mGlonass_L1OC           GLONASS L1OC 软件接收机
    ./mGlonass_L2OC           GLONASS L2OCp 软件接收机
    ./mGlonass_L3OC           GLONASS L3OC 软件接收机
    ./mBDS_B1I                北斗B1I软件接收机
    ./mBDS_B3I                北斗B3I软件接收机
    ./mBDS-3_B1C              北斗三号B1C软件接收机
    ./mBDS-3_B2a              北斗三号B2a软件接收机
    ./mBDS_B2b                北斗B2b软件接收机

每个接收机目录的通用结构

    ./init.m                  接收机启动脚本：设置路径、读取配置、显示中频信号的
                              时域波形和频谱，经用户确认后启动处理。
    ./initSettings.m          当前信号接收机的参数配置文件。
    ./postProcessing.m        顶层处理脚本：调度捕获、通道初始化、跟踪、导航解算、
                              结果保存和绘图。
    ./acquisitionCPU.m        CPU/MATLAB 捕获实现。
    ./acquisitionGPU.m        MATLAB GPU-array 捕获实现。
    ./trkChannelsSerial.m     通道串行跟踪实现。
    ./trkChannelsParallel.m   通道并行跟踪实现。
    ./postNavigation.m        导航电文解码和定位解算。
                              ./mGlonass_L3OC 的该文件位于 ./include/postNavigation.m。
    ./include                 与当前信号相关的辅助函数、本地码生成、MATLAB 相关器、
                              导航电文解码、绘图和通道初始化函数。
                              请保留 mGalileo_E1/include 中的 E1b.dat/E1c.dat，
                              以及 mGalileo_E6/include 中的 HAS Reed-Solomon 生成矩阵。
    ./include/ldpcDecoder     LDPC 解码相关函数；仅存在于需要 LDPC 导航解码的接收机，
                              例如 mBDS-3_B1C、mBDS-3_B2a、mBDS_B2b 和 mGPS_L1C。
    ./Common                  本接收机使用的通用数学与定位函数，例如伪距计算、
                              坐标转换、最小二乘定位和环路参数计算。

./native_Correlators 中的共用相关器

    corrX86SIMDSerial*.cpp       SIMD 通道串行跟踪相关器源码。
    corrX86SIMDParallel*.cpp     SIMD 通道并行跟踪相关器源码。
    corrNeonSIMD*.cpp            ARM64 NEON 通道串行/并行相关器源码。
    corrGPUSerial*.cu            CUDA 通道串行跟踪相关器源码。
    corrGPUParallel*.cu          CUDA 通道并行跟踪相关器源码。
    corrGPUParallelFused*.cu     CUDA 融合内核通道并行跟踪相关器源码。
    *.mexw64 / *.mexa64          Windows x64 / Linux x86-64 的 MEX 文件。
    *.mexmaca64                  Apple silicon 原生 MEX 文件（NEON）。

相关器选择

    settings.correlatorType = 0
                              使用各接收机 ./include 目录中的 MATLAB 脚本相关器，
                              便于阅读和调试。
    settings.correlatorType = 1
                              使用 ./native_Correlators 中的 CPU SIMD MEX 相关器，
                              程序根据当前 MATLAB 架构自动选择 x86 或 ARM64 NEON，
                              不依赖 CUDA 即可加速跟踪。
    settings.correlatorType = 2
                              使用 ./native_Correlators 中的 CUDA GPU MEX 相关器，
                              用于高吞吐的通道串行或通道并行跟踪。

Software Dependencies
-------------------------------------------------------------------------------
* 需要 MATLAB。MEX 文件应针对当前 MATLAB 所在操作系统和处理器架构编译。
* IF 频谱绘图使用 MATLAB Signal Processing Toolbox 中的 pwelch。
* 导航电文解码按信号需要使用 MATLAB Communications Toolbox 中的部分函数：
  -- 创建有限域数组：gf()
  -- BCH 解码：bchdec()
  -- 创建 CRC 检测对象：comm.CRCDetector()
  -- 将卷积码生成多项式转换为网格结构：poly2trellis()
  -- Viterbi 解码：vitdec()
  -- 执行 CRC 检测：step()
* GPU-array 捕获需要 MATLAB Parallel Computing Toolbox 和支持 CUDA 的 GPU。
  直接运行 CUDA MEX 相关器需要兼容的 NVIDIA GPU、驱动和 CUDA 运行库。
* x86 SIMD MEX 需要支持AVX2和FMA的CPU；ARM64使用独立的NEON源码。
  编译需要 MATLAB 支持的 C++ 编译器。
  重新编译 CUDA MEX 需要 CUDA Toolkit、兼容的编译器和 mexcuda
  （Parallel Computing Toolbox）。
* 如果只需要 CPU 运行，可在各接收机 initSettings.m 中关闭 GPU 捕获，并选择
  MATLAB 相关器或 SIMD 相关器。
  
  

Compile correlators
-------------------------------------------------------------------------------
本节说明 mSoftGNSS 的 MATLAB MEX 编译方法，请按操作系统和 CPU 架构选择。
命令在 MATLAB 命令窗口执行，初始目录为 mSoftGNSS 根目录。进入本工程的
native_Correlators 目录后，按对应平台章节编译 SIMD 和 CUDA 相关器。

下文的 cudaArch = '89' 以 RTX 4080 Laptop GPU 为例。换用其他 GPU 时，
应按其计算能力修改，并确认所用 CUDA Toolkit 支持该目标。

|   平台  | SIMD 要求 | 编译产物 |
| ------ | -------- | ------ |
| Windows    x86-64 | 支持 AVX2 和 FMA 的 x86-64 CPU | .mexw64 |
| Linux      x86-64 | 支持 AVX2 和 FMA 的 x86-64 CPU | .mexa64 |

### Windows x86-64

SIMD 以 MinGW-w64/G++ 为例。在 mex -setup C++ 中选择 MATLAB 支持的MinGW64 Compiler (C++)，编译 6 个 x86 SIMD 相关器。

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

然后切换到 MATLAB 和 CUDA Toolkit 支持的 Microsoft Visual C++，编译 CUDA 相关器。

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

完成 SIMD 和 CUDA 编译后，本目录中生成 .mexw64 文件。

### Linux x86-64

使用 MATLAB 支持的 GCC/G++，先将其加入启动 MATLAB 时的 PATH，再通过mex -setup C++ 配置。
先编译x86 SIMD相关器：

```matlab
cd('native_Correlators');
mex -setup C++
sources = dir('corrX86SIMD*.cpp');
for k = 1:numel(sources)
    mex('-R2017b', '-O', ...
        'CXXFLAGS=$CXXFLAGS -std=c++14 -mavx2 -mfma', sources(k).name);
end
```

然后在同一目录编译 CUDA 相关器。主机 GCC/G++ 必须同时受到 MATLAB 和 CUDA
Toolkit 支持。

```matlab
cudaArch = '89';
cudaFlags = sprintf('NVCCFLAGS=-gencode=arch=compute_%s,code=sm_%s', ...
    cudaArch, cudaArch);
sources = dir('corrGPU*.cu');
for k = 1:numel(sources)
    mexcuda('-R2017b', '-O', cudaFlags, sources(k).name);
end
```

完成 SIMD 和 CUDA 编译后，本目录中生成 .mexa64 文件。

生成的 MEX 文件应与对应源文件同名，并放在本工程 native_Correlators 目录中，可在 MATLAB 命令窗口检查：

```matlab
mexext
dir(['*.' mexext])
```

How to use
-------------------------------------------------------------------------------
* Step 1: 进入目标信号对应的接收机目录，例如 "mBDS-3_B1C" 或 "mGPS_L1CA"。
* Step 2: 将 IF 数据放入 "IF_Data_Set"，或在 settings.fileName 中设置数据文件
          的绝对路径。
* Step 3: 在 "initSettings.m" 中配置处理时长、IF 数据、捕获、跟踪、相关器和
          导航解算参数。skipNumberOfSamples 按采样点计数，一个实数值或一对 I/Q
          分量计为一个采样点。
* Step 4: 通过 settings.gpuACQflag 选择 CPU/GPU 捕获，通过 settings.trkMode
          选择通道串行/并行跟踪，通过 settings.correlatorType 选择脚本、SIMD
          或 GPU 相关器。
* Step 5: 在当前接收机目录下运行 "init.m"，查看中频数据图后，
          输入 1 开始处理，或输入 0 退出。
* Step 6: 处理完成后，查看捕获、跟踪和导航定位结果及相关图形。



Implementation details
-------------------------------------------------------------------------------
本软件接收机套件沿用了 SoftGNSS 风格的 MATLAB SDR 组织方式，并扩展了多种
现代化 GNSS 信号、CPU/GPU 捕获、通道串行/并行跟踪，以及共用 SIMD/CUDA MEX
相关器。

通道并行跟踪以历元为外层循环，通过共享中频数据缓冲区和批量相关器调用，集中完成
各活动通道的相关运算，再分别更新各通道的跟踪环路。该架构在保留各通道独立状态
与环路反馈的同时，减少中频数据的重复读取和接口调用开销，并支持 SIMD/GPU
多通道批量加速。

BDS-3 B1C/B2a 接收机设计可参考：
Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and implementation of an open-source
BDS-3 B1C/B2a SDR receiver. GPS Solut (2019) 23: 60.
https://doi.org/10.1007/s10291-019-0853-z


​	   
Test signal (collected by LimeSDR) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C.bin
       -- for GPS L1 C/A, Galileo E1 and BDS B1C 
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1NdZeVrEdEjDnJiV91Cg0qg?pwd=bvh5 password: bvh5
    * GPS L1C.iq
       -- for GPS L1C 
       -- dataType: int8
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int8）
       -- IF: 0 Hz
       -- samplingFreq: 25e6 Hz
       -- link: https://pan.baidu.com/s/1hZ4S8R4VWB5OtPr9VyPCYQ?pwd=tv8b password: tv8b
    * GPS_L2.bin
       -- for GPS L2C
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1KqdyWyuVPAgx0EoUdoLwag?pwd=5tya password: 5tya
    * L5_E5a_B2a.bin
       -- for GPS L5, Galileo E5a and BDS B2a
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1rX2lguJgMmgjAdZIJK1FNg?pwd=bnzm password: bnzm
    
    * E5b_B2b.bin   
       -- for Galileo E5b, BDS B2b
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/128dwHIKUG9iM0MigS_7vJw?pwd=x4uq password: x4uq
    * GAL_E6.bin  
       -- for Galileo E6
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1OrHtpNG-dQxSmmVc9rX6tA?pwd=qfx8 password: qfx8
    * GLO_L1OC_L1.bin 
       -- for GLONASS L1OF  
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 1.005e6 Hz
       -- samplingFreq: 30e6 Hz
       -- GLONASS L1OC
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1w-fYtmM3oW_pxR4ZH_HGyw?pwd=js93 password: js93
    * GLO_L2OC_L2.bin  
       -- for GLONASS L2OF 
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: -2.06e6 Hz
       -- samplingFreq: 30e6 Hz
    
       -- for GLONASS L2OC 
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1jefXiCkeI9qHSqwZ_RXr2g?pwd=33w6 password: 33w6
    * GLO_L3OC.bin  
       -- for GLONASS L3OC	
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1ySNHn3W66Q6n9fxERR7nOQ?pwd=289m password: 289m
    * BDS_B1I.bin
       -- for BDS B1I 
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1qV-TCtROGHL8iJZlqw5mZA?pwd=39vd password: 39vd
    * BDS_B3I.bin       
       -- for BDS B3I    
       -- dataType: int16
       -- fileType: 2（I/Q 交织复数采样，每个分量为 int16）
       -- IF: 0 Hz
       -- samplingFreq: 30e6 Hz
       -- link: https://pan.baidu.com/s/1T_g1to1e8OW5OFtHlQsZeg?pwd=6w1u password: 6w1u 



Test signal (collected by NUT4NT) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C_real.bin
       -- for GPS L1 C/A, Galileo E1 and BDS B1C
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -14.58e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/1cFGCBwcWA0vHBeiL7WU_2Q?pwd=r8au password: r8au
    * B1I_real.bin
       -- for BDS B1I
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -28.902e6 Hz
       -- samplingFreq: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1y4LfFZSO-Qm3Tvshvg_Gzg?pwd=86ii password: 86ii
    * GLO_L1OF_L1OC_real.bin
       -- for GLONASS L1OF
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +12e6 Hz
       -- samplingFreq: 53e6 Hz
    
       -- for GLONASS L1OC
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +10.995e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/1-fxXe6VW36izl2ve3dMQsw?pwd=j7p9 password: j7p9
    * GLO_L2OF_L2OC_real.bin
       -- for GLONASS L2OF
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -14e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for GLONASS L2OC
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -11.94e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/1P-xWe-aAUu3esymnbsPd-w?pwd=pnr8 password: pnr8
    * L2C_real.bin
       -- for GPS L2C (L2CM and L2CL)
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -7.4e6 Hz
       -- samplingFreq: 53e6 Hz
       -- link: https://pan.baidu.com/s/16OZPphh6SAObSKf1t2ykPA?pwd=ss7k password: ss7k
    * E5b_B2b_GLO_L3OC_real.bin
       -- for Galileo E5b, BDS B2b
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +17.14e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for GLONASS L3OC
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +12.025e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/1Z-LK5UZ_DQgQk654cY0LgA?pwd=puzw password: puzw
    * L5_E5a_B2a_real.bin
       -- for GPS L5, Galileo E5a and BDS B2a
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: -13.55e6 Hz
       -- samplingFreq: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1xlftENWf76QwEVzzmxSUDw?pwd=m74w password: m74w
    * B3I_E6_real.bin
       -- for BDS B3I
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +8.52e6 Hz
       -- samplingFreq: 70e6 Hz
    
       -- for Galileo E6B
       -- dataType: int8
       -- fileType: 1（实数采样）
       -- IF: +18.75e6 Hz
       -- samplingFreq: 70e6 Hz
       -- link: https://pan.baidu.com/s/14q3K1MKNnAtxz64mUu9DKw?pwd=unn8 password: unn8
