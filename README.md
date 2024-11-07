# GDA: GNSS data analysis software

## Summary

GDA(GNSS Data Analysis)is a comprehensive software for **GNSS data file conversion and data analysis**. The core data processing module of the software is written in **Rust** with reference to the Rtklib variable framework and the Anubis XML file structure, the auxiliary decoding libraries are developed based on the Rtklib open source code, and the auxiliary drawing scripts are written in Python.The main functions of GDA include.

1. RTCM3 decode
2. RINEX/RTCM3 files post-processing quality analysis
3. NTRIP/ZMQ real-time data stream quality analysis
4. NRTIP real-time data stream reception

*Code packaging requires locally installed rust language and its package management compilation tool cargo, you can cargo build --release for the release version of the program packaging, third-party dependency libraries required will be automatically downloaded in the network environment (may need to be mirrored source configuration), packaged executables call description is as follows.*

## Package Catalog Structure

```tex
GDA
├── /decode                Decode Catalog
│   ├── /include           header file directory
│   ├── /src               rtklib source code and decoding library c files
│   └── CMakeLists.txt     CMake build configuration file
├── /qc                    Quality Analysis Catalog
│   ├── /src               source code directory
│   ├── Cargo.toml         Rust Project Configuration Files
│   └── decode.dll         shared library (computing)
├── /qcplot                Drawing Script Catalog
│   ├── qc_plot.py         quality analysis mapping scripts
│   └── test.py            test script
```

## Help Information

- Enter the command -h/--help to refer to the program call help, the program has the following functional modules
  - ***decode*** rtcm3 files decode
  - ***binstore*** ntrip real-time data stream reception
  - ***ntrip*** ntrip real-time data stream quality analysis
  - ***zmq*** ZMQ real-time data stream quality analysis
  - ***obsqc*** rinex files post-processing quality analysis
  - ***rtcmqc*** rtcm3 file post-processing quality analysis

```bash
D:\code\rust\gda\target\release>gda.exe -h
A GNSS toolbox for data processing like decoding and quality analysis.

Usage: gda.exe <COMMAND>

Commands:
  decode    Decode RTCM3 to RINEX
  binstore  Store Bin Data to File
  ntrip     QC by NTRIP Message to *File or Redis
  zmq       QC by ZMQ to File or *Redis
  obsqc     QC by RINEX File
  rtcmqc    QC by RTCM3 File
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help

If you have any questions, please reach out to: ghk40041@whu.edu.cn
```

## Decode

- decode -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe decode -h
Decode RTCM3 to RINEX

Usage: gda.exe decode [OPTIONS] --inp <IFILE>

Options:
  -i, --inp <IFILE>        Input RTCM3 File
  -o, --obs <OBSFILE>      Output OBS File
  -n, --nav <NAVFILE>      Output NAV File
  -s, --ts <TS>            Time of Start
  -e, --te <TE>            Time of End
  -t, --ti <TI>            Sample
  -r, --tr <TR>            Approximate Time
  -v, --version <VERSION>  Output File Version [default: 3.04]
  -y, --sys <NAVSYS>       System Selection
  -h, --help               Print help
```

### parameter description

- -i/--inp : Input file, rtcm3 raw file, no requirement for file extension, **must pass**
- -o/--obs : Output obs file, customize the output name, need to include the path, otherwise output in the relative path in the project directory
- -n/--nav : Output nav file, customize the output name, need to include the path, otherwise output in the relative path in the project directory

> Note that at least one output file should be specified; the program will exit if no output file is specified.

- -s/--ts : Start time, not mandatory, in the format of "yyyy/mm/dd hh:mm:ss", with semicolon quotes, can be used to filter data.
- -e/--te : End time, not mandatory, same format as start time, here is a specific example *-e "2024/01/11 00:00:00"*
- -r/--tr : Approximate time, need to be with the time of the file in the same GPS week, if not passed will automatically extract the current time of the program running
- -v/--version : Output file version, can not pass the default value of 3.04, other commonly used version 2.10/3.02
- -y/--sys : Satellite system selection, you can use the default value (GPS|GLO|GAL|BDS|QZS|IRN), if you want to specify the system, multiple systems are separated by a half comma, such as specifying the use of GPS and BDS, the parameter is -y G,C , the corresponding characters of the system are shown in the following table

| System name | Abbreviated character |
| :---------: | :-------------------: |
|     GPS     |           G           |
|   GLONASS   |           R           |
|   Galileo   |           E           |
|     BDS     |           C           |
|    QZSS     |           Q           |
|     IRN     |           I           |

## Binstore

- binstore -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe binstore -h
Store Bin Data to File

Usage: gda.exe binstore [OPTIONS] --host <HOST> --port <PORT> --mountpoint <MOUNTPOINT> --username <USERNAME> --password <PASSWORD>

Options:
  -s, --host <HOST>              NTRIP Server Hosting
  -p, --port <PORT>              NTRIP Server Port
  -m, --mountpoint <MOUNTPOINT>  NTRIP Mount Point
  -u, --username <USERNAME>      NTRIP Username
  -w, --password <PASSWORD>      NTRIP Password
  -o, --oup <OUTFILE>            Output File, defaults to the name of the mount point in the project directory
  -h, --help                     Print help
```

### parameter description

- -s/--host : ntrip server host, must pass
- -p/--port : ntrip server port, must pass
- -m/--mountpoint : ntrip mount point, must pass
- -u/--username : ntrip username
- -w/--password : ntrip password
- -o/--oup : The output file, if not passed, is saved in the project directory by default with the name of the mount point.

## Ntrip

- ntrip -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe ntrip -h
QC by NTRIP Message to *File or Redis

Usage: gda.exe ntrip [OPTIONS] --host <HOST> --port <PORT> --mountpoint <MOUNTPOINT> --username <USERNAME> --password <PASSWORD>

Options:
  -s, --host <HOST>              NTRIP Server Hosting
  -p, --port <PORT>              NTRIP Server Port
  -m, --mountpoint <MOUNTPOINT>  NTRIP Mount Point
  -u, --username <USERNAME>      NTRIP Username
  -w, --password <PASSWORD>      NTRIP Password
  -q, --qcpath <QCPATH>          QC File Output Path, saved in the current directory by default with the name of the mount point
  -i, --redisip <REDISIP>        Redis Server Address
  -o, --redisport <REDISPORT>    Redis Server Port
  -a, --auth <REDISAUTH>         Redis Server Authentication Commands
  -d, --addition <ADDITION>      Additional Documents for Ephemeris Supplement
  -h, --help                     Print help
```

### parameter description

- -s/--host : ntrip server host, must pass
- -p/--port : ntrip server port, must pass
- -m/--mountpoint : ntrip mount point, must pass
- -u/--username : ntrip username
- -w/--password : ntrip password
- -q/--qcpath : Analyze the output path of the file, the file is named after the name of the mount point, if you do not pass the default saved in the project directory
- -i/--redisip : redis server address
- -o/--redisport : redis server port
- -a/--auth : redis server authentication commands

> ntrip_qc mode is based on file output, redis-related parameters are not transmitted

- -d/--addition : Pass in a supplemental ephemeris file to be used as a substitute for ephemeris data if it has not yet been received, not mandatory

## Zmq

- zmq -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe zmq -h
QC by ZMQ to File or *Redis

Usage: gda.exe zmq [OPTIONS] --host <HOST> --port <PORT>

Options:
  -s, --host <HOST>            Subscription Address
  -p, --port <PORT>            Server Port
  -t, --stations <STATION>     Site Configuration
  -i, --redisip <REDISIP>      Redis Server Address
  -o, --redisport <REDISPORT>  Redis Server Port
  -a, --auth <REDISAUTH>       Redis Server Authentication Commands
  -q, --qcpath <QCPATH>        Default Output Path for Analysis Files
  -d, --addition <ADDITION>    Supplementary ephemeris
  -h, --help                   Print help
```

### parameter description

- -s/--host : Subscription address, must pass on
- -p/--port : Server port, must be passed
- -t/--stations : Subscribe to the site name, site name separated by a comma, if you do not pass the default subscription queue in all sites
- -i/--redisip : redis server address
- -o/--redisport : redis server port
- -a/--auth : redis server authentication commands
- -q/--qcpath : Analyze the output path of the file, the file is named after the name of the mount point, if you do not pass the default saved in the project directory

> zmq_qc mode is based on redis output, the file parameter may not be passed

- -d/--addition : Supplementary ephemeris files are passed in to replace ephemeris data for processing when it has not been received, and are not required.

## Obsqc

- obsqc -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe obsqc -h
QC by RINEX File

Usage: gda.exe obsqc [OPTIONS] --output <OFILE>

Options:
  -i, --input <IFILE>   Input file
  -o, --output <OFILE>  Output file
  -s, --ts <TS>         Time of Start
  -e, --te <TE>         Time of End
  -t, --ti <TI>         Approximate Time
  -q, --detail <Q>      Output Details (0, 1, 2) [default: 0]
  -h, --help            Print help
```

### parameter description

- -i/--input : RINEX input file, obs file must be passed in, can pass in both obs observation data file and nav ephemeris file, separated by half commas
- -o/--output : Quality analysis results output file, must be transferred
- -s/--ts : Start time, not mandatory, in the format of "yyyy/mm/dd hh:mm:ss", with semicolon quotes, can be used to filter data.
- -e/--te : End time, not mandatory, same format as start time
- -t/--ti : The time interval, i.e. the sampling rate, is not mandatory; when not passed in, the program will use the multitude of time intervals between the first 60 calendar elements as the sampling rate.
- -q/--detail : Detailed program parameters for quality analysis, set up at three levels as follows, default 0
  - 0 : Streamlining of the analysis so that the results contain only overall and systemic information
  - 1 : Detailed analysis, adding satellite information and signal information
  - 2 (Greater than 1 is sufficient) : Complete analysis for outputting data required for height angle distribution, sky plot, and time series plotting

## Rtcmqc

- rtcmqc -h/--help

```bash
D:\code\rust\gda\target\release>gda.exe rtcmqc -h
QC by RTCM3 File

Usage: gda.exe rtcmqc [OPTIONS] --input <IFILE> --output <OFILE>

Options:
  -i, --input <IFILE>   Input File
  -o, --output <OFILE>  Output File
  -s, --ts <TS>         Time of Start
  -e, --te <TE>         Time of End
  -t, --ti <TI>         Sample
  -r, --tr <TR>         Approximate Time
  -q, --detail <Q>      Output Details (0, 1, 2) [default: 0]
  -h, --help            Print help
```

### parameter description

- -i/--input : RTCM3 input file, must pass
- -o/--output : Quality analysis results output file, must be transferred
- -s/--ts : Start time, not mandatory, in the format of "yyyy/mm/dd hh:mm:ss", with semicolon quotes, can be used to filter data.
- -e/--te : End time, not mandatory, same format as start time
- -t/--ti : The time interval, i.e. the sampling rate, is not mandatory; when not passed in, the program will use the multitude of time intervals between the first 60 calendar elements as the sampling rate.
- -r/--tr : Approximate time, need to be with the time of the file in the same GPS week, if not passed will automatically extract the current time of the program running
- -q/--detail : Detailed program parameters for quality analysis, set up at three levels as follows, default 0
  - 0 : Streamlining of the analysis so that the results contain only overall and systemic information
  - 1 : Detailed analysis, adding satellite information and signal information
  - 2 (Greater than 1 is sufficient) : Complete analysis for outputting data required for height angle distribution, sky plot, and time series plotting

## QC Document Description

|                      Overall indicators                      |
| :----------------------------------------------------------: |
|     First_Epoch（Epoch of commencement of observations）     |
|         Last_Epoch（Observation termination epoch）          |
|             Hours（Observation duration）(hour)              |
|                  Sample（Sampling rate）(s)                  |
|             MinEle（Minimum elevation angle）(°)             |
|              Epoch_%Rt（Epoch completeness）(%)              |
|          Obs_%Rt（Observed signal completeness）(%)          |
| Obs_%Rt(>10)（Observed signal completeness/elevation angle greater than 10°）(%) |
|                 Slip（Number of cycle-slip）                 |
|                   o/slps（obs/cycle-slip）                   |
|             SPP_STD（SPP Positioning Error）(m)              |

|                      System indicators                       |
| :----------------------------------------------------------: |
|          ExpObs（Number of expected observations）           |
|           HavObs（Number of actual observations）            |
|          Obs_%Rt（Completeness of observations）(%)          |
|                 Sat（Number of satellites）                  |
|             Sig（Number of pseudorange signals）             |
| xCoSv（Number containing only single-frequency pseudoranges） |
|   xPhSv（Number containing only single-frequency phases）    |
|            Slip（Number of cycle-slip observed）             |
|                   o/slps（obs/cycle-slip）                   |
|            MPx（pseudorange multipath effect）(m)            |

|                     Satellite indicators                     |
| :----------------------------------------------------------: |
|               Obsnum（Number of observations）               |
|               Slipnum（Number of cycle-slip）                |
|      nGap（Number of interruptions greater than 10min）      |
|                   o/slps（obs/cycle-slip）                   |
|            Sat_%Rt（Satellite completeness）（%）            |
| Sat_%Rt>10（Satellite completeness/elevation angle greater than 10°）（%） |
|            MPx（pseudorange multipath effect）(m)            |

|                       Signal indicator                       |
| :----------------------------------------------------------: |
|                       Freq.Band（band)                       |
|                     CODE（signal code）                      |
|          ExpObs（Number of expected observations）           |
|           HavObs（Number of actual observations）            |
|                 %Ratio（completeness）（%）                  |
| Exp>10（Number of expected observations/elevation angle greater than 10°） |
| Hav>10（Number of actual observations/height angle greater than 10°） |
|  %Rt>10（Completeness/height angle greater than 10°）（%）   |
|                SNR（Carrier-to-noise ratio）                 |

|        Remaining information         |
| :----------------------------------: |
|         prn（Satellite PRN）         |
| epochnum（calendar element number）  |
|    azimuth（azimuth angle）（°）     |
|  elevation（elevation angle）（°）   |
| obsflag（observation start markers） |
|    slipflag（cycle-slip markers）    |



## Decoder Library

Program by the functional modules involved in rtcm raw file processing will be used **decoding dynamic link library** , the library by  **rtklib open source code for secondary development** access to the detailed code in the open source directory, according to the system requirements for their own compilation, in the program through the cargo for packaging, will be **dynamic link library placed in the src file under the same level of the directory under**, the program will be The program will be packaged with the dynamic library content, the program will be run for the first time in the project directory to create a lib folder and release the creation of the previously stored decoding libraries for decoding calls.

## QC Document Mapping

The results of the post hoc quality analysis are output in the form of a .qc file. The python plotting code is given in the open source directory, and the plotting code can be packaged as an executable file by using the third-party python packaging library pyinstaller, which is a program called with the following instructions

- Enter the command -h/--help to see the program call help.

```bash
D:\code\python\qcplot\dist>qc_plot_pro.exe -h
usage: qc_plot_pro.exe [-h] [-f {pos,sys,sig,skyel,ts,all}] [-i INPUT] [-o OUTPUT]

Plot qc_file

options:
  -h, --help            show this help message and exit
  -f {pos,sys,sig,skyel,ts,all}, --function {pos,sys,sig,skyel,ts,all}
                        Function to plot
  -i INPUT, --input INPUT
                        Input qc_file
  -o OUTPUT, --output OUTPUT
                        Output png path
```

### parameter description

- -f/--function : Drawing type, must be passed, type is as follows
  - pos : dops distribution map and SPP coordinate sequence
  - sys : Statistical chart of quality information by system
  - sig : Signal quality information distribution map
  - skyel : Satellite sky map and altitude angle distribution
  - ts : chronology
  - all : All diagrams above
- -i/--input_file : Input qc file path, must pass
- -o/--output_file_path : Output png image save path

### running example

```bash
D:\code\rust\qc\target\release>qc.exe rtcmqc -i XNXA.2024197binRTCM3 -o XNXA.qc -r "2024/07/14 00:00:00" -q 2
input file:XNXA.2024197binRTCM3
output file:XNXA.qc
/ [00:00:03] Finish scanning, totally 49494 epochs, 2340409 sat obs and 1234 sat nav                                          	[00:00:26] [########################################] 49494/49494 (0s) QC completed
D:\code\rust\qc\target\release>qc_plot.exe -f all -i XNXA.qc -o .\
finish plot
```

- The example rtcm3 file XNXA.2024197binRTCM3 was analyzed for quality after the first act to obtain a qc file, followed by a plotting process to obtain the picture :.
  - Chronology XNXA_ts.png
  - Statistical chart of system quality information XNXA_sys.png
  - Signal Quality Information Statistical Chart XNXA_sig.png
  - Satellite Sky Map XNXA_skymap.png
  - Elevation Angle Distribution XNXA_elmap.png
  - Dops Distribution Map XNXA_dops.png
  - SPP Coordinate Sequence Diagram XNXA_pos.png

> Images are auto-named site_f.png with site name site and drawing type f

