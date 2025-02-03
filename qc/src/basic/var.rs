use crate::basic::func::sqvar;
use crate::basic::time::time2str;
use std::io::{Write, Result};
use crate::convbin::c_var;

// 常量定义
pub const MAXANT: usize = 64;
pub const MAXSAT: usize = 221;
pub const NFREQ: usize = 3;
pub const NEXOBS: usize = 3;
pub const NFREOBS: usize = NFREQ + NEXOBS;
pub const MAXOBSTYPE: usize = 64;
pub const MAXRNXLEN: usize = 16 * MAXOBSTYPE + 4;
pub const MAXCOMMENT: usize = 100;
pub const MAXFILE: usize = 16;
pub const MAXINFILE: usize = 16;
pub const MAXPOSHEAD: usize = 1024;
pub const MAXOBS: usize = 96;
pub const MAXLEAPS: usize = 64;
pub const NINCOBS: usize = 24 * 3600 * 3;
pub const MAXRCV: usize = 64;
pub const MAXNIGP: usize = 201;
pub const MAXBAND: usize = 10;
pub const CODE_NONE: u8 = 0; // obs code: none or unknown
pub const CODE_L1C: u8 = 1; // obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS)
pub const CODE_L1P: u8 = 2; // obs code: L1P,G1P,B1P (GPS,GLO,BDS)
pub const CODE_L1W: u8 = 3; // obs code: L1 Z-track (GPS)
pub const CODE_L1Y: u8 = 4; // obs code: L1Y (GPS)
pub const CODE_L1M: u8 = 5; // obs code: L1M (GPS)
pub const CODE_L1N: u8 = 6; // obs code: L1codeless,B1codeless (GPS,BDS)
pub const CODE_L1S: u8 = 7; // obs code: L1C(D) (GPS,QZS)
pub const CODE_L1L: u8 = 8; // obs code: L1C(P) (GPS,QZS)
pub const CODE_L1E: u8 = 9; // (not used)
pub const CODE_L1A: u8 = 10; // obs code: E1A,B1A (GAL,BDS)
pub const CODE_L1B: u8 = 11; // obs code: E1B (GAL)
pub const CODE_L1X: u8 = 12; // obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS)
pub const CODE_L1Z: u8 = 13; // obs code: E1A+B+C,L1S (GAL,QZS)
pub const CODE_L2C: u8 = 14; // obs code: L2C/A,G1C/A (GPS,GLO)
pub const CODE_L2D: u8 = 15; // obs code: L2 L1C/A-(P2-P1) (GPS)
pub const CODE_L2S: u8 = 16; // obs code: L2C(M) (GPS,QZS)
pub const CODE_L2L: u8 = 17; // obs code: L2C(L) (GPS,QZS)
pub const CODE_L2X: u8 = 18; // obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS)
pub const CODE_L2P: u8 = 19; // obs code: L2P,G2P (GPS,GLO)
pub const CODE_L2W: u8 = 20; // obs code: L2 Z-track (GPS)
pub const CODE_L2Y: u8 = 21; // obs code: L2Y (GPS)
pub const CODE_L2M: u8 = 22; // obs code: L2M (GPS)
pub const CODE_L2N: u8 = 23; // obs code: L2codeless (GPS)
pub const CODE_L5I: u8 = 24; // obs code: L5I,E5aI (GPS,GAL,QZS,SBS)
pub const CODE_L5Q: u8 = 25; // obs code: L5Q,E5aQ (GPS,GAL,QZS,SBS)
pub const CODE_L5X: u8 = 26; // obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS)
pub const CODE_L7I: u8 = 27; // obs code: E5bI,B2bI (GAL,BDS)
pub const CODE_L7Q: u8 = 28; // obs code: E5bQ,B2bQ (GAL,BDS)
pub const CODE_L7X: u8 = 29; // obs code: E5bI+Q,B2bI+Q (GAL,BDS)
pub const CODE_L6A: u8 = 30; // obs code: E6A,B3A (GAL,BDS)
pub const CODE_L6B: u8 = 31; // obs code: E6B (GAL)
pub const CODE_L6C: u8 = 32; // obs code: E6C (GAL)
pub const CODE_L6X: u8 = 33; // obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS)
pub const CODE_L6Z: u8 = 34; // obs code: E6A+B+C,L6D+E (GAL,QZS)
pub const CODE_L6S: u8 = 35; // obs code: L6S (QZS)
pub const CODE_L6L: u8 = 36; // obs code: L6L (QZS)
pub const CODE_L8I: u8 = 37; // obs code: E5abI (GAL)
pub const CODE_L8Q: u8 = 38; // obs code: E5abQ (GAL)
pub const CODE_L8X: u8 = 39; // obs code: E5abI+Q,B2abD+P (GAL,BDS)
pub const CODE_L2I: u8 = 40; // obs code: B1_2I (BDS)
pub const CODE_L2Q: u8 = 41; // obs code: B1_2Q (BDS)
pub const CODE_L6I: u8 = 42; // obs code: B3I (BDS)
pub const CODE_L6Q: u8 = 43; // obs code: B3Q (BDS)
pub const CODE_L3I: u8 = 44; // obs code: G3I (GLO)
pub const CODE_L3Q: u8 = 45; // obs code: G3Q (GLO)
pub const CODE_L3X: u8 = 46; // obs code: G3I+Q (GLO)
pub const CODE_L1I: u8 = 47; // obs code: B1I (BDS) (obsolete)
pub const CODE_L1Q: u8 = 48; // obs code: B1Q (BDS) (obsolete)
pub const CODE_L5A: u8 = 49; // obs code: L5A SPS (IRN)
pub const CODE_L5B: u8 = 50; // obs code: L5B RS(D) (IRN)
pub const CODE_L5C: u8 = 51; // obs code: L5C RS(P) (IRN)
pub const CODE_L9A: u8 = 52; // obs code: SA SPS (IRN)
pub const CODE_L9B: u8 = 53; // obs code: SB RS(D) (IRN)
pub const CODE_L9C: u8 = 54; // obs code: SC RS(P) (IRN)
pub const CODE_L9X: u8 = 55; // obs code: SB+C (IRN)
pub const CODE_L1D: u8 = 56; // obs code: B1D (BDS)
pub const CODE_L5D: u8 = 57; // obs code: L5D(L5S),B2aD (QZS,BDS)
pub const CODE_L5P: u8 = 58; // obs code: L5P(L5S),B2aP (QZS,BDS)
pub const CODE_L5Z: u8 = 59; // obs code: L5D+P(L5S) (QZS)
pub const CODE_L6E: u8 = 60; // obs code: L6E (QZS)
pub const CODE_L7D: u8 = 61; // obs code: B2bD (BDS)
pub const CODE_L7P: u8 = 62; // obs code: B2bP (BDS)
pub const CODE_L7Z: u8 = 63; // obs code: B2bD+P (BDS)
pub const CODE_L8D: u8 = 64; // obs code: B2abD (BDS)
pub const CODE_L8P: u8 = 65; // obs code: B2abP (BDS)
pub const CODE_L4A: u8 = 66; // obs code: G1aL1OCd (GLO)
pub const CODE_L4B: u8 = 67; // obs code: G1aL1OCd (GLO)
pub const CODE_L4X: u8 = 68; // obs code: G1al1OCd+p (GLO)
pub const MAXCODE: u8 = 68; // max number of obs code

pub const PI: f64 = std::f64::consts::PI;
pub const D2R: f64 = PI / 180.0;
pub const R2D: f64 = 180.0/PI;
pub const PMODE_SINGLE: i32 = 0;

pub const MAX_VAR_EPH: f64 = 300.0 * 300.0;
pub const NAVEXP:&str = "D";

pub const SYS_GPS: usize = 1 << 0;
pub const SYS_GLO: usize = 1 << 1;
pub const SYS_GAL: usize = 1 << 2;
pub const SYS_QZS: usize = 1 << 3;
pub const SYS_CMP: usize = 1 << 4;
pub const SYS_IRN: usize = 1 << 5;
pub const SYS_NONE: usize = 0;
pub const SYS_ALL: usize = (1 << 6) - 1;
pub const NUMSYS: usize = 6;

pub const TSYS_GPS: i32 = 0;
pub const TSYS_UTC: i32 = 1;
pub const TSYS_GLO: i32 = 2;
pub const TSYS_GAL: i32 = 3;
pub const TSYS_QZS: i32 = 4;
pub const TSYS_CMP: i32 = 5;
pub const TSYS_IRN: i32 = 6;

pub const NSATGPS: usize = 32;
pub const NSATGLO: usize = 27;
pub const NSATGAL: usize = 36;
pub const NSATQZS: usize = 10;
pub const NSATCMP: usize = 63;
pub const NSATIRN: usize = 10;
pub const MINPRNGPS: usize = 1;
pub const MAXPRNGPS: usize = 32;
pub const MINPRNGLO: usize = 1;
pub const MAXPRNGLO: usize = 27;
pub const MINPRNGAL: usize = 1;
pub const MAXPRNGAL: usize = 36;
pub const MINPRNQZS: usize = 193;
pub const MAXPRNQZS: usize = 202;
pub const MINPRNCMP: usize = 1;
pub const MAXPRNCMP: usize = 63;
pub const MINPRNIRN: usize = 1;
pub const MAXPRNIRN: usize = 14;

pub const MAXFREQ: usize = 7;
pub const FREQ1: f64 = 1.57542E9; // L1/E1/B1C frequency (Hz)
pub const FREQ2: f64 = 1.22760E9; // L2 frequency (Hz)
pub const FREQ5: f64 = 1.17645E9; // L5/E5a/B2a frequency (Hz)
pub const FREQ6: f64 = 1.27875E9; // E6/L6 frequency (Hz)
pub const FREQ7: f64 = 1.20714E9; // E5b frequency (Hz)
pub const FREQ8: f64 = 1.191795E9; // E5a+b frequency (Hz)
pub const FREQ9: f64 = 2.492028E9; // S frequency (Hz)
pub const FREQ1_GLO: f64 = 1.60200E9; // GLONASS G1 base frequency (Hz)
pub const DFRQ1_GLO: f64 = 0.56250E6; // GLONASS G1 bias frequency (Hz/n)
pub const FREQ2_GLO: f64 = 1.24600E9; // GLONASS G2 base frequency (Hz)
pub const DFRQ2_GLO: f64 = 0.43750E6; // GLONASS G2 bias frequency (Hz/n)
pub const FREQ3_GLO: f64 = 1.202025E9; // GLONASS G3 frequency (Hz)
pub const FREQ1A_GLO: f64 = 1.600995E9; // GLONASS G1a frequency (Hz)
pub const FREQ2A_GLO: f64 = 1.248060E9; // GLONASS G2a frequency (Hz)
pub const FREQ1_CMP: f64 = 1.561098E9; // BDS B1I frequency (Hz)
pub const FREQ2_CMP: f64 = 1.20714E9; // BDS B2I/B2b frequency (Hz)
pub const FREQ3_CMP: f64 = 1.26852E9; // BDS B3 frequency (Hz)
pub const MINFREQ_GLO: i32 = -7;
pub const MAXFREQ_GLO: i32 = 13;

pub const IONOOPT_OFF: i32 = 0;
pub const IONOOPT_BRDC: i32 = 1;
pub const IONOOPT_SBAS: i32 = 2;
pub const IONOOPT_IFLC: i32 = 3;
pub const IONOOPT_EST: i32 = 4;
pub const IONOOPT_TEC: i32 = 5;
pub const IONOOPT_QZS: i32 = 6;
pub const IONOOPT_STEC: i32 = 8;
pub const IONOOPT_POLY: i32 = 9;

pub const SNR_UNIT: f64 = 0.001;

pub const DTTOL: f64 = 0.025;

pub const LLI_SLIP: i32 = 0x01;
pub const LLI_HALFC: i32 = 0x02;
pub const LLI_BOCTRK: i32 = 0x04;

pub const OMGE: f64 = 7.2921151467e-5;
pub const CLIGHT: f64 = 299792458.0;

pub const MODE_OBSQC: i32 = 1;
pub const MODE_RTCMQC: i32 = 2;

pub const OBS_CODES: [&str; 70] = [
    "", "1C", "1P", "1W", "1Y", "1M", "1N", "1S", "1L", "1E", // 0-9
    "1A", "1B", "1X", "1Z", "2C", "2D", "2S", "2L", "2X", "2P", // 10-19
    "2W", "2Y", "2M", "2N", "5I", "5Q", "5X", "7I", "7Q", "7X", // 20-29
    "6A", "6B", "6C", "6X", "6Z", "6S", "6L", "8L", "8Q", "8X", // 30-39
    "2I", "2Q", "6I", "6Q", "3I", "3Q", "3X", "1I", "1Q", "5A", // 40-49
    "5B", "5C", "9A", "9B", "9C", "9X", "1D", "5D", "5P", "5Z", // 50-59
    "6E", "7D", "7P", "7Z", "8D", "8P", "4A", "4B", "4X", "", // 60-69
];

pub const CODEPRIS: [[&str; MAXFREQ]; 6] = [
    [
        /* GPS */
        "CPYWMNSL",
        "PYWCMNDLSX",
        "IQX",
        "",
        "",
        "",
        "",
    ],
    [/* GLO */ "CPABX", "PCABX", "IQX", "", "", "", ""],
    [/* GAL */ "CABXZ", "IQX", "IQX", "ABCXZ", "IQX", "", ""],
    [/* QZS */ "CLSXZ", "LSX", "IQXDPZ", "LSXEZ", "", "", ""],
    [
        /* BDS */
        "IQXDPAN", "IQXDPZ", "DPX", "IQXA", "DPX", "", "",
    ],
    [/* IRN */ "ABCX", "ABCX", "", "", "", "", ""],
];

pub const LEAPS: [[f64; 7]; 19] = [
    [2017.0, 1.0, 1.0, 0.0, 0.0, 0.0, -18.0],
    [2015.0, 7.0, 1.0, 0.0, 0.0, 0.0, -17.0],
    [2012.0, 7.0, 1.0, 0.0, 0.0, 0.0, -16.0],
    [2009.0, 1.0, 1.0, 0.0, 0.0, 0.0, -15.0],
    [2006.0, 1.0, 1.0, 0.0, 0.0, 0.0, -14.0],
    [1999.0, 1.0, 1.0, 0.0, 0.0, 0.0, -13.0],
    [1997.0, 7.0, 1.0, 0.0, 0.0, 0.0, -12.0],
    [1996.0, 1.0, 1.0, 0.0, 0.0, 0.0, -11.0],
    [1994.0, 7.0, 1.0, 0.0, 0.0, 0.0, -10.0],
    [1993.0, 7.0, 1.0, 0.0, 0.0, 0.0, -9.0],
    [1992.0, 7.0, 1.0, 0.0, 0.0, 0.0, -8.0],
    [1991.0, 1.0, 1.0, 0.0, 0.0, 0.0, -7.0],
    [1990.0, 1.0, 1.0, 0.0, 0.0, 0.0, -6.0],
    [1988.0, 1.0, 1.0, 0.0, 0.0, 0.0, -5.0],
    [1985.0, 7.0, 1.0, 0.0, 0.0, 0.0, -4.0],
    [1983.0, 7.0, 1.0, 0.0, 0.0, 0.0, -3.0],
    [1982.0, 7.0, 1.0, 0.0, 0.0, 0.0, -2.0],
    [1981.0, 7.0, 1.0, 0.0, 0.0, 0.0, -1.0],
    [0.0; 7],
];

pub const GLO_FCN: [i32; 32] = [
    1, -4, 5, 6, 1, -4, 5, 6,
    -2, -7, 0, -1, -2, -7, 0, -1,
    4, -3, 3, 2, 4, -3, 3, 2,
    0, 0, 0, 0, 0, 0, 0, 0
];

#[derive(Debug, Clone, Copy, Default)]
pub struct GTime {
    pub time: i64, // 使用 i64 来表示 time_t
    pub sec: f64,
}

impl From<c_var::GTime> for GTime {
    fn from(value: c_var::GTime) -> Self {
        GTime {
            time: value.time,
            sec: value.sec,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SnrMask {
    pub ena: [i32; 2],
    pub mask: [[f64; 9]; NFREQ],
}

#[derive(Debug, Clone, Copy)]
pub struct Pcv {
    pub sat: usize,
    pub r#type: [char; MAXANT],
    pub ts: GTime,
    pub te: GTime,
    pub off: [[f64; 3]; NFREQ],
    pub var: [[f64; 19]; NFREQ],
}

#[derive(Debug, Clone)]
pub struct Pcvs {
    pub n: i32,
    pub nmax: i32,
    pub pcv: Vec<Pcv>, // 使用 Vec 来表示动态数组
}

#[derive(Debug, Clone, Copy)]
pub struct PrcOpt {
    pub mode: i32,
    pub soltype: i32,
    pub nf: i32,
    pub navsys: usize,
    pub elmin: f64,
    pub snrmask: SnrMask,
    pub sateph: i32,
    pub modear: i32,
    pub glomodear: i32,
    pub bdsmodear: i32,
    pub maxout: i32,
    pub minlock: i32,
    pub minfix: i32,
    pub armaxiter: i32,
    pub ionoopt: i32,
    pub tropopt: i32,
    pub dynamics: i32,
    pub tidecorr: i32,
    pub niter: i32,
    pub codesmooth: i32,
    pub intpref: i32,
    pub sbascorr: i32,
    pub sbassatsel: i32,
    pub rovpos: i32,
    pub refpos: i32,
    pub eratio: [f64; NFREQ],
    pub err: [f64; 5],
    pub std: [f64; 3],
    pub prn: [f64; 6],
    pub sclkstab: f64,
    pub thresar: [f64; 8],
    pub elmaskar: f64,
    pub elmaskhold: f64,
    pub thresslip: f64,
    pub maxtdiff: f64,
    pub maxinno: f64,
    pub maxgdop: f64,
    pub baseline: [f64; 2],
    pub ru: [f64; 3],
    pub rb: [f64; 3],
    pub anttype: [[char; MAXANT]; 2],
    pub antdel: [[f64; 3]; 2],
    // pub pcvr: [Pcv; 2],
    pub exsats: [u8; MAXSAT],
    pub maxaveep: i32,
    pub initrst: i32,
    pub outsingle: i32,
    pub rnxopt: [[char; 256]; 2],
    pub posopt: [i32; 6],
    pub syncsol: i32,
    pub odisp: [[f64; 6 * 11]; 2],
    pub freqopt: i32,
    pub pppopt: [char; 256],
}

impl Default for PrcOpt {
    fn default() -> Self {
        PrcOpt {
            mode: PMODE_SINGLE,
            soltype: 2,
            nf: 4,
            navsys: SYS_GPS | SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP,
            elmin: 15.0 * D2R,
            snrmask: SnrMask {
                ena: [0, 0],
                mask: [[0.0; 9]; NFREQ],
            },
            sateph: 0,
            modear: 3,
            glomodear: 0,
            bdsmodear: 1,
            maxout: 50,
            minlock: 0,
            minfix: 50,
            armaxiter: 1,
            ionoopt: 1,
            tropopt: 1,
            dynamics: 1,
            tidecorr: 0,
            niter: 1,
            codesmooth: 0,
            intpref: 1,
            sbascorr: 0,
            sbassatsel: 0,
            rovpos: 0,
            refpos: 1,
            eratio: [150.0, 150.0, 0.0],
            err: [100.0, 0.003, 0.003, 0.0, 1.0],
            std: [30.0, 0.03, 0.3],
            prn: [1E-4, 1E-3, 1E-4, 0.3, 0.7, 0.0],
            sclkstab: 5E-12,
            thresar: [3.0, 0.9999, 0.0, 0.1, 0.05, 0.0, 0.0, 0.0],
            elmaskar: 15.0 * D2R,
            elmaskhold: 15.0 * D2R,
            thresslip: 0.05,
            maxtdiff: 30.0,
            maxinno: 30.0,
            maxgdop: 30.0,
            baseline: [0.0, 0.0],
            ru: [0.0, 0.0, 0.0],
            rb: [0.0, 0.0, 0.0],
            anttype: [['\0'; MAXANT]; 2],
            antdel: [[0.0; 3]; 2],
            exsats: [0; MAXSAT],
            maxaveep: 0,
            initrst: 0,
            outsingle: 0,
            rnxopt: [['\0'; 256]; 2],
            posopt: [0; 6],
            syncsol: 0,
            odisp: [[0.0; 6 * 11]; 2],
            freqopt: 0,
            pppopt: ['\0'; 256],
        }
    }
}

#[derive(Debug, Clone)]
pub struct RnxOpt {
    pub ts: GTime,
    pub te: GTime,
    pub tint: i32,
    pub ttol: f64,
    pub tunit: f64,
    pub rnxver: i32,
    pub navsys: usize,
    pub obstype: i32,
    pub freqtype: i32,
    pub mask: [[char; 64]; NUMSYS],
    pub staid: String,
    pub prog: String,
    pub runby: String,
    pub marker: String,
    pub markerno: String,
    pub markertype: String,
    pub name: [String; 2],
    pub rec: [String; 3],
    pub ant: [String; 3],
    pub apppos: [f64; 3],
    pub antdel: [f64; 3],
    pub glo_cp_bias: [f64; 4],
    pub comment: Box<[String]>,
    pub rcvopt: String,
    pub exsats: [u8; MAXSAT],
    pub glofcn: [i32; 32],
    pub outiono: i32,
    pub outtime: i32,
    pub outleaps: i32,
    pub autopos: i32,
    pub phshift: bool,
    pub halfcyc: i32,
    pub sep_nav: i32,
    pub tstart: GTime,
    pub tend: GTime,
    pub trtcm: GTime,
    pub tobs: [[[char; 4]; MAXOBSTYPE]; NUMSYS],
    pub shift: [[f64; MAXOBSTYPE]; NUMSYS],
    pub nobs: [i32; NUMSYS],
    pub lastouttime: GTime,
    pub iffirstoutput: i32,
    pub ifscan: i32,
}

impl RnxOpt {
    pub fn new() -> Self {
        RnxOpt {
            ts: GTime::default(),
            te: GTime::default(),
            tint: 0,
            ttol: 0.0,
            tunit: 0.0,
            rnxver: 0,
            navsys: 0,
            obstype: 0,
            freqtype: 0,
            mask: [['\0'; 64]; NUMSYS],
            staid: String::new(),
            prog: String::new(),
            runby: String::new(),
            marker: String::new(),
            markerno: String::new(),
            markertype: String::new(),
            name: [String::new(),String::new()],
            rec: [String::new(),String::new(),String::new()],
            ant: [String::new(),String::new(),String::new()],
            apppos: [0.0; 3],
            antdel: [0.0; 3],
            glo_cp_bias: [0.0; 4],
            comment: vec![String::new(); MAXCOMMENT].into_boxed_slice(),
            rcvopt: String::new(),
            exsats: [0; MAXSAT],
            glofcn: [0; 32],
            outiono: 0,
            outtime: 0,
            outleaps: 0,
            autopos: 0,
            phshift: false,
            halfcyc: 0,
            sep_nav: 0,
            tstart: GTime::default(),
            tend: GTime::default(),
            trtcm: GTime::default(),
            tobs: [[['\0'; 4]; MAXOBSTYPE]; NUMSYS],
            shift: [[0.0; MAXOBSTYPE]; NUMSYS],
            nobs: [0; NUMSYS],
            lastouttime: GTime::default(),
            iffirstoutput: 0,
            ifscan: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct SolOpt {
    pub q: i32,
    pub rr: [f64; 6],
    pub qr: [f64; 6],
    pub qv: [f64; 6],
    pub type_: u8,
    pub ns: u8,
    pub time: GTime,
    pub dtr: [f64; 6],
    pub age: f64,
    pub stat: u8,
    pub ratio: f64,
}

impl SolOpt {
    pub fn outecef<W: Write>(&self, writer: &mut W) -> Result<()> {
        let time = time2str(self.time);
        let ecef = format!(
            "{time}{:14.4}{:14.4}{:14.4}{:3}{:3}{:8.4}{:8.4}{:8.4}{:8.4}{:8.4}{:8.4}{:6.2}{:6.1}\n",
            self.rr[0], self.rr[1], self.rr[2],
            self.stat, self.ns,
            self.qr[0].sqrt(), self.qr[1].sqrt(), self.qr[2].sqrt(),
            sqvar(self.qr[3]), sqvar(self.qr[4]), sqvar(self.qr[5]),
            self.age, self.ratio
        );
        writer.write_all(ecef.as_bytes())?;
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Eph {
    pub sat: usize,  // satellite number
    pub iode: i32, // IODE
    pub iodc: i32, // IODC
    pub sva: i32,  // SV accuracy (URA index)
    pub svh: i32,  // SV health (0: ok)
    pub week: i32, // GPS/QZS: gps week, GAL: galileo week
    pub code: i32, // GPS/QZS: code on L2
    // GAL: data source defined as rinex 3.03
    // BDS: data source (0: unknown, 1: B1I, 2: B1Q, 3: B2I, 4: B2Q, 5: B3I, 6: B3Q)
    pub flag: i32, // GPS/QZS: L2 P data flag
    // BDS: nav type (0: unknown, 1: IGSO/MEO, 2: GEO)
    pub toe: GTime, // Toe
    pub toc: GTime, // Toc
    pub ttr: GTime, // T_trans
    // SV orbit parameters
    pub a: f64,
    pub e: f64,
    pub i0: f64,
    pub omg0: f64,
    pub omg: f64,
    pub m0: f64,
    pub deln: f64,
    pub omgd: f64,
    pub idot: f64,
    pub crc: f64,
    pub crs: f64,
    pub cuc: f64,
    pub cus: f64,
    pub cic: f64,
    pub cis: f64,
    pub toes: f64, // Toe (s) in week
    pub fit: f64,  // fit interval (h)
    pub f0: f64,
    pub f1: f64,
    pub f2: f64,       // SV clock parameters (af0, af1, af2)
    pub tgd: [f64; 6], // group delay parameters
    // GPS/QZS: tgd[0] = TGD
    // GAL: tgd[0] = BGD_E1E5a, tgd[1] = BGD_E1E5b
    // CMP: tgd[0] = TGD_B1I, tgd[1] = TGD_B2I/B2b, tgd[2] = TGD_B1Cp
    //      tgd[3] = TGD_B2ap, tgd[4] = ISC_B1Cd, tgd[5] = ISC_B2ad
    pub adot: f64, // Adot for CNAV
    pub ndot: f64, // ndot for CNAV
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Geph {
    pub sat: usize,      // satellite number
    pub iode: i32,     // IODE (0-6 bit of tb field)
    pub frq: i32,      // satellite frequency number
    pub svh: i32,      // satellite health
    pub sva: i32,      // satellite accuracy
    pub age: i32,      // age of operation
    pub toe: GTime,    // epoch of ephemerides (gpst)
    pub tof: GTime,    // message frame time (gpst)
    pub pos: [f64; 3], // satellite position (ecef) (m)
    pub vel: [f64; 3], // satellite velocity (ecef) (m/s)
    pub acc: [f64; 3], // satellite acceleration (ecef) (m/s^2)
    pub taun: f64,     // SV clock bias (s)
    pub gamn: f64,     // relative frequency bias
    pub dtaun: f64,    // delay between L1 and L2 (s)
}

#[derive(Debug, Clone)]
pub struct Nav {
    pub n: usize,
    pub ng: usize,
    pub eph: Vec<Eph>,
    pub geph: Vec<Geph>,
    pub glo_fcn: [i32; 32],
    pub utc_gps: [f64; 8],
    pub utc_glo: [f64; 8],
    pub utc_gal: [f64; 8],
    pub utc_qzs: [f64; 8],
    pub utc_cmp: [f64; 8],
    pub utc_irn: [f64; 9],
    pub ion_gps: [f64; 8],
    pub ion_gal: [f64; 4],
    pub ion_qzs: [f64; 8],
    pub ion_cmp: [f64; 8],
    pub ion_irn: [f64; 8],
    pub cbias: [[f64; 3]; MAXSAT],
}

impl Nav {
    pub fn new() -> Self {
        Nav {
            n: 0,
            ng: 0,
            eph: Vec::new(),
            geph: Vec::new(),
            glo_fcn: [0; 32],
            utc_gps: [0.0; 8],
            utc_glo: [0.0; 8],
            utc_gal: [0.0; 8],
            utc_qzs: [0.0; 8],
            utc_cmp: [0.0; 8],
            utc_irn: [0.0; 9],
            ion_gps: [0.0; 8],
            ion_gal: [0.0; 4],
            ion_qzs: [0.0; 8],
            ion_cmp: [0.0; 8],
            ion_irn: [0.0; 8],
            cbias: [[0.0; 3]; MAXSAT],
        }
    }
}

impl Nav {
    pub fn add_eph(&mut self, eph: Eph) -> bool {
        self.eph.push(eph);
        self.n += 1;
        true
    }
    pub fn add_geph(&mut self, geph: Geph) -> bool {
        self.geph.push(geph);
        self.ng += 1;
        true
    }
}

#[derive(Debug, Clone)]
pub struct Obss {
    pub n: usize,
    pub nmax: usize,
    pub obs: Vec<Obs>,
}

impl Obss {
    pub fn new() -> Self {
        Obss {
            n: 0,
            nmax: 0,
            obs: Vec::new(),
        }
    }

    pub fn add_obs_data(&mut self, data: &Obs) -> bool {
        if data.sat <= 0 {
            return false;
        }
        if self.n >= self.nmax {
            self.nmax = if self.nmax == 0 {
                NINCOBS
            } else {
                self.nmax * 2
            };
        }
        self.obs.push(data.clone());
        self.n += 1;
        return true;
    }
}

#[derive(Debug, Clone, Default)]
pub struct Obs {
    pub time: GTime,
    pub sat: usize,
    pub rcv: u8,
    pub snr: [i32; NFREQ + NEXOBS],
    pub lli: [i32; NFREQ + NEXOBS],
    pub code: [u8; NFREQ + NEXOBS],
    pub l: [f64; NFREQ + NEXOBS],
    pub p: [f64; NFREQ + NEXOBS],
    pub d: [f64; NFREQ + NEXOBS],
}

#[derive(Debug, Clone, Default)]
pub struct Sta {
    pub name: String,          // marker name
    pub marker: String,        // marker number
    pub antdes: String,        // antenna descriptor
    pub antsno: String,        // antenna serial number
    pub rectype: String,       // receiver type descriptor
    pub recver: String,        // receiver firmware version
    pub recsno: String,        // receiver serial number
    pub antsetup: i32,         // antenna setup id
    pub itrf: i32,             // ITRF realization year
    pub deltype: i32,          // antenna delta type (0: enu, 1: xyz)
    pub pos: [f64; 3],         // station position (ecef) (m)
    pub del: [f64; 3],         // antenna position delta (e/n/u or x/y/z) (m)
    pub hgt: f64,              // antenna height (m)
    pub glo_cp_align: i32,     // GLONASS code-phase alignment (0: no, 1: yes)
    pub glo_cp_bias: [f64; 4], // GLONASS code-phase biases {1C, 1P, 2C, 2P} (m)
}

#[derive(Debug, Clone, Default)]
pub struct Stas {
    pub staid: i32,
    pub ts: GTime,
    pub te: GTime,
    pub sta: Sta,
}

#[derive(Debug, Clone)]
pub struct Sigidx {
    pub n: usize,
    pub idx: [i32; MAXOBSTYPE],
    pub pos: [i32; MAXOBSTYPE],
    pub type_: [i32; MAXOBSTYPE],
    pub code: [u8; MAXOBSTYPE],
    pub shift: [f64; MAXOBSTYPE],
    pub pri: [i32; MAXOBSTYPE],
}

impl Default for Sigidx {
    fn default() -> Self {
        Sigidx {
            n: 0,
            idx: [0; MAXOBSTYPE],
            pos: [0; MAXOBSTYPE],
            type_: [0; MAXOBSTYPE],
            code: [0; MAXOBSTYPE],
            shift: [0.0; MAXOBSTYPE],
            pri: [0; MAXOBSTYPE],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Ssat {
    pub sys: i32,
    pub vs: u8,
    pub azel: [f64; 2],
    pub resp: [f64; NFREQ],
    pub resc: [f64; NFREQ],
    pub vsat: [u8; NFREQ],
    pub snr: [i32; NFREQ],
    pub fix: [u8; NFREQ],
    pub slip: [u8; NFREQ + NEXOBS],
    pub half: [u8; NFREQ],
    pub lock: [i32; NFREQ],
    pub outc: [i32; NFREQ],
    pub slipc: [i32; NFREQ],
    pub rejc: [i32; NFREQ],
    pub gf: [f64; NFREQ + NEXOBS],
    pub mw: [f64; NFREQ + NEXOBS],
    pub dl: [f64; NFREQ + NEXOBS],
    pub phw: f64,
    pub pt: [[GTime; NFREQ]; 2],
    pub ph: [[f64; NFREQ]; 2],
    pub dop: [f64; 4],
}

impl Default for Ssat {
    fn default() -> Self {
        Ssat {
            sys: 0,
            vs: 0,
            azel: [0.0; 2],
            resp: [0.0; NFREQ],
            resc: [0.0; NFREQ],
            vsat: [0; NFREQ],
            snr: [0; NFREQ],
            fix: [0; NFREQ],
            slip: [0; NFREQ + NEXOBS],
            half: [0; NFREQ],
            lock: [0; NFREQ],
            outc: [0; NFREQ],
            slipc: [0; NFREQ],
            rejc: [0; NFREQ],
            gf: [0.0; NFREQ + NEXOBS],
            mw: [0.0; NFREQ + NEXOBS],
            dl: [0.0; NFREQ + NEXOBS],
            phw: 0.0,
            pt: [[GTime::default(); NFREQ]; 2],
            ph: [[0.0; NFREQ]; 2],
            dop: [0.0; 4],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Rtk {
    pub sol: SolOpt,
    pub rb: [f32; 6],
    pub nx: i32,
    pub na: i32,
    pub tt: f32,
    pub x: Vec<f32>,
    pub p: Vec<f32>,
    pub xa: Vec<f32>,
    pub pa: Vec<f32>,
    pub nfix: usize,
    pub ssat: Vec<Ssat>,
    pub opt: PrcOpt,
}

impl Rtk {
    pub fn init(opt: &PrcOpt) -> Self {
        Self {
            sol: SolOpt::default(),
            rb: [0.0; 6],
            nx: 9,
            na: 9,
            tt: 0.0,
            x: vec![0.0; 9],
            p: vec![0.0; 9 * 9],
            xa: vec![0.0; 9],
            pa: vec![0.0; 9 * 9],
            nfix: 0,
            ssat: vec![Ssat::default(); MAXSAT],
            opt: opt.clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RtcmRes {
    pub r#type: i32,
    pub staid: i32,
    pub stah: i32,
    pub seqno: i32,
    pub time: GTime,
    pub time_s: GTime,
    pub obss: Obss,
    pub eph: Eph,
    pub geph: Geph,
    pub obsflag: bool,
    pub ephsat: usize,
    pub len: usize,
    pub sta: Sta,
}

impl RtcmRes {
    pub fn new() -> Self {
        Self {
            r#type: 0,
            staid: 0,
            stah: 0,
            seqno: 0,
            time: GTime::default(),
            time_s: GTime::default(),
            obss: Obss::new(),
            eph: Eph::default(),
            geph: Geph::default(),
            obsflag: false,
            ephsat: 0,
            len: 0,
            sta: Sta::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct QcOpt {
    pub ts: GTime,
    pub te: GTime,
    pub tr: GTime,
    pub ti: f64,
    pub detail: i32,
    pub popt: PrcOpt,
    pub rnxver:i32,
    pub navsys:usize,
}