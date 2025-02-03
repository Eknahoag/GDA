use crate::basic::var::*;
use std::collections::VecDeque;

#[allow(dead_code)]
const QUEUE_SIZE: usize = 50;
#[allow(dead_code)]
const THRES_MW_JUMP: f32 = 10.0;
#[allow(dead_code)]
const LIMIT_EL: f32 = 10.0;

#[derive(Debug, Clone)]
pub struct Mp {
    pub satlist: Vec<Vec<VecDeque<f32>>>,
    pub sat_std: [[Vec<f32>; NFREOBS]; MAXSAT],
    pub satnum: [[usize; NFREOBS]; MAXSAT],
}

impl Mp {
    pub fn new() -> Self {
        let empty_vec = Vec::new();
        let empty_vec_array: [Vec<f32>; NFREOBS] = array_init::array_init(|_| empty_vec.clone());
        let data: [[Vec<f32>; NFREOBS]; MAXSAT] = array_init::array_init(|_| empty_vec_array.clone());
        Mp {
            satlist: (0..MAXSAT).map(|_| (0..NFREOBS).map(|_| VecDeque::new()).collect()).collect(),
            sat_std: data,
            satnum: [[0; NFREOBS]; MAXSAT],
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ExpandHav {
    pub totsatobs: i32,
    pub totsatobs_el: i32,
    pub totobs: i32,
    pub totobs_el: i32,
}

#[derive(Debug, Clone, Default)]
pub struct Noise{
    pub single_diff:[f64;NFREOBS],
    pub double_diff:[f64;NFREOBS],
    pub triple_diff:[Vec<f64>;NFREOBS],
    pub form:[f64;NFREOBS],
    pub form_single_diff:[f64;NFREOBS],
    pub form_double_diff:[f64;NFREOBS]
}

#[derive(Debug, Clone, Default)]
pub struct TotInfo {
    pub firstepoch: GTime,
    pub lastepoch: GTime,
    pub sample: i32,
    pub minele: f32,
    pub ep_rt: f32,
    pub sat_rt: f32,
    pub sat_el_rt: f32,
    pub obs_rt: f32,
    pub obs_el_rt: f32,
    pub oslps: i32,
    pub enu_std: [f64; 3],
    pub mark: f32,
    pub exp: ExpandHav,
    pub hav: ExpandHav,
    pub totobs: i32,
    pub totslip: i32,
}

#[derive(Debug, Clone, Default, Copy)]
pub struct SysInfo {
    pub exp_obs: i32,
    pub hav_obs: i32,
    pub satnum: i32,
    pub epnum: i32,
    pub meansat: f32,
    pub codenum: i32,
    pub x_co_sv: i32,
    pub x_ph_sv: i32,
    pub slipnum: i32,
    pub oslps: i32,
    pub totmp: [f32; NFREQ + NEXOBS],
    pub code: [u8; NFREQ + NEXOBS],
}

#[derive(Debug, Clone)]
pub struct SatInfo {
    pub exp_obs: i32,
    pub exp_obs_el: i32,
    pub hav_obs: i32,
    pub hav_obs_el: i32,
    pub sat_rt: f32,
    pub sat_el_rt: f32,
    pub sat_mp: [f32; NFREQ + NEXOBS],
    pub azel: [Vec<f32>; 2],
    pub resp: [f32; 2],
    pub snr: [f32; NFREQ + NEXOBS],
    pub eptime: Vec<GTime>,
    pub aetime: Vec<GTime>,
    pub slip_flag: Vec<i32>,
    pub obs_flag: Vec<i32>,
    pub l1snr: Vec<f32>,
    pub freq: [i32; NFREQ + NEXOBS],
    pub code: [u8; NFREQ + NEXOBS],
    pub verify: bool,
    pub pretime: GTime,
    pub ngap: i32,
    pub slip: [u8; NFREQ + NEXOBS],
    pub gf: [f64; NFREQ + NEXOBS],
    pub mw: [f64; NFREQ + NEXOBS],
    pub slipnum: i32,
    pub obsnum: i32,
    pub oslps: i32,
    pub pnoise:Noise,
    pub lnoise:Noise
}

impl SatInfo {
    pub fn new() -> Self {
        SatInfo {
            exp_obs: 0,
            exp_obs_el: 0,
            hav_obs: 0,
            hav_obs_el: 0,
            sat_rt: 0.0,
            sat_el_rt: 0.0,
            sat_mp: [0.0; NFREQ + NEXOBS],
            azel: [Vec::new(), Vec::new()],
            resp: [0.0; 2],
            snr: [0.0; NFREQ + NEXOBS],
            eptime: Vec::new(),
            aetime: Vec::new(),
            slip_flag: Vec::new(),
            obs_flag: Vec::new(),
            l1snr: Vec::new(),
            freq: [0; NFREQ + NEXOBS],
            code: [0; NFREQ + NEXOBS],
            verify: false,
            pretime: GTime::default(),
            ngap: 0,
            slip: [0; NFREQ + NEXOBS],
            gf: [0.0; NFREQ + NEXOBS],
            mw: [0.0; NFREQ + NEXOBS],
            slipnum: 0,
            obsnum: 0,
            oslps: 0,
            pnoise:Noise::default(),
            lnoise:Noise::default()
        }
    }
}

#[derive(Debug, Clone, Default, Copy)]
pub struct SigInfo {
    pub exp_obs: i32,
    pub exp_obs_el: i32,
    pub hav_obs: i32,
    pub hav_obs_el: i32,
    pub snr: f32,
    pub snrnum: usize,
    pub pnoise:f64,
    pub lnoise:f64,
    pub pnoise_box:[f64;5],
    pub lnoise_box:[f64;5],
}

#[derive(Debug, Clone)]
pub struct Pos {
    pub coor: [Vec<f64>; 3],
    pub dop:[Vec<f64>;4],
    pub time:Vec<GTime>,
    pub ave_pos:[f64;3]
}

impl Pos {
    pub fn new() -> Self {
        Pos {
            coor: [Vec::new(), Vec::new(), Vec::new()],
            dop: [Vec::new(), Vec::new(), Vec::new(),Vec::new()],
            time:Vec::new(),
            ave_pos:[0.0;3]
        }
    }
}

#[derive(Debug, Clone)]
pub struct AutoTi {
    pub pre: GTime,
    pub tiarr: Vec<i32>,
    pub ti: i32,
}

impl AutoTi {
    pub fn new() -> Self {
        AutoTi {
            pre: GTime::default(),
            tiarr: Vec::new(),
            ti: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct QC {
    pub totinfo: TotInfo,
    pub epnum: usize,
    pub sysinfo: Vec<SysInfo>,
    pub satinfo: Vec<SatInfo>,
    pub siginfo: Vec<Vec<SigInfo>>,
    pub mp: Box<Mp>,
    pub ati: AutoTi,
    pub pos: Box<Pos>,
    pub satlist: Vec<usize>,
}

impl QC {
    pub fn init() -> Self {
        Self {
            totinfo: TotInfo::default(),
            epnum: 0,
            sysinfo: vec![SysInfo::default(); 6],
            satinfo: vec![SatInfo::new(); MAXSAT],
            siginfo: vec![vec![SigInfo::default(); NFREOBS]; NUMSYS],
            mp: Box::new(Mp::new()),
            ati: AutoTi::new(),
            pos: Box::new(Pos::new()),
            satlist: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum Outtype {
    File(String),
    Redis(RedisParams),
    All(String, RedisParams),
}

#[derive(Debug, Clone)]
pub struct RedisParams {
    pub ip: Option<String>,
    pub port: Option<String>,
    pub auth: Option<String>,
}


