use std::cell::RefCell;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, BufWriter, Error, Write,self};
use std::panic;
use std::path::Path;
use array_init::array_init;
use indicatif::{ProgressBar, ProgressStyle};
use nalgebra::{min, Vector3};
use num_traits::pow;
use crate::basic::code::code2idx;
use crate::basic::code::{code2obs, sat2freq};
use crate::basic::eph::{ephclk, satpos};
use crate::basic::func::{cal_box, rms, sqr, Sqrt};
use crate::basic::pos::{ecef2enu, ecef2pos, geodist, inputobs, pntpos, satazel};
use crate::basic::sat::{satno2id, satsys};
use crate::basic::time::{time2str, timeadd, timediff, timeget};
use crate::basic::var::*;
use crate::convbin::decode::{decode_epoch, DecodeError, Decoder};
use crate::convbin::decode::DecodeError::{DataError, EndofFile, IoError};
use crate::basic::read::{readobsnav, uniqnav};
use crate::qc::cors::qc2redis;
use crate::qc::qc_var::*;

const THRES_MW_JUMP: f64 = 10.0;
const QUEUE_SIZE: usize = 50;
pub const LIMIT_EL: f64 = 10.0;

// Use thread-localized storage for `data`.
thread_local! {
    static PANIC_DATA: RefCell<Option<Vec<u8>>> = RefCell::new(None);
}

fn setup_panic_hook() {
    panic::set_hook(Box::new(|info| {
        // Getting panic information
        let panic_message = info.payload().downcast_ref::<&str>().unwrap_or(&"Unknown panic message");

        // Print panic message
        println!("Panic occurred: {}", panic_message);

        // Printing data during panic
        PANIC_DATA.with(|data| {
            let data = data.borrow();
            for d in data.iter().enumerate() {
                println!("Data causing panic: {:?}", d);
            }
        });
    }));
}

// ----------------------------
// Data Processing Related Functions
// ----------------------------

// calculate the queue mean
fn que_mean(que: &VecDeque<f32>) -> f32 {
    let sum: f32 = que.iter().sum();
    sum / que.len() as f32
}

// calculate the queue std
fn que_std(que: &VecDeque<f32>) -> f32 {
    let mean = que_mean(que);
    let variance: f32 = que.iter().map(|&x| (x - mean).powi(2)).sum::<f32>() / que.len() as f32;
    variance.sqrt()
}

// determine if outlier of new value in queue
fn que_outlier(que: &mut VecDeque<f32>) -> bool {
    if que_std(&que) > 5.0 { return true; }
    return false;
}

// get mp std
fn que_make(satidx: usize, j: usize, mp_value: f32, mp: &mut Mp) {
    let mean_pre = que_mean(&mp.satlist[satidx][j]);
    let std_pre = que_std(&mp.satlist[satidx][j]);
    mp.satlist[satidx][j].push_back(mp_value);
    let std_enq = que_std(&mp.satlist[satidx][j]);

    if !mean_pre.is_nan() && que_outlier(&mut mp.satlist[satidx][j]) {
        if std_pre != 0.0 {
            mp.sat_std[satidx][j].push(std_pre);
        }
        mp.satlist[satidx][j] = VecDeque::new();
        mp.satlist[satidx][j].push_back(mp_value);
    } else {
        if mp.satlist[satidx][j].len() == QUEUE_SIZE {
            mp.satlist[satidx][j].pop_front();
            mp.sat_std[satidx][j].push(std_enq);
        }
    }
}

// Find mode of array
fn find_mode(arr: &[i32]) -> i32 {
    let mut mode = arr[0];
    let mut max_count = 1;

    for &num in arr.iter() {
        if num == 0 {
            continue;
        }
        let count = arr.iter().filter(|&&x| x == num).count();
        if count > max_count {
            mode = num;
            max_count = count;
        }
    }

    mode
}

// Initialize epoch
fn epoch_init(satinfo: &mut [SatInfo], obs: &Vec<Obs>, n: usize) {
    for i in 0..n {
        let sat = obs[i].sat - 1;

        // Verify if obs exists
        for j in 0..NFREOBS {
            if obs[i].p[j] != 0.0 || obs[i].l[j] != 0.0 {
                satinfo[sat].verify = true;
                break;
            }
        }

        // Initialize slip
        for j in 0..(NFREQ + NEXOBS) {
            satinfo[sat].slip[j] = 0;
        }

        // Set ngap, reset gf and mw if satobs interrupt for a period
        if satinfo[sat].pretime.time != 0 && obs[i].time.time - satinfo[sat].pretime.time >= 600 {
            satinfo[sat].ngap += 1; // ngap count
            for j in 1..(NFREQ + NEXOBS) {
                satinfo[sat].gf[j] = 0.0;
                satinfo[sat].mw[j] = 0.0;
            }
        }

        // Update pretime
        satinfo[sat].pretime = obs[i].time;
    }
}

// Initialize epoch for real-time
fn epoch_init_rt(
    satinfo: &mut [SatInfo],
    sysinfo: &mut Vec<SysInfo>,
    siginfo: &mut Vec<Vec<SigInfo>>,
    mp: &mut Mp,
    obs: &Vec<Obs>,
    n: usize) {
    for i in 0..n {
        let satidx = obs[i].sat - 1;

        // Verify if obs exists
        for j in 0..NFREOBS {
            if obs[i].p[j] != 0.0 || obs[i].l[j] != 0.0 {
                satinfo[satidx].verify = true;
                break;
            }
        }

        // Initialize slip
        for j in 0..(NFREQ + NEXOBS) {
            satinfo[satidx].slip[j] = 0;
        }

        // Set ngap, reset gf and mw if satobs interrupt for a period
        if satinfo[satidx].pretime.time != 0 && obs[i].time.time - satinfo[satidx].pretime.time >= 600 {
            satinfo[satidx].ngap += 1; // ngap count
            for j in 1..(NFREQ + NEXOBS) {
                satinfo[satidx].gf[j] = 0.0;
                satinfo[satidx].mw[j] = 0.0;
            }
        }

        // Update pretime
        satinfo[satidx].pretime = obs[i].time;

        // Init obs count
        satinfo[satidx].hav_obs = 0;
        satinfo[satidx].hav_obs_el = 0;
        satinfo[satidx].exp_obs = 0;
        satinfo[satidx].exp_obs_el = 0;
        satinfo[satidx].obsnum = 0;
        satinfo[satidx].slipnum = 0;
        // Init azel
        satinfo[satidx].azel = [Vec::new(), Vec::new()];
    }
    for i in 0..NUMSYS {
        for j in 0..NFREOBS {
            siginfo[i][j].exp_obs = 0;
            siginfo[i][j].exp_obs_el = 0;
            siginfo[i][j].hav_obs = 0;
            siginfo[i][j].hav_obs_el = 0;
            // Init mp
            sysinfo[i].totmp[j] = 0.0;
            satinfo[i].sat_mp[j] = 0.0;
            mp.sat_std[i][j] = Vec::new();
        }
    }
}

// time check
pub fn set_qctime(qc: &mut QC, obstime: GTime) -> bool {
    if qc.totinfo.firstepoch.time == 0 {
        qc.totinfo.firstepoch = obstime;
    }
    if obstime.time < qc.totinfo.firstepoch.time {
        return false;
    }
    if obstime.time < qc.totinfo.lastepoch.time {
        return false;
    }
    if qc.totinfo.lastepoch.time != obstime.time {
        qc.epnum += 1;
    }
    qc.totinfo.lastepoch = obstime;

    true
}

// scan rtcm file
fn scan_rtcm(
    decoder: &Decoder,
    rtcmfile: &str,
    navs: &mut Nav,
    qc: &mut QC,
    navsys: usize,
    tr: GTime,
    loop_num: &mut u64,
) -> Result<(), DecodeError> {
    let path = Path::new(rtcmfile);
    let file = File::open(path).map_err(|_| IoError)?;
    let mut reader = BufReader::new(file);

    let mut trt = tr;

    let mut count: u64 = 0;
    let mut obsnum = 0;
    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap()
        .tick_strings(&["-", "\\", "|", "/"]));

    loop {
        let mut obss = Obss::new();
        match decode_epoch(decoder, &mut reader, &mut obss, Some(navs), trt) {
            Err(e) => match e {
                IoError => break,
                EndofFile => break,
                _ => continue
            }
            Ok(true) => {
                for i in 0..obss.n {
                    let obs = obss.obs[i].clone();
                    let sys = satsys(obs.sat);
                    let syslog = sys.ilog2() as usize;
                    if (sys & navsys) == 0 {
                        continue;
                    }

                    // restore satlist
                    let mut satflag = true;
                    for j in 0..qc.satlist.len() {
                        if qc.satlist[j] == obs.sat {
                            satflag = false;
                            break;
                        }
                    }
                    if satflag {
                        qc.satlist.push(obs.sat)
                    }

                    // restore code information
                    for j in 0..NFREOBS {
                        if obs.code[j] > 0 && qc.sysinfo[syslog].code[j] == 0 {
                            qc.sysinfo[syslog].code[j] = obs.code[j];
                            qc.sysinfo[syslog].codenum += 1;
                        }

                        if qc.satinfo[obs.sat - 1].code[j] == 0 &&
                            (obs.l[j] != 0.0 || obs.p[j] != 0.0) {
                            qc.satinfo[obs.sat - 1].code[j] = obs.code[j];
                        }
                    }
                }

                // calculate sample
                cal_sample(&mut qc.ati, obss.obs[0].time);

                trt = obss.obs[0].time;

                count += 1;
                obsnum += obss.obs.len();
                let message = format!("Scanning RTCM:  Epoch:{},Obs:{},Nav:{}",
                                      count, obsnum, navs.eph.len() + navs.geph.len());
                pb.set_message(message.to_string());
                pb.tick();
            }
            Ok(false) => continue
        }
    }

    // check nav
    uniqnav(navs);

    let message = format!("Finish scanning, totally {} epochs, {} sat obs and {} sat nav",
                          count, obsnum, navs.eph.len() + navs.geph.len());
    pb.finish_with_message(message.to_string());

    if qc.satlist.len() == 0 || qc.ati.ti <= 0 {
        return Err(DataError);
    }

    *loop_num = count;

    Ok(())
}


// ----------------------------
// Signal Processing Related Functions
// ----------------------------

// Extracts the signal frequency band based on the system and observation code.
fn getfreqband(sys: usize, code: u8) -> &'static str {
    let obs_code = code2obs(code);

    match sys {
        SYS_GPS => match obs_code.chars().next() {
            Some('1') => "L1",
            Some('2') => "L2",
            Some('5') => "L5",
            _ => "",
        },
        SYS_GLO => match obs_code.chars().next() {
            Some('1') => "G1",
            Some('4') => "G1a",
            Some('2') => "G2",
            Some('6') => "G2a",
            Some('3') => "G3",
            _ => "",
        },
        SYS_GAL => match obs_code.chars().next() {
            Some('1') => "E1",
            Some('5') => "E5a",
            Some('7') => "E5b",
            Some('8') => "E5",
            Some('6') => "E6",
            _ => "",
        },
        SYS_CMP => match obs_code {
            "1D" | "1P" | "1X" => "B1C",
            "1S" | "1L" | "1Z" => "B1A",
            "2I" | "2Q" | "2X" => "B1",
            "5D" | "5P" | "5X" => "B2a",
            "7I" | "7Q" | "7X" => "B2",
            "7D" | "7P" | "7Z" => "B2b",
            "8L" | "8Q" | "8X" => "B2a+B2b",
            "6I" | "6Q" | "6X" => "B3",
            "6D" | "6P" | "6Z" => "B3A",
            _ => "",
        },
        SYS_QZS => match obs_code.chars().next() {
            Some('1') => "L1",
            Some('2') => "L2",
            Some('5') => "L5",
            Some('6') => "L6",
            _ => "",
        },
        _ => "",
    }
}

// Extracts the signal frequency field from a given code
pub fn codefreq(code: u8) -> i32 {
    let cf = code2obs(code);
    if cf == "" {
        return 0;
    }
    cf.chars().next().unwrap().to_digit(10).unwrap() as i32
}

// check bds signal type
fn bds_sig_check(code: usize) -> i32 {
    if code <= 0 || code > MAXCODE as usize {
        return 0;
    }

    if OBS_CODES[code].chars().next().and_then(|c| c.to_digit(10)).map(|d| d as i32).unwrap() == 2 {
        return 1;
    }
    if OBS_CODES[code].chars().next().and_then(|c| c.to_digit(10)).map(|d| d as i32).unwrap() == 1 {
        return 3;
    }
    if OBS_CODES[code].chars().next().and_then(|c| c.to_digit(10)).map(|d| d as i32).unwrap() == 5 {
        return 3;
    }
    if OBS_CODES[code].chars().next().and_then(|c| c.to_digit(10)).map(|d| d as i32).unwrap() == 8 {
        return 3;
    }
    if ["7I", "7Q", "7X"].contains(&OBS_CODES[code]) {
        return 2;
    }
    if ["7D", "7P", "7Z"].contains(&OBS_CODES[code]) {
        return 3;
    }
    if ["6I", "6Q", "6X"].contains(&OBS_CODES[code]) {
        return 1;
    }
    if ["6D", "6P", "6Z"].contains(&OBS_CODES[code]) {
        return 3;
    }
    0
}

// check if gps signal is l5
fn gps_l5_check(sat: usize) -> bool {
    const GPS_NO_L5: [usize; 14] = [5, 7, 12, 15, 17, 29, 31, 2, 13, 16, 19, 20, 21, 22];

    // 假设 satsys 函数已经定义，返回系统类型
    let sys = satsys(sat);

    // 假设 SYS_GPS 常量已经定义
    if sys != SYS_GPS {
        return false;
    }

    for &gps_no in GPS_NO_L5.iter() {
        if gps_no == sat {
            return false;
        }
    }

    true
}

// ----------------------------
// Data calculation/statistics related functions
// ----------------------------

// Melbourne-Wubbena linear combination
fn mwmeas(obs: &Obs, nav: &Nav, f2: usize) -> f64 {
    let freq1 = sat2freq(obs.sat, obs.code[0], Some(nav));
    let freq2 = sat2freq(obs.sat, obs.code[f2], Some(nav));

    if freq1 == 0.0
        || freq2 == 0.0
        || obs.l[0] == 0.0
        || obs.l[f2] == 0.0
        || obs.p[0] == 0.0
        || obs.p[f2] == 0.0
    {
        return 0.0;
    }

    (obs.l[0] - obs.l[f2]) * CLIGHT / (freq1 - freq2)
        - (freq1 * obs.p[0] + freq2 * obs.p[f2]) / (freq1 + freq2)
}

// geometry-free phase measurement
fn gfmeas(obs: &Obs, nav: &Nav, f2: usize) -> f64 {
    let freq1 = sat2freq(obs.sat, obs.code[0], Some(nav));
    let freq2 = sat2freq(obs.sat, obs.code[f2], Some(nav));

    if freq1 == 0.0 || freq2 == 0.0 || obs.l[0] == 0.0 || obs.l[f2] == 0.0 {
        return 0.0;
    }

    (obs.l[0] / freq1 - obs.l[f2] / freq2) * CLIGHT
}

#[allow(dead_code)]
// detect cycle slip by LLI
fn detslp_ll(rtk: &mut Rtk, obs: &[Obs], n: usize) {
    for i in 0..n.min(MAXOBS) {
        for j in 0..(NFREQ + NEXOBS) {
            let obs_i = &obs[i];
            let ssat = &mut rtk.ssat[obs_i.sat - 1];

            if obs_i.l[j] == 0.0
                || ssat.pt[0][j].time == 0 // 假设 Time {} 是无效时间
                || (obs_i.lli[j] & 3) == 0
                || timediff(obs_i.time, ssat.pt[0][j]) < DTTOL
            {
                continue;
            }
        }
    }
}

// detect cycle slip by geometry free phase jump
pub fn detslp_gf(satinfo: &mut [SatInfo], obs: &[Obs], n: usize, nav: &Nav) {
    for i in 0..n.min(MAXOBS) {
        let mut f2 = 0;
        for k in 1..NFREQ + NEXOBS {
            if obs[i].code[k] != 0 {
                f2 = k;

                let g1 = gfmeas(&obs[i], nav, f2);
                if g1 == 0.0 {
                    continue;
                }

                let freq1 = sat2freq(obs[i].sat, obs[i].code[0], Some(nav));
                let freq2 = sat2freq(obs[i].sat, obs[i].code[f2], Some(nav));

                let g0 = satinfo[obs[i].sat - 1].gf[f2];
                satinfo[obs[i].sat - 1].gf[f2] = g1;

                let thresslip = (CLIGHT / freq1 + CLIGHT / freq2) * 0.5 * 0.5;
                if g0 != 0.0 && (g1 - g0).abs() > thresslip {
                    satinfo[obs[i].sat - 1].slip[0] |= 1;
                    satinfo[obs[i].sat - 1].slip[f2] |= 1;
                }
            }
        }
        if f2 == 0 {
            for j in 1..NFREQ + NEXOBS {
                satinfo[obs[i].sat - 1].gf[j] = 0.0;
            }
        }
    }
}

// detect slip by Melbourne-Wubbena linear combination jump
pub fn detslp_mw(satinfo: &mut [SatInfo], obs: &[Obs], n: usize, nav: &Nav) {
    for i in 0..n.min(MAXOBS) {
        let mut f2 = 0;
        for k in 1..NFREQ + NEXOBS {
            if obs[i].code[k] != 0 {
                f2 = k;

                let w1 = mwmeas(&obs[i], nav, f2);
                if w1 == 0.0 {
                    satinfo[obs[i].sat - 1].mw[f2] = 0.0;
                    continue;
                }

                let w0 = satinfo[obs[i].sat - 1].mw[f2];
                satinfo[obs[i].sat - 1].mw[f2] = w1;

                if w0 != 0.0 && (w1 - w0).abs() > THRES_MW_JUMP {
                    satinfo[obs[i].sat - 1].slip[0] |= 1;
                    satinfo[obs[i].sat - 1].slip[f2] |= 1;
                }
            }
        }
        if f2 == 0 {
            for j in 1..NFREQ + NEXOBS {
                satinfo[obs[i].sat - 1].mw[j] = 0.0;
            }
        }
    }
}

// Cal mp by epoch
pub fn get_multipath(mp: &mut Mp, satinfo: &[SatInfo], nav: &Nav, obs: &Vec<Obs>, n: usize) {
    for i in 0..n {
        let sat = obs[i].sat;
        let satidx = sat - 1;
        let sys = satsys(sat);
        let mut mp_values = [0.0; NFREQ + NEXOBS];
        let mut refresh_flag = false;
        let mut cancal = [0; NFREQ + NEXOBS];

        for j in 0..NFREQ + NEXOBS {
            if satinfo[satidx].slip[j] != 0 {
                refresh_flag = true;
                break;
            }
        }
        /* init mp queue when slip */
        if refresh_flag {
            for j in 0..NFREQ + NEXOBS {
                if mp.satlist[satidx][j].len() > 1 && mp.satlist[satidx][j].len() < 50 {
                    let std = que_std(&mp.satlist[satidx][j]);
                    mp.sat_std[satidx][j].push(std);
                }
                mp.satlist[satidx][j] = VecDeque::new();
            }
        }

        if obs[i].p[0] == 0.0 || obs[i].l[0] == 0.0 {
            continue;
        }
        /* count cancal num */
        for j in 0..NFREQ + NEXOBS {
            if obs[i].p[j] != 0.0 && obs[i].l[j] != 0.0 {
                cancal[j] = 1;
            }
        }

        if cancal.iter().filter(|&&x| x == 1).count() < 2 {
            continue;
        }

        let mut f = [-1; NFREQ + NEXOBS];
        for j in 0..NFREQ + NEXOBS {
            if obs[i].code[j] != 0 {
                f[j] = code2idx(sys, obs[i].code[j]);
            }
            let mut k = j as i32 - 1;
            while k >= 0 {
                if f[j] == f[k as usize] {
                    f[j] = -1; /* Consider case like where C2W is absent but C2X is present */
                }
                k -= 1;
            }
        }

        let mut freq = [0.0; NFREQ + NEXOBS];
        for j in 0..NFREQ + NEXOBS {
            if f[j] >= 0 {
                freq[j] = sat2freq(sat, obs[i].code[j], Some(nav));
            }
        }

        for j in 1..NFREQ + NEXOBS {
            if f[j] == 1 {
                if cancal[j] == 0 {
                    break;
                }
                let freq1 = freq[0];
                let freq2 = freq[j];
                let c = sqr(freq1 / freq2);
                let a = 2.0 / (c - 1.0);
                let b = 2.0 * c / (c - 1.0);
                if mp_values[0] == 0.0 {
                    mp_values[0] = obs[i].p[0] - (1.0 + a) * obs[i].l[0] * CLIGHT / freq1
                        + a * obs[i].l[j] * CLIGHT / freq2;
                    que_make(satidx, 0, mp_values[0] as f32, mp);
                }
                mp_values[j] = obs[i].p[j] - b * obs[i].l[0] * CLIGHT / freq1
                    + (b - 1.0) * obs[i].l[j] * CLIGHT / freq2;
                que_make(satidx, j, mp_values[j] as f32, mp);
            }
        }

        for j in 1..NFREQ + NEXOBS {
            if f[j] > 1 {
                if cancal[j] == 0 {
                    continue;
                }
                let freq1 = freq[0];
                let freq2 = freq[j];
                let c = sqr(freq1 / freq2);
                let a = 2.0 / (c - 1.0);
                let b = 2.0 * c / (c - 1.0);
                mp_values[j] = obs[i].p[j] - b * obs[i].l[0] * CLIGHT / freq1
                    + (b - 1.0) * obs[i].l[j] * CLIGHT / freq2;
                que_make(satidx, j, mp_values[j] as f32, mp);
                if mp_values[0] == 0.0 {
                    mp_values[0] = obs[i].p[0] - (1.0 + a) * obs[i].l[0] * CLIGHT / freq1
                        + a * obs[i].l[j] * CLIGHT / freq2;
                    que_make(satidx, 0, mp_values[0] as f32, mp);
                }
            }
        }
    }
}

// Count each signal observation
fn sigcount(
    siginfo: &mut Vec<Vec<SigInfo>>,
    obs: &Obs,
    code: &[u8; NFREQ + NEXOBS],
    ifel: bool,
) {
    let mut hav = [0; NFREQ + NEXOBS];
    let mut exp = [0; NFREQ + NEXOBS];

    let sys = satsys(obs.sat);
    let syslog = (sys as f64).log2() as usize;

    for i in 0..NFREQ + NEXOBS {
        if code[i] != 0 {
            let mut found = false;
            for j in 0..NFREQ + NEXOBS {
                if code[i] == obs.code[j] {
                    /* find current obs code */
                    if obs.p[j] != 0.0 || obs.l[j] != 0.0 {
                        hav[i] += 1;
                    }
                    exp[i] += 1;
                    found = true;
                    break;
                }
            }
            if !found {
                /* The current code is missing */
                if (code2idx(SYS_GPS, code[i]) == 2 && gps_l5_check(obs.sat))
                    || (bds_sig_check(code[i] as usize) == 3 && obs.sat > 121)
                    || (bds_sig_check(code[i] as usize) == 2 && obs.sat <= 121)
                {
                    exp[i] += 1;
                }
            }
        }
        if ifel {
            siginfo[syslog][i].hav_obs_el += hav[i];
            siginfo[syslog][i].exp_obs_el += exp[i];
        } else {
            siginfo[syslog][i].hav_obs += hav[i];
            siginfo[syslog][i].exp_obs += exp[i];
        }
    }
}

// Count obsnum, slipnum and  starting point to end point of ts
fn osfcount(satinfo: &mut Vec<SatInfo>, obs: &Vec<Obs>, n: usize, ti: i64) {
    for i in 0..n {
        let satidx = obs[i].sat - 1;
        let epnum = satinfo[satidx].eptime.len();
        for j in 0..NFREQ + NEXOBS {
            /* slip occur */
            if satinfo[satidx].slip[j] != 0 {
                satinfo[satidx].slipnum += 1;
                satinfo[satidx].slip_flag.push(1);
            } else { satinfo[satidx].slip_flag.push(0) }
            /* Count obsnum */
            if obs[i].p[j] != 0.0 {
                satinfo[satidx].obsnum += 1;
            }
        }
        if satinfo[satidx].verify {
            let time_d;
            let pen_flag;
            if epnum > 2 {
                time_d =
                    satinfo[satidx].eptime[epnum - 1].time - satinfo[satidx].eptime[epnum - 2].time;
                pen_flag = satinfo[satidx].obs_flag[epnum - 2];
            } else {
                time_d = 0;
                pen_flag = 0;
            }
            if epnum == 1 {
                /* first epoch */
                satinfo[satidx].obs_flag.push(1);
            } else if epnum == 2 {
                satinfo[satidx].obs_flag.push(2);
            } else if (time_d <= ti) && pen_flag == 1 {
                /* previous obs_flag==1 and epoch is continuous */
                satinfo[satidx].obs_flag.push(2);
            } else if (time_d <= ti) && pen_flag == 2 {
                /* previous obs_flag==2 and epoch is continuous */
                satinfo[satidx].obs_flag.push(2);
                satinfo[satidx].obs_flag[epnum - 2] = 0;
            } else if time_d > ti {
                /* epoch not continuous */
                satinfo[satidx].obs_flag.push(1);
            }
        }
    }
}

// Cal auto sample(ti)
pub fn cal_sample(ati: &mut AutoTi, ref_time: GTime) {
    if ati.tiarr.len() < 60 && ref_time.time != 0 {
        if ati.pre.time != 0 { ati.tiarr.push(timediff(ref_time, ati.pre) as i32) }
        if ati.tiarr.len() > 0 {
            ati.ti = find_mode(&ati.tiarr);
        }
        ati.pre = ref_time;
    }
}

// Cal the enu of each spp result
fn get_pos(pos: &mut Pos, rtk: &Rtk) {
    if rtk.sol.rr[0] != 0.0 {
        let mut pass_flag = false;
        for i in 0..3 {
            if !pos.coor[i].is_empty() &&
                (pos.coor[i].last().unwrap() - rtk.sol.rr[i]).abs() > 50.0 { pass_flag = true }
        }
        if !pass_flag {
            for i in 0..3 { pos.coor[i].push(rtk.sol.rr[i]) }

            let t = if rtk.sol.time.sec > 0.0 {
                GTime {
                    time: rtk.sol.time.time + 1,
                    sec: 0.0,
                }
            } else { rtk.sol.time };
            // Storage time and dops
            pos.time.push(t);
            for j in 0..MAXSAT {
                if rtk.ssat[j].dop[1] > 0.0
                {
                    for k in 0..4 {
                        pos.dop[k].push(rtk.ssat[j].dop[k])
                    }
                    break;
                }
            }
            // Ensure consistent length
            if pos.time.len()>0 && pos.time.len()!=pos.dop[0].len(){
                for k in 0..4 {
                    pos.dop[k].push(0.0)
                }
            }
            // Calculate the mean value of the coordinates
            for i in 0..3 {
                pos.ave_pos[i] = (pos.ave_pos[i] * (pos.time.len() - 1) as f64 + rtk.sol.rr[i]) / pos.time.len() as f64
            }
        }
    }
}

// Calculate the enu of spp result
fn cal_sppneu(pos: &mut Pos, enu: &mut [f64]) {
    // Convert xyz sequence to enu sequence with coordinate averaging
    for i in 0..pos.time.len() as usize {
        let mut delta = [0.0; 3];
        let mut blh = [0.0; 3];
        for j in 0..3 {
            delta[j] = pos.coor[j][i] - pos.ave_pos[j];
        }
        ecef2pos(&[pos.coor[0][i], pos.coor[1][i], pos.coor[2][i]], &mut blh);
        let enu = ecef2enu(&blh, &delta);
        for j in 0..3 {
            pos.coor[j][i] = enu[j];
        }
    }

    let mut mean = [0.0_f64; 3];
    for j in 0..3 {
        mean[j] = pos.coor[j].iter().sum::<f64>() / pos.coor[j].len() as f64;
    }

    let mut var = [0.0_f64; 3];
    for j in 0..3 {
        var[j] = pos.coor[j]
            .iter()
            .map(|&x| (x - mean[j]).powi(2))
            .sum::<f64>()
            / pos.coor[j].len() as f64;
    }

    for j in 0..3 {
        enu[j] = var[j].sqrt();
    }
}

// Calculate the noise of pseudorange or carrier
fn cal_noise(
    noise: &mut Noise,
    data: [f64; NFREOBS],
    time_change: bool,
    slip: [u8; NFREQ + NEXOBS],
) {
    for j in 0..NFREOBS {
        // Data discontinuity, initialization of stored predecessor variables
        if data[j] == 0.0 || time_change || slip[j] !=0 {
            noise.single_diff[j] = 0.0;
            noise.double_diff[j] = 0.0;
            noise.form[j] = 0.0;
            noise.form_single_diff[j] = 0.0;
            noise.form_double_diff[j] = 0.0;
            continue;
        }
        // singleton processing
        if noise.form[j] != 0.0 {
            noise.form_single_diff[j] = noise.single_diff[j];
            noise.single_diff[j] = data[j] - noise.form[j];
        }
        noise.form[j] = data[j];
        // double differential processing
        if noise.form_single_diff[j] != 0.0 {
            noise.form_double_diff[j] = noise.double_diff[j];
            noise.double_diff[j] = noise.single_diff[j] - noise.form_single_diff[j];
        }
        // cope with the three evils
        if noise.form_double_diff[j] != 0.0 {
            noise.triple_diff[j].push(noise.double_diff[j] - noise.form_double_diff[j]);
            // Triple Difference Constant Value Judgment
            if noise.triple_diff[j].last().unwrap_or(&0.0).abs() > 100.0 {
                noise.single_diff[j] = 0.0;
                noise.double_diff[j] = 0.0;
                noise.form_single_diff[j] = 0.0;
                noise.form_double_diff[j] = 0.0;
                noise.triple_diff[j].pop();
            }
        }
    }
}

// get information per epoch
pub fn getepdata(
    sysinfo: &mut Vec<SysInfo>,
    satinfo: &mut Vec<SatInfo>,
    siginfo: &mut Vec<Vec<SigInfo>>,
    rtk: &Rtk,
    obs: &Vec<Obs>,
    ti: i32,
    min_ele: &mut f32,
    n: usize,
) {
    let mut satnum = [0.0; NUMSYS];
    /* Store ephtime */
    /* Store azel dop resp snr xCosv xPhsv */
    for i in 0..n {
        let satidx = obs[i].sat - 1;
        let syslog = satsys(satidx + 1).ilog2() as usize;
        let mut pnum = 0;
        let mut lnum = 0;

        if satinfo[satidx].verify {
            /* Store MinEle */
            if rtk.ssat[satidx].azel[1] > 0.0
                && (*min_ele == 0.0 || rtk.ssat[satidx].azel[1] * 180.0 / PI < *min_ele as f64)
            {
                *min_ele = (rtk.ssat[satidx].azel[1] * 180.0 / PI) as f32;
            }

            /* Store azel for each epoch by satellite */
            if !rtk.ssat[satidx].azel[1].is_nan()
            {
                satinfo[satidx].azel[0].push((rtk.ssat[satidx].azel[0] * 180.0 / PI) as f32);
                satinfo[satidx].azel[1].push((rtk.ssat[satidx].azel[1] * 180.0 / PI) as f32);
                satinfo[satidx].aetime.push(obs[i].time);
                /* Store L1 snr for each satellite epoch */
                satinfo[satidx].l1snr.push(obs[i].snr[0] as f32 / 1000.0);
            }

            /* Cal Pseudo-range diff */
            let time_change = ti > 0 && obs[i].time.time - satinfo[satidx].pretime.time > ti as i64;
            let slip=satinfo[satidx].slip.clone();
            cal_noise(&mut satinfo[satidx].pnoise, obs[i].p, time_change.clone(),slip);
            /* Cal Carrier diff */
            let mut l = [0.0; NFREOBS];
            for j in 0..NFREOBS {
                if obs[i].l[j] == 0.0 { continue; }
                l[j] = obs[i].l[j] * CLIGHT / sat2freq(obs[i].sat, obs[i].code[j], None);
            }
            cal_noise(&mut satinfo[satidx].lnoise, l, time_change,slip);

            /* Store resp averages by satellite, 0:resp 1:respnum */
            satinfo[satidx].resp[0] = (satinfo[satidx].resp[0] * satinfo[satidx].resp[1]
                + rtk.ssat[satidx].resp[0] as f32)
                / (satinfo[satidx].resp[1] + 1.0);
            satinfo[satidx].resp[1] += 1.0;

            satinfo[satidx].eptime.push(obs[i].time);

            for j in 0..(NFREQ + NEXOBS) {
                if obs[i].p[j] != 0.0 {
                    pnum += 1;
                }
                if obs[i].l[j] != 0.0 {
                    lnum += 1;
                }
                if obs[i].snr[j] != 0 {
                    /* Store mean snr of each satellite signal */
                    siginfo[syslog][j].snr = (siginfo[syslog][j].snr
                        * siginfo[syslog][j].snrnum as f32
                        + obs[i].snr[j] as f32 / 1000.0)
                        / (siginfo[syslog][j].snrnum as f32 + 1.0);
                    siginfo[syslog][j].snrnum += 1;
                }
            }

            /* Store xCosv、xPhsv */
            if pnum == 1 {
                sysinfo[syslog].x_co_sv += 1;
            }
            if lnum == 1 {
                sysinfo[syslog].x_ph_sv += 1;
            }

            satnum[syslog] += 1.0;
        }
    }

    for i in 0..NUMSYS {
        if satnum[i] == 0.0 { continue; }
        sysinfo[i].epnum += 1;
        sysinfo[i].meansat = (sysinfo[i].meansat * (sysinfo[i].epnum - 1) as f32 + satnum[i]) / sysinfo[i].epnum as f32;
    }
}

// Count Exp/Hav data for follow calculation
pub fn count_exp_hav(
    sysinfo: &mut Vec<SysInfo>,
    satinfo: &mut Vec<SatInfo>,
    siginfo: &mut Vec<Vec<SigInfo>>,
    rtk: &Rtk,
    obs: &Vec<Obs>,
    nav: &mut Nav,
    satlist: &Vec<usize>,
    satcount: usize,
    n: usize,
    el: f64,
) {
    for i in 0..satcount {
        let sys = satsys(satlist[i]);
        let syslog = sys.ilog2() as usize;
        let satidx = satlist[i] - 1;

        let mut found_obs = false;

        for j in 0..n {
            if satlist[i] == obs[j].sat {
                /* obs exist in the current epoch of the satellite */
                satinfo[satidx].exp_obs += 1;
                satinfo[satidx].hav_obs += 1;
                /* Statistical pseudo-distance, carrier and snr by satellite and frequency */
                sigcount(siginfo, &obs[j], &sysinfo[syslog].code, false);

                if rtk.ssat[satidx].azel[1] >= el * PI / 180.0 {
                    satinfo[satidx].exp_obs_el += 1;
                    satinfo[satidx].hav_obs_el += 1;
                    sigcount(siginfo, &obs[j], &sysinfo[syslog].code, true);
                }
                found_obs = true;
                break;
            }
        }

        if !found_obs {
            let mut dt = 0.0;
            if !ephclk(obs[0].time, obs[0].time, satlist[i], nav, &mut dt) {
                continue;
            }

            let corrtime = timeadd(obs[0].time, -dt);
            let mut rs = [0.0; 6];
            let mut dts = [0.0; 2];
            let mut var = 0.0;
            let mut svh = 0;
            let mut pos = [0.0; 3];
            let mut e = [0.0; 3];
            let mut azel = [0.0; 2];

            if !satpos(
                corrtime,
                obs[0].time,
                satlist[i],
                nav,
                &mut rs,
                &mut dts,
                &mut var,
                &mut svh,
            ) {
                continue;
            }

            ecef2pos(&rtk.sol.rr, &mut pos);
            geodist(&Vector3::new(rs[0], rs[1], rs[2]), &rtk.sol.rr, &mut e);
            satazel(&pos, &e, Some(&mut azel));

            if azel[1] > 0.0 {
                /* sat nav exist but no sat obs */
                satinfo[satidx].exp_obs += 1;
                /* corresponding sys sigobs judgment */
                for k in 0..(NFREQ + NEXOBS) {
                    if (sys == SYS_GPS && code2idx(SYS_GPS, sysinfo[syslog].code[k]) == 2 && !gps_l5_check(satlist[i])) /* This GPS sat does not have L5 freq.b */
                        || ((satlist[i] <= 121 && satlist[i] > 105) &&
                        bds_sig_check(sysinfo[syslog].code[k] as usize) == 3)  /* This BDS sat is not a BDS3 sat */
                        || ((satlist[i] > 121) && bds_sig_check(sysinfo[syslog].code[k] as usize) == 2)
                    /* This BDS sat is not a BDS2 sat */
                    {
                        continue;
                    }
                    if sysinfo[syslog].code[k] != 0 {
                        siginfo[syslog][k].exp_obs += 1;
                    }
                }
            }
            if azel[1] >= el * PI / 180.0 {
                satinfo[satidx].exp_obs_el += 1;
                for k in 0..(NFREQ + NEXOBS) {
                    if (sys == SYS_GPS
                        && code2idx(SYS_GPS, sysinfo[syslog].code[k]) == 2
                        && !gps_l5_check(satlist[i]))
                        || ((satlist[i] <= 121 && satlist[i] > 105)
                        && bds_sig_check(sysinfo[syslog].code[k] as usize) == 3)
                        || ((satlist[i] > 121)
                        && bds_sig_check(sysinfo[syslog].code[k] as usize) == 2)
                    {
                        continue;
                    }
                    if sysinfo[syslog].code[k] != 0 {
                        siginfo[syslog][k].exp_obs_el += 1;
                    }
                }
            }
        }
    }
}

// calculate some statistics
pub fn cal_symbol(
    totinfo: &mut TotInfo,
    sysinfo: &mut Vec<SysInfo>,
    satinfo: &mut Vec<SatInfo>,
    siginfo: &mut Vec<Vec<SigInfo>>,
    navsys: usize,
) {
    let mut pnoise_list: [[Vec<f64>; NFREOBS]; NUMSYS] = array_init(|_| array_init(|_| Vec::new()));
    let mut lnoise_list= pnoise_list.clone();
    for i in 0..MAXSAT {
        let sys = satsys(i + 1);
        let syslog = (sys as f64).log2() as usize;

        // Count obs/slip
        if satinfo[i].slipnum > 0 || satinfo[i].obsnum > 0 {
            satinfo[i].sat_rt = if satinfo[i].exp_obs > 0 {
                satinfo[i].hav_obs as f32 / satinfo[i].exp_obs as f32 * 100.0
            } else {
                0.0
            };
            satinfo[i].sat_el_rt = if satinfo[i].exp_obs_el > 0 {
                satinfo[i].hav_obs_el as f32 / satinfo[i].exp_obs_el as f32 * 100.0
            } else {
                0.0
            };
            sysinfo[syslog].satnum += 1;
            sysinfo[syslog].slipnum += satinfo[i].slipnum;
            totinfo.totobs += satinfo[i].obsnum;
            totinfo.totslip += satinfo[i].slipnum;

            satinfo[i].oslps = if satinfo[i].slipnum > 0 {
                satinfo[i].obsnum / satinfo[i].slipnum
            } else {
                satinfo[i].obsnum
            };
        }

        // Count effective sat
        if satinfo[i].exp_obs > 0 {
            totinfo.exp.totsatobs += satinfo[i].exp_obs;
            totinfo.hav.totsatobs += satinfo[i].hav_obs;
        }
        if satinfo[i].exp_obs_el > 0 {
            totinfo.exp.totsatobs_el += satinfo[i].exp_obs_el;
            totinfo.hav.totsatobs_el += satinfo[i].hav_obs_el;
        }

        // Cal the signal pseudo-range noise

        if satinfo[i].pnoise.triple_diff.iter().any(|vec| !vec.is_empty()) {
            for j in 0..NFREOBS {
                let sum = satinfo[i].pnoise.triple_diff[j].iter().map(|&x| x * x).sum::<f64>();
                if sum == 0.0 { continue; }
                let prange_noise = (sum / satinfo[i].pnoise.triple_diff[j].len() as f64 / 8.0_f64).sqrt();
                pnoise_list[syslog][j].push(prange_noise);
            }
        }
        // Cal the signal carrier noise
        if satinfo[i].lnoise.triple_diff.iter().any(|vec| !vec.is_empty()) {
            for j in 0..NFREOBS {
                let sum = satinfo[i].lnoise.triple_diff[j].iter().map(|&x| x * x).sum::<f64>();
                if sum == 0.0 { continue; }
                // for k in &satinfo[i].lnoise.triple_diff[j]{
                //     println!("{:.2} ",k)
                // }
                // println!("sum:{:.2}",sum);
                let carrier_noise = (sum / satinfo[i].lnoise.triple_diff[j].len() as f64 / 8.0_f64).sqrt();
                lnoise_list[syslog][j].push(carrier_noise);
            }
        }
    }

    // Count sys obs-Exp/Hav
    for i in 0..6 {
        for j in 0..(NFREQ + NEXOBS) {
            if siginfo[i][j].exp_obs > 0 {
                sysinfo[i].exp_obs += siginfo[i][j].exp_obs;
                sysinfo[i].hav_obs += siginfo[i][j].hav_obs;
            }
        }

        // Cal sys slip
        if sysinfo[i].exp_obs > 0 {
            sysinfo[i].oslps = if sysinfo[i].slipnum > 0 {
                sysinfo[i].hav_obs / sysinfo[i].slipnum
            } else {
                sysinfo[i].hav_obs
            };
        }
    }

    // To be total obs data
    for i in 0..6 {
        let sys = 1 << i;
        for j in 0..(NFREQ + NEXOBS) {
            if (sys & navsys) != 0 && sysinfo[i].code[j] != 0 {
                totinfo.hav.totobs += siginfo[i][j].hav_obs;
                totinfo.hav.totobs_el += siginfo[i][j].hav_obs_el;
                totinfo.exp.totobs += siginfo[i][j].exp_obs;
                totinfo.exp.totobs_el += siginfo[i][j].exp_obs_el;
            }
        }
    }

    // Get statistic needed for the noise box plot
    for i in 0..NUMSYS{
        for j in 0..NFREOBS{
            siginfo[i][j].pnoise_box=cal_box(&mut pnoise_list[i][j]).unwrap_or([0.0;5]);
            siginfo[i][j].lnoise_box=cal_box(&mut lnoise_list[i][j]).unwrap_or([0.0;5]);
            siginfo[i][j].pnoise=Some(pnoise_list[i][j].iter().sum::<f64>()/pnoise_list[i][j].len() as f64).unwrap_or(0.0);
            siginfo[i][j].lnoise=Some(lnoise_list[i][j].iter().sum::<f64>()/lnoise_list[i][j].len() as f64).unwrap_or(0.0);
        }
    }
}

// count the different levels of mp
fn sum_mp(sysinfo: &mut Vec<SysInfo>, satinfo: &mut Vec<SatInfo>) {
    let mut sys_count = [0; 6];

    for i in 0..MAXSAT {
        for j in 0..(NFREQ + NEXOBS) {
            if satinfo[i].sat_mp[j] != 0.0 {
                satinfo[i].freq[j] = codefreq(satinfo[i].code[j]);

                let mut k = j as isize - 1;
                while k >= 0 {
                    if satinfo[i].freq[j] == satinfo[i].freq[k as usize] {
                        satinfo[i].freq[j] = 0; // Avoid output of different signal mp at the same frequency
                        break;
                    }
                    k -= 1;
                }

                let sys = satsys(i + 1);
                let syslog = (sys as f64).log2() as usize;
                sysinfo[syslog].totmp[j] = (sysinfo[syslog].totmp[j] * sys_count[syslog] as f32
                    + satinfo[i].sat_mp[j] * satinfo[i].obsnum as f32)
                    / (sys_count[syslog] as f32 + satinfo[i].obsnum as f32);
                sys_count[syslog] += satinfo[i].obsnum;
            }
        }
    }
}

// ----------------------------
// Other functions
// ----------------------------
pub fn namebysyslog(sys: usize) -> &'static str {
    match sys {
        0 => "GPS",
        1 => "GLO",
        2 => "GAL",
        3 => "QZS",
        4 => "BDS",
        _ => ""
    }
}

fn get_dist_k(que: &VecDeque<f32>, idx: usize) -> f32 {
    // 动态设定k值
    let k = que.len().sqrt();
    // 取出待估点
    let mut new_que = que.clone();
    let point = new_que.remove(idx).unwrap();
    // 计算可达距离
    let mut distances: Vec<f32> = new_que.iter().map(|&x| (x - point).abs()).collect();
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let neighbors: Vec<f32> = distances.iter().take(k).cloned().collect();
    *neighbors.last().unwrap()
}

fn get_index_k(que: &VecDeque<f32>, idx: usize) -> Vec<usize> {
    // 动态设定k值
    let k = que.len().sqrt();
    // 取出待估点
    let mut new_que = que.clone();
    let point = new_que.remove(idx).unwrap();
    // 提取k邻域内的值下标
    let mut distances: Vec<(usize, f32)> = new_que.iter().enumerate()
        .map(|(i, &x)| (i, (x - point).abs()))
        .collect();
    distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    distances.iter().take(k).map(|&(index, _)| index).collect()
}
#[allow(dead_code)]
fn lof(que: &VecDeque<f32>) -> bool {
    // 动态设定k值
    let k = que.len().sqrt();
    // 提取k邻域内的下标以获得值
    let k_index = get_index_k(&que, que.len() - 1);
    // 统计新入mp值k邻域内的局部可达密度
    let mut sum_reach_dist = 0.0;
    for i in 0..k {
        let k_dist = get_dist_k(&que, k_index[i]);
        sum_reach_dist += k_dist.max((que[k_index[i]] - que[que.len() - 1]).abs());
    }
    sum_reach_dist += 1E-10;
    let lrd_p = k as f32 / sum_reach_dist;

    //统计mp值k邻域内的所有点的平均可达密度
    let mut sum_lrd = 0.0;
    for i in 0..k {
        let o_k_index = get_index_k(&que, k_index[i]);
        sum_reach_dist = 0.0;
        for j in 0..k {
            let k_dist = get_dist_k(&que, o_k_index[j]);
            sum_reach_dist += k_dist.max((que[o_k_index[j]] - que[k_index[i]]).abs());
        }
        sum_reach_dist += 1E-10;
        sum_lrd += k as f32 / sum_reach_dist;
    }

    let lof = sum_lrd / k as f32 / lrd_p;

    println!("{}", lof);

    lof > 2.0
}

// update nav when data from real-time stream
fn update_rtnav(navs: &mut Nav, ref_time: GTime) {
    let eph_to_remove: Vec<usize> = navs
        .eph
        .iter()
        .enumerate()
        .filter(|(_, eph)| timediff(eph.toe, ref_time) > 4.0 * 3600.0)
        .map(|(index, _)| index)
        .collect();

    for index in eph_to_remove.iter().rev() {
        navs.eph.remove(*index);
        navs.n -= 1;
    }

    let geph_to_remove: Vec<usize> = navs
        .geph
        .iter()
        .enumerate()
        .filter(|(_, geph)| timediff(geph.toe, ref_time) > 4.0 * 3600.0)
        .map(|(index, _)| index)
        .collect();

    for index in geph_to_remove.iter().rev() {
        navs.geph.remove(*index);
        navs.ng -= 1;
    }
}


// ----------------------------
// output function
// ----------------------------
// output total information
fn out_totsum<W: Write>(totinfo: &TotInfo, nav_flag: bool, ti: i32, el: f64, writer: &mut W) -> Result<(), Error> {
    let firstepoch = time2str(totinfo.firstepoch);
    let lastepoch = time2str(totinfo.lastepoch);
    writeln!(
        writer,
        "{:<7}{:>21}{:>21}{:>7}{:>8}{:>8}{:>11}{:>9}{:>10}{:<2}{:>8}{:>9}{:>16}",
        "#TOTSUM", "First_Epoch________", "Last_Epoch_________", "Hours", "Sample", "MinEle",
        "Epoch_%Rt", "Obs_%Rt", "Obs_%Rt>", el, "Slip", "o/slps", "SPP_STD"
    )?;
    if nav_flag {
        writeln!(
            writer,
            "=TOTSUM{:>21}{:>21}{:7.2}{:8}{:8.2}{:11.2}{:9.2}{:12.2}{:8.0}{:9}  E:{:<4.2}  N:{:<4.2}  U:{:<4.2}",
            firstepoch, lastepoch,
            (totinfo.lastepoch.time - totinfo.firstepoch.time) as f32 / 3600.0, ti, totinfo.minele,
            totinfo.ep_rt, totinfo.obs_rt, totinfo.obs_el_rt, totinfo.totslip, totinfo.oslps,
            totinfo.enu_std[0], totinfo.enu_std[1], totinfo.enu_std[2]
        )?;
    } else {
        writeln!(
            writer,
            "=TOTSUM{:>21}{:>21}{:>7.2}{:>8}{:>8}{:11.2}{:9.2}{:>12}{:8}{:9}{:>16}",
            firstepoch, lastepoch,
            (totinfo.lastepoch.time - totinfo.firstepoch.time) as f32 / 3600.0, ti, "-",
            totinfo.ep_rt, totinfo.obs_rt, "-", totinfo.totslip, totinfo.oslps, "-"
        )?;
    }
    Ok(())
}

// output system information
pub fn out_syssum<W: Write>(sysinfo: &Vec<SysInfo>, writer: &mut W) -> Result<(), Error> {
    let mut codecopy = [[0; NFREQ + NEXOBS]; 6];

    writeln!(writer, "\n{:<7}{:>8}{:>8}{:>9}{:>9}{:>5}{:>7}{:>7}{:>6}{:>8}{:>22}",
             "#SYSSUM", "ExpObs", "HavObs", "Obs_%Rt", "Meansat", "Sig", "xCoSv", "xPhSv", "Slip", "o/slps", "MPx")?;

    for i in 0..6 {
        for j in 0..(NFREQ + NEXOBS) {
            codecopy[i][j] = sysinfo[i].code[j];
        }
    }

    for i in 0..6 {
        let mut mpx_flag = 0;
        let sysname = namebysyslog(i);

        if sysinfo[i].exp_obs > 0 {
            write!(writer, "={}SUM{:8.0}{:8.0}{:9.2}{:9.2}{:5}{:7}{:7}{:6}{:8.0}    ",
                   sysname, sysinfo[i].exp_obs, sysinfo[i].hav_obs,
                   sysinfo[i].hav_obs as f32 / sysinfo[i].exp_obs as f32 * 100.0,
                   sysinfo[i].meansat, sysinfo[i].codenum, sysinfo[i].x_co_sv,
                   sysinfo[i].x_ph_sv, sysinfo[i].slipnum,
                   sysinfo[i].oslps)?;
            let mut j = 0;
            while j < NFREOBS {
                let mut loop_flag = true;
                if codecopy[i][j] != 0 {
                    let mut k = j + 1;
                    let current_freq = codefreq(codecopy[i][j]);
                    while k < (NFREQ + NEXOBS) {
                        if codefreq(codecopy[i][k]) > 0 &&
                            current_freq > codefreq(codecopy[i][k]) {
                            break;
                        }
                        k += 1;
                    }
                    if k == (NFREQ + NEXOBS) {
                        if current_freq > mpx_flag {
                            if sysinfo[i].totmp[j] < 0.001 {
                                write!(writer, "mp{}:{:<8}", current_freq, "-")?;
                            } else {
                                write!(writer, "mp{}:{:<8.3}", current_freq, sysinfo[i].totmp[j])?;
                            }
                            mpx_flag = current_freq;
                        }
                    } else {
                        if codefreq(codecopy[i][k]) > mpx_flag {
                            if i != 5 || code2obs(codecopy[i][k]).chars().next().unwrap() != '1' {
                                if sysinfo[i].totmp[k] < 0.001 {
                                    write!(writer, "mp{}:{:<8}", codefreq(codecopy[i][k]), "-")?;
                                } else {
                                    write!(writer, "mp{}:{:<8.3}", codefreq(codecopy[i][k]), sysinfo[i].totmp[k])?;
                                }
                                mpx_flag = codefreq(codecopy[i][k]);
                            }
                        }
                        codecopy[i][k] = 0;
                        loop_flag = false;
                    }
                }
                if loop_flag {
                    j += 1;
                }
            }
            writeln!(writer, "")?;
        }
    }
    Ok(())
}

// output satellite information
fn out_sat<W: Write>(satinfo: &mut Vec<SatInfo>, el: f64, writer: &mut W) -> Result<(), Error> {
    write!(
        writer,
        "\n{:<7}{:>8}{:>9}{:>6}{:>8}{:>9}{:>10}{:<2}{:>17}",
        "#SatPRN", "Obsnum", "Slipnum", "nGap", "o/slps", "Sat_%Rt", "Sat_%Rt>", el, "MPx"
    )?;

    for i in 0..MAXSAT {
        if satinfo[i].obsnum != 0 || satinfo[i].slipnum != 0 {
            let id = satno2id(i + 1);

            write!(
                writer,
                "\n{:7}{:8}{:9}{:6}{:8}",
                id,
                satinfo[i].obsnum,
                satinfo[i].slipnum,
                satinfo[i].ngap,
                satinfo[i].oslps
            )?;

            if satinfo[i].sat_rt.is_nan() {
                write!(writer, "{:>9}", "-")?;
            } else {
                write!(writer, "{:9.2}", satinfo[i].sat_rt)?;
            }

            if satinfo[i].sat_el_rt.is_nan() {
                write!(writer, "{:>12}", "-")?;
            } else {
                write!(writer, "{:12.2}", satinfo[i].sat_el_rt)?;
            }

            write!(writer, "  ")?;

            let mut j: i32 = 0;
            while j < NFREOBS as i32 {
                if satinfo[i].freq[j as usize] != 0 {
                    let mut k = j as usize + 1;
                    while k < (NFREQ + NEXOBS) {
                        if satinfo[i].freq[j as usize] > satinfo[i].freq[k] && satinfo[i].freq[k] > 0 {
                            break;
                        }
                        k += 1;
                    }
                    if k == (NFREQ + NEXOBS) {
                        write!(writer, "mp{}:{:<6.3}", satinfo[i].freq[j as usize], satinfo[i].sat_mp[j as usize])?;
                    } else {
                        write!(writer, "mp{}:{:<6.3}", satinfo[i].freq[k], satinfo[i].sat_mp[k])?;
                        satinfo[i].freq[k] = 0;
                        j -= 1;
                    }
                }
                j += 1;
            }
        }
    }
    Ok(())
}

// output signal information
fn out_signal<W: Write>(
    sysinfo: &Vec<SysInfo>,
    siginfo: &Vec<Vec<SigInfo>>,
    navsys: usize,
    el: f64,
    writer: &mut W,
) -> Result<(), Error> {
    write!(
        writer,
        "\n\n{:<4}{:>10}{:>5}{:>8}{:>8}{:>8}{:>6}{:<2}{:>6}{:<2}{:>6}{:<2}{:>5}{:>5}{:>24}{:>39}",
        "#SYS", "Freq.Band", "CODE", "ExpObs", "HavObs", "%Ratio",
        "Exp>", el, "Hav>", el, "%Rt>", el, "SNR", "MP", "PNoise", "LNoise"
    )?;
    for i in 0..6 {
        let sysname = namebysyslog(i);
        let sys = pow(2, i);
        for j in 0..(NFREQ + NEXOBS) {
            if (sys & navsys) != 0 && sysinfo[i].code[j] != 0 {
                write!(
                    writer,
                    "\n={}{:>10}{:>5}{:8}{:8}{:8.2}",
                    sysname,
                    getfreqband(sys, sysinfo[i].code[j]),
                    OBS_CODES[sysinfo[i].code[j] as usize],
                    siginfo[i][j].exp_obs,
                    siginfo[i][j].hav_obs,
                    siginfo[i][j].hav_obs as f32 / siginfo[i][j].exp_obs as f32 * 100.0
                )?;
                if siginfo[i][j].exp_obs_el != 0 {
                    write!(
                        writer,
                        "{:8}{:8}{:8.2}",
                        siginfo[i][j].exp_obs_el,
                        siginfo[i][j].hav_obs_el,
                        siginfo[i][j].hav_obs_el as f32 / siginfo[i][j].exp_obs_el as f32 * 100.0
                    )?;
                } else {
                    write!(writer, "{:>8}{:>8}{:>8}", "-", "-", "-")?;
                }
                write!(writer, "{:5.1}{:6.3}",siginfo[i][j].snr, sysinfo[i].totmp[j])?;
                write!(
                    writer," [{:<5.3}{:6.3}{:6.3}{:6.3}{:6.3}{:6.3}] [{:<5.3}{:6.3}{:6.3}{:6.3}{:6.3}{:6.3}]",
                    siginfo[i][j].pnoise_box[0],siginfo[i][j].pnoise_box[1],
                    siginfo[i][j].pnoise_box[2],siginfo[i][j].pnoise_box[3],
                    siginfo[i][j].pnoise_box[4],siginfo[i][j].pnoise,
                    siginfo[i][j].lnoise_box[0],siginfo[i][j].lnoise_box[1],
                    siginfo[i][j].lnoise_box[2],siginfo[i][j].lnoise_box[3],
                    siginfo[i][j].lnoise_box[4],siginfo[i][j].lnoise
                )?;
            }
        }
    }
    Ok(())
}

// output information of azimuth/elevation angle
fn out_skymap<W: Write>(
    satinfo: &Vec<SatInfo>,
    firstepoch: GTime,
    ti: i32,
    writer: &mut W,
) -> Result<(), Error> {
    writeln!(writer, "\n\n{}{:>5}{:>10}{:>9}{:>11}{:>5}",
             "#Sky", "prn", "epochnum", "azimuth", "elevation", "snr")?;

    for i in 0..MAXSAT {
        let aetime = satinfo[i].aetime.clone();
        let azel = satinfo[i].azel.clone();
        let aenum = aetime.len();

        // if i != 3 {continue}
        // let id = satno2id(i + 1);
        // for j in 0..aenum {
        //     let epochnum = (aetime[j].time - firstepoch.time) as i32 / ti;
        //     if !azel[1][j].is_nan() && azel[0][j] > 0.0 && azel[1][j] > 0.0 {
        //         { write!(writer, "   ")?; }
        //         writeln!(writer, "{:>6}{:>10}{:>9.2}{:>11.2}{:>5}",
        //                  id, epochnum, azel[0][j], azel[1][j], satinfo[i].l1snr[j])?;
        //     }
        // }
        // break;


        if aenum == 0 {
            continue;
        }
        let id = satno2id(i + 1);
        let mut formout = 0;
        let mut arcout = true;

        if aenum <= 80 {
            for j in 0..aenum {
                let epochnum = (aetime[j].time - firstepoch.time) as i32 / ti;
                if formout == 0 { formout = epochnum }
                if !azel[1][j].is_nan() && azel[0][j] > 0.0 && azel[1][j] > 0.0 {
                    if epochnum - formout > 10 * 60 / ti { arcout = true }  /* epoch break off more than 10min */
                    if arcout {
                        write!(writer, "arc")?;  /* output arc to indicate one arc begin */
                        arcout = false;
                    } else { write!(writer, "   ")?; }
                    writeln!(writer, "{:>6}{:>10}{:>9.2}{:>11.2}{:>7.3}",
                             id, epochnum, azel[0][j], azel[1][j], satinfo[i].l1snr[j])?;
                    formout = epochnum;
                }
            }
        } else {
            // Stores the starting point of each segment
            let mut stage = vec![0];
            // An interruption of more than 10 minutes is considered a new segment
            for j in 0..aenum - 1 {
                if timediff(aetime[j + 1], aetime[j]) >= 10.0 * 60.0 {
                    stage.push(j + 1);
                }
            }

            for n in 0..stage.len() {
                // get end of each segment
                let end = if n == stage.len() - 1 { aenum } else { stage[n + 1] };
                // get pumping step of each segment
                let step = if end - stage[n] <= 80 { 1.0 } else { (end - stage[n]) as f32 / 80.0 };
                // check if this segment all 0
                let mut zero_flag = true;
                for j in 0..min(80, end - stage[n]) {
                    let k = stage[n] + (j as f32 * step) as usize;
                    if azel[0][k] > 0.0 && azel[1][k] > 0.0 {
                        zero_flag = false;
                        break;
                    }
                }
                if zero_flag { continue; }
                arcout = true;
                // output
                for j in 0..min(80, end - stage[n]) {
                    let k = stage[n] + (j as f32 * step) as usize;
                    let epochnum = (aetime[k].time - firstepoch.time) as i32 / ti;
                    if !azel[1][k].is_nan() && azel[0][k] > 0.0 && azel[1][k] > 0.0 {
                        if arcout {
                            write!(writer, "arc")?;  /* output arc to indicate one arc begin */
                            arcout = false;
                        } else { write!(writer, "   ")?; }
                        writeln!(writer, "{:>6}{:>10}{:>9.2}{:>11.2}{:>7.3}",
                                 id, epochnum, azel[0][k], azel[1][k], satinfo[i].l1snr[k])?;
                    }
                    // output the last one
                    if j == 79 && k < end - 1 && !arcout {
                        let epochnum = (aetime[end - 1].time - firstepoch.time) as i32 / ti;
                        if !azel[1][end - 1].is_nan() && azel[0][end - 1] > 0.0 && azel[1][end - 1] > 0.0 {
                            writeln!(writer, "{:>9}{:>10}{:>9.2}{:>11.2}{:>7.3}",
                                     id, epochnum, azel[0][end - 1], azel[1][end - 1], satinfo[i].l1snr[end - 1])?;
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

// output information of position
fn out_pos<W: Write>(
    pos: &Pos,
    firstepoch: GTime,
    ti: i32,
    writer: &mut W,
) -> Result<(), Error> {
    writeln!(
        writer, "\n{:<6}{}{:<15.4}{}{:<15.4}{}{:<15.4}", "#POS",
        "x:", pos.ave_pos[0], "y:", pos.ave_pos[1], "z:", pos.ave_pos[2]
    )?;
    writeln!(
        writer, "{:<8}{:>7}{:>7}{:>7}{:>7}{:>8}{:>8}{:>8}",
        "epochnum", "GDOP", "PDOP", "HDOP", "VDOP",
        "SPP_E", "SPP_N", "SPP_U"
    )?;
    let epnum = pos.time.len();
    /* output directly without pumping as epnum<400 */
    if epnum <= 800 {
        for i in 0..epnum {
            let epochnum = (pos.time[i].time - firstepoch.time) as i32 / ti;
            if pos.dop[1][i] != 0.0 {
                writeln!(
                    writer, "{:<8}{:>7.2}{:>7.2}{:>7.2}{:>7.2}{:>8.3}{:>8.3}{:>8.3}",
                    epochnum, pos.dop[0][i], pos.dop[1][i], pos.dop[2][i], pos.dop[3][i],
                    pos.coor[0][i], pos.coor[1][i], pos.coor[2][i]
                )?;
            }
        }
    } else {
        /* Thinning to 400 */
        let step = epnum as f32 / 800.0;
        for i in 0..800 {
            let k = (i as f32 * step) as usize;
            let epochnum = (pos.time[k].time - firstepoch.time) as i32 / ti;
            if pos.dop[1][k] != 0.0 {
                writeln!(
                    writer, "{:<8}{:>7.2}{:>7.2}{:>7.2}{:>7.2}{:>8.3}{:>8.3}{:>8.3}",
                    epochnum, pos.dop[0][k], pos.dop[1][k], pos.dop[2][k], pos.dop[3][k],
                    pos.coor[0][k], pos.coor[1][k], pos.coor[2][k]
                )?;
            }
        }
    }
    Ok(())
}

// output information of time series
fn out_ts<W: Write>(
    satinfo: &Vec<SatInfo>,
    firstepoch: GTime,
    ti: i32,
    writer: &mut W,
) -> Result<(), Error> {
    writeln!(writer, "\n{}{:>5}{:>10}{:>9}{:>10}",
             "#list", "prn", "epochnum", "obsflag", "slipflag")?;
    for i in 0..MAXSAT {
        let epnum = satinfo[i].eptime.len();
        if epnum > 0 {
            let id = satno2id(i + 1);
            for j in 0..epnum {
                if satinfo[i].obs_flag[j] != 0 || satinfo[i].slip_flag[j] != 0 {
                    if j == 0 { write!(writer, "arc")?; } else { write!(writer, "   ")?; }
                    writeln!(writer, "{:>7}{:>10}{:>9}{:>10}", id,
                             (satinfo[i].eptime[j].time - firstepoch.time) as i32 / ti,
                             satinfo[i].obs_flag[j], satinfo[i].slip_flag[j]
                    )?;
                }
            }
        }
    }
    Ok(())
}

// output total information[real-time]
fn out_totsum_rt<W: Write>(totinfo: &TotInfo, nav_flag: bool, ti: i32, el: f64, pos: &[f64], writer: &mut W) -> Result<(), Error> {
    let firstepoch = time2str(totinfo.firstepoch);
    let lastepoch = time2str(totinfo.lastepoch);
    writeln!(
        writer,
        "{:<7}{:>21}{:>21}{:>7}{:>8}{:>8}{:>11}{:>9}{:>10}{:<2}{:>8}{:>9}{:>25}",
        "#TOTSUM", "First_Epoch________", "Last_Epoch_________", "Hours", "Sample", "MinEle",
        "Epoch_%Rt", "Obs_%Rt", "Obs_%Rt>", el, "Slip", "o/slps", "SPP"
    )?;
    if nav_flag {
        writeln!(
            writer,
            "=TOTSUM{:>21}{:>21}{:7.2}{:8}{:8.2}{:11.2}{:9.2}{:12.2}{:8.0}{:9}  X:{:<.2}  Y:{:<.2}  Z:{:<.2}",
            firstepoch, lastepoch,
            (totinfo.lastepoch.time - totinfo.firstepoch.time) as f32 / 3600.0, ti, totinfo.minele,
            totinfo.ep_rt, totinfo.obs_rt, totinfo.obs_el_rt, totinfo.totslip, totinfo.oslps,
            pos[0], pos[1], pos[2]
        )?;
    } else {
        writeln!(
            writer,
            "=TOTSUM{:>21}{:>21}{:>7.2}{:>8}{:>8}{:11.2}{:9.2}{:>12}{:8}{:9}{:>25}",
            firstepoch, lastepoch,
            (totinfo.lastepoch.time - totinfo.firstepoch.time) as f32 / 3600.0, ti, "-",
            totinfo.ep_rt, totinfo.obs_rt, "-", totinfo.totslip, totinfo.oslps, "-"
        )?;
    }
    Ok(())
}

// output satellite information[real-time]
fn out_sat_rt<W: Write>(qc: &QC, writer: &mut W) -> Result<(), Error> {
    writeln!(writer, "\n{:<4}{:>4}{:>9}{:>11}", "#PRN", "nf", "azimuth", "elevation")?;
    for i in 0..MAXSAT {
        if qc.satinfo[i].obsnum > 0 || qc.satinfo[i].slipnum > 0 {
            let id = satno2id(i + 1);
            let mut nf = 0;
            for code in qc.satinfo[i].code {
                if code > 0 { nf += 1 }
            }
            let mut azel = [0.0; 2];
            if qc.satinfo[i].azel[0].len() > 0 {
                azel[0] = qc.satinfo[i].azel[0][0];
                azel[1] = qc.satinfo[i].azel[1][0];
            }
            writeln!(writer, " {}{:>4}{:>9.2}{:>11.2}", id, nf, azel[0], azel[1])?;
        }
    }
    Ok(())
}

// output signal information[real-time]
fn out_signal_rt<W: Write>(
    sysinfo: &Vec<SysInfo>,
    siginfo: &Vec<Vec<SigInfo>>,
    navsys: usize,
    writer: &mut W,
) -> Result<(), Error> {
    write!(writer, "\n{:<4}{:>5}{:>5}", "#SYS", "Sig", "SNR")?;
    for i in 0..6 {
        let sysname = namebysyslog(i);
        let sys = pow(2, i);
        for j in 0..(NFREQ + NEXOBS) {
            if (sys & navsys) != 0 && sysinfo[i].code[j] != 0 {
                write!(writer, "\n {}{:>2}S{:<2}",
                       sysname, "", OBS_CODES[sysinfo[i].code[j] as usize])?;
                write!(writer, "{:5.1}", siginfo[i][j].snr)?;
            }
        }
    }

    Ok(())
}

// ----------------------------
// qc function
// ----------------------------
/* Statistics while traversing obs */

/// Quality check (QC) for observation (obs) data
///
/// This function performs quality checks on observation data, including time synchronization,
/// data initialization, single-point positioning (SPP), cycle slip detection, and multipath
/// calculation. It prepares the statistics for subsequent processing and analysis.
///
/// # Arguments
/// * `qc` - A mutable reference to the QC struct for storing quality information
/// * `rtk` - A mutable reference to the Rtk struct for storing position information
/// * `obs` - A reference to a vector of Obs structs containing observation data
/// * `nav` - A reference to the Nav struct containing navigation data
pub fn qc_obs(qc: &mut QC, rtk: &mut Rtk, obs: &Vec<Obs>, nav: &Nav) {
    let n = obs.len();
    /* time init and count */
    if set_qctime(qc, obs[0].time)==false { return }
    /* Cal auto sample(ti) */
    cal_sample(&mut qc.ati, obs[0].time);
    /* epoch information init*/
    epoch_init(&mut qc.satinfo, obs, n);
    /* SPP */
    if nav.n + nav.ng > 0 {
        let _ = pntpos(
            obs,
            n,
            &mut nav.clone(),
            &rtk.opt,
            &mut rtk.sol,
            None,
            Some(&mut rtk.ssat),
        );
    }
    /* Cal the neu of each spp */
    get_pos(&mut qc.pos, &rtk);

    /* Count Exp/Hav data for follow calculation */
    count_exp_hav(
        &mut qc.sysinfo,
        &mut qc.satinfo,
        &mut qc.siginfo,
        rtk,
        obs,
        &mut nav.clone(),
        &qc.satlist,
        qc.satlist.len(),
        n,
        LIMIT_EL,
    );
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(&mut qc.satinfo, obs, n, nav);
    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(&mut qc.satinfo, obs, n, nav);
    /* Cal MP by sliding window, std are stored to cal rms */
    get_multipath(&mut qc.mp, &qc.satinfo, nav, obs, n);

    getepdata(
        &mut qc.sysinfo,
        &mut qc.satinfo,
        &mut qc.siginfo,
        rtk,
        obs,
        qc.ati.ti,
        &mut qc.totinfo.minele,
        n,
    );
    /* Count obsnum, slipnum and  starting point to end point of ts */
    osfcount(&mut qc.satinfo, obs, n, qc.ati.ti as i64);

    /* Only the satellites in which the epoch appears are involved */
    for i in 0..n {
        /* store last epoch time */
        qc.satinfo[obs[i].sat - 1].pretime = obs[i].time;
        /* reset */
        qc.satinfo[obs[i].sat - 1].verify = false;
    }
}

/// Performs quality check calculations.
///
/// This function carries out a series of calculations to assess the quality of data received from satellites,
/// including signal residuals, satellite performance, system integrity, and overall observation quality statistics.
///
/// # Arguments
/// - `qc`: A mutable reference to a `QC` struct, which contains all the quality check data and results.
/// - `navsys`: An `usize` value representing the navigation system identifier, used to specify the navigation system.
pub fn qc_cal(qc: &mut QC, navsys: usize) {
    /* Cal rms of mp-std */
    for i in 0..MAXSAT {
        for j in 0..NFREOBS {
            if qc.mp.sat_std[i][j].len() > 0 {
                qc.satinfo[i].sat_mp[j] = rms(&qc.mp.sat_std[i][j]);
            }
        }
    }
    /* Cal spp(neu) std */
    cal_sppneu(&mut qc.pos, &mut qc.totinfo.enu_std);
    /* Summaize Exp/Hav data to cal integrity */
    cal_symbol(
        &mut qc.totinfo,
        &mut qc.sysinfo,
        &mut qc.satinfo,
        &mut qc.siginfo,
        navsys,
    );
    /* Count sys mp */
    sum_mp(&mut qc.sysinfo, &mut qc.satinfo);
    /* Cal total integrity and oslps in symbol */
    qc.totinfo.ep_rt = qc.epnum as f32
        / ((qc.totinfo.lastepoch.time - qc.totinfo.firstepoch.time) as f32 / qc.ati.ti as f32 + 1.0)
        * 100.0;
    qc.totinfo.sat_rt = qc.totinfo.hav.totsatobs as f32 / qc.totinfo.exp.totsatobs as f32 * 100.0;
    qc.totinfo.sat_el_rt = qc.totinfo.hav.totsatobs_el as f32 / qc.totinfo.exp.totsatobs_el as f32 * 100.0;
    qc.totinfo.obs_rt = qc.totinfo.hav.totobs as f32 / qc.totinfo.exp.totobs as f32 * 100.0;
    qc.totinfo.obs_el_rt = qc.totinfo.hav.totobs_el as f32 / qc.totinfo.exp.totobs_el as f32 * 100.0;
    qc.totinfo.oslps = if qc.totinfo.totslip > 0 {
        qc.totinfo.totobs / qc.totinfo.totslip
    } else {
        qc.totinfo.totobs
    };
}

/// Outputs quality check results to a writer.
///
/// This function outputs the quality check results to a writer, including total summary, system summary, satellite summary,
/// signal summary, sky map, and position information. It also handles the case where navigation data is not available.
///
/// # Arguments
/// - `qc`: A reference to a `QC` struct containing the quality check results.
/// - `navsys`: An `usize` value representing the navigation system identifier, used to specify the navigation system.
/// - `q`: An `i32` value indicating the level of output required.
pub fn qc_out<W: Write>(
    qc: &QC,
    navsys: usize,
    q: i32,
    mut nav_flag: bool,
    writer: &mut W,
) {
    if qc.totinfo.obs_el_rt.is_nan() { nav_flag = false }
    out_totsum(&qc.totinfo, nav_flag, qc.ati.ti, LIMIT_EL, writer).expect("Error in total information out!");
    out_syssum(&qc.sysinfo, writer).expect("Error in sys information out");
    if q != 0 {
        if nav_flag {
            out_sat(&mut qc.satinfo.clone(), LIMIT_EL, writer).expect("Error in sat information out!");
            out_signal(&qc.sysinfo, &qc.siginfo, navsys, LIMIT_EL, writer).expect("Error in sig information out");
        }
        if q > 1 {
            if nav_flag {
                out_skymap(&qc.satinfo, qc.totinfo.firstepoch, qc.ati.ti, writer).expect("Error in skymap information out!");
                out_pos(&qc.pos, qc.totinfo.firstepoch, qc.ati.ti, writer).expect("Error in pdop information out!");
            }
            out_ts(&qc.satinfo, qc.totinfo.firstepoch, qc.ati.ti, writer).expect("Error in ts information out!");
        }
    }
}

/// Performs quality check calculations for real-time observations.
///
/// This function performs quality check calculations for real-time observations, including initializing time,
/// counting epoch information, calculating sample time, and performing SPP calculations. It also handles the case
/// where navigation data is not available.
///
/// # Arguments
/// - `qc`: A mutable reference to a `QC` struct, which contains all the quality check data and results.
/// - `rtk`: A mutable reference to a `Rtk` struct, which contains all thereal-time tracking data and results.
pub fn qc_obs_rt(qc: &mut QC, rtk: &mut Rtk, _popt: &PrcOpt, obs: &Vec<Obs>, nav: &Nav) {
    let n = obs.len();
    /* time init and count */
    if set_qctime(qc, obs[0].time)==false { return }
    /* Cal auto sample(ti) */
    cal_sample(&mut qc.ati, obs[0].time);
    /* epoch information init*/
    epoch_init_rt(&mut qc.satinfo, &mut qc.sysinfo, &mut qc.siginfo, &mut qc.mp, obs, n);
    /* SPP */
    if nav.n + nav.ng > 0 {
        let _ = pntpos(
            obs,
            n,
            &mut nav.clone(),
            &rtk.opt,
            &mut rtk.sol,
            None,
            Some(&mut rtk.ssat),
        );
    }

    /* Count Exp/Hav data for follow calculation */
    getepdata(
        &mut qc.sysinfo,
        &mut qc.satinfo,
        &mut qc.siginfo,
        rtk,
        obs,
        qc.ati.ti,
        &mut qc.totinfo.minele,
        n,
    );
    /* Count Exp/Hav data for follow calculation */
    count_exp_hav(
        &mut qc.sysinfo,
        &mut qc.satinfo,
        &mut qc.siginfo,
        rtk,
        obs,
        &mut nav.clone(),
        &qc.satlist,
        qc.satlist.len(),
        n,
        crate::qc::qc_fun::LIMIT_EL,
    );
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(&mut qc.satinfo, obs, n, nav);
    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(&mut qc.satinfo, obs, n, nav);
    /* Cal MP by sliding window, std are stored to cal rms */
    get_multipath(&mut qc.mp, &qc.satinfo, nav, obs, n);


    /* Only the satellites in which the epoch appears are involved */
    for i in 0..n {
        let satidx = obs[i].sat - 1;
        /* store last epoch time */
        qc.satinfo[satidx].pretime = obs[i].time;
        /* reset */
        qc.satinfo[satidx].verify = false;
        /* Count obsnum, slipnum */
        for j in 0..NFREQ + NEXOBS {
            /* slip occur */
            if qc.satinfo[satidx].slip[j] != 0 {
                qc.satinfo[satidx].slipnum += 1;
            }
            /* Count obsnum */
            if obs[i].p[j] != 0.0 {
                qc.satinfo[satidx].obsnum += 1;
            }
        }
    }
}

/// Calculates quality check results for real-time observations.
///
/// This function calculates quality check results for real-time observations, including calculating total summary,
/// system summary, satellite summary, signal summary, sky map, and position information. It also handles the case
/// where navigation data is not available.
///
/// # Arguments
/// - `qc`: A mutable reference to a `QC` struct, which contains all the quality check data and results.
/// - `navsys`: An `usize` value representing the navigation system identifier, usedto specify the navigation system.
pub fn qc_cal_rt(qc: &mut QC, navsys: usize) {
    let mut sys_count = [0.0; 6];
    /* Cal rms of mp-std for sat_mp */
    for i in 0..MAXSAT {
        for j in 0..NFREOBS {
            if qc.mp.sat_std[i][j].len() > 0 {
                qc.satinfo[i].sat_mp[j] = qc.mp.sat_std[i][j][0];
            }
        }
    }

    /* Cal sys_mp*/
    for i in 0..MAXSAT {
        for j in 0..NFREOBS {
            if qc.satinfo[i].sat_mp[j] == 0.0 { continue; }
            let syslog = satsys(i + 1).ilog2() as usize;
            qc.sysinfo[syslog].totmp[j] = (qc.sysinfo[syslog].totmp[j] * sys_count[syslog]
                + qc.satinfo[i].sat_mp[j]) / (sys_count[syslog] + 1.0);
            sys_count[syslog] += 1.0;
        }
    }


    for i in 0..MAXSAT {
        let sys = satsys(i + 1);
        let syslog = (sys as f64).log2() as usize;

        // Count obs/slip
        if qc.satinfo[i].slipnum > 0 || qc.satinfo[i].obsnum > 0 {
            qc.sysinfo[syslog].satnum += 1;
            qc.sysinfo[syslog].slipnum += qc.satinfo[i].slipnum;
            qc.totinfo.totobs += qc.satinfo[i].obsnum;
            qc.totinfo.totslip += qc.satinfo[i].slipnum;
        }
    }

    for i in 0..6 {
        // Count sys obs-Exp/Hav
        for j in 0..(NFREQ + NEXOBS) {
            if qc.siginfo[i][j].exp_obs > 0 {
                qc.sysinfo[i].exp_obs += qc.siginfo[i][j].exp_obs;
                qc.sysinfo[i].hav_obs += qc.siginfo[i][j].hav_obs;
            }
        }
        // Cal sys slip
        if qc.sysinfo[i].exp_obs > 0 {
            qc.sysinfo[i].oslps = if qc.sysinfo[i].slipnum > 0 {
                qc.sysinfo[i].hav_obs / qc.sysinfo[i].slipnum
            } else {
                qc.sysinfo[i].hav_obs
            };
        }
        // To be total obs data
        let sys = 1 << i;
        for j in 0..(NFREQ + NEXOBS) {
            if (sys & navsys) != 0 && qc.sysinfo[i].code[j] != 0 {
                qc.totinfo.hav.totobs += qc.siginfo[i][j].hav_obs;
                qc.totinfo.hav.totobs_el += qc.siginfo[i][j].hav_obs_el;
                qc.totinfo.exp.totobs += qc.siginfo[i][j].exp_obs;
                qc.totinfo.exp.totobs_el += qc.siginfo[i][j].exp_obs_el;
            }
        }
    }

    /* Cal total integrity and oslps in symbol */
    qc.totinfo.ep_rt = if qc.ati.ti > 0 {
        qc.epnum as f32
            / ((qc.totinfo.lastepoch.time - qc.totinfo.firstepoch.time) as f32 / qc.ati.ti as f32 + 1.0)
            * 100.0
    } else { 0.0 };
    qc.totinfo.obs_rt = qc.totinfo.hav.totobs as f32 / qc.totinfo.exp.totobs as f32 * 100.0;
    qc.totinfo.obs_el_rt = qc.totinfo.hav.totobs_el as f32 / qc.totinfo.exp.totobs_el as f32 * 100.0;
    qc.totinfo.oslps = if qc.totinfo.totslip > 0 {
        qc.totinfo.totobs / qc.totinfo.totslip
    } else {
        qc.totinfo.totobs
    };
    qc.totinfo.sample = qc.ati.ti;
}

/// Outputs real-time quality check results to a writer.
///
/// This function outputs real-time quality check results to a writer, including total summary, system summary,
/// satellite summary, signal summary, sky map, and position information. It also handles the case where navigation
/// data is not available.
///
/// # Arguments
/// - `qc`: A reference to a `QC` struct, which contains all the quality check data and results.
/// - `navsys`: An `usize` value representing the navigation system identifier, used to specify the navigation system.
pub fn qc_out_rt<W: Write>(
    qc: &QC,
    navsys: usize,
    mut nav_flag: bool,
    rtk: &Rtk,
    writer: &mut W,
) {
    if qc.totinfo.obs_el_rt.is_nan() { nav_flag = false }
    out_totsum_rt(&qc.totinfo, nav_flag, qc.ati.ti, LIMIT_EL, &rtk.sol.rr, writer).expect("Error in total information out!");
    out_syssum(&qc.sysinfo, writer).expect("Error in sys information out");
    out_sat_rt(&qc, writer).expect("Error in realtime sat information out");
    out_signal_rt(&qc.sysinfo, &qc.siginfo, navsys, writer).expect("Error in realtime sig information out");
}

/// Exports real-time quality check results to the specified output type.
///
/// # Arguments
///
/// * `site` - A string slice that holds the site identifier.
/// * `qc` - A mutable reference to a `QC` struct for quality check data.
/// * `rtk` - A mutable reference to an `Rtk` struct for position information.
/// * `popt` - A reference to a `PrcOpt` struct for processing options.
/// * `obss` - A reference to an `Obss` struct containing observation data.
/// * `navs` - A reference to a `Nav` struct containing navigation data.
/// * `outtype` - A reference to an `Outtype` enum specifying the output type.
///
/// # Output Types
///
/// The function supports three types of outputs:
/// - `File`: Writes QC data to a file.
/// - `Redis`: Sends QC data to a Redis server.
/// - `All`: Writes QC data to both a file and a Redis server.
pub fn rt_qc_export(
    site: &str,
    qc: &mut QC,
    rtk: &mut Rtk,
    popt: &PrcOpt,
    obss: &Obss,
    navs: &Nav,
    outtype: &Outtype,
) {
    let mut outflag = true;
    let navflag = if navs.n == 0 || navs.ng == 0 { false } else { true };
    qc_obs_rt(qc, rtk, &popt, &obss.obs, navs);
    qc_cal_rt(qc, popt.navsys);
    match outtype {
        Outtype::File(qcfile) => {
            if !&qcfile.is_empty() {
                let file = File::create(qcfile).expect("Fail to create output file!");
                let mut writer = BufWriter::new(file);
                qc_out_rt(qc, popt.navsys, navflag, &rtk, &mut writer);
            }
        }
        Outtype::Redis(redis_params) => {
            let r = qc2redis(site, qc, obss, redis_params.ip.clone(), redis_params.port.clone(), redis_params.auth.clone());
            if r.is_err() { outflag = false }
        }
        Outtype::All(qcfile, redis_params) => {
            if !&qcfile.is_empty() {
                let file = File::create(qcfile).expect("Fail to create output file!");
                let mut writer = BufWriter::new(file);
                qc_out_rt(qc, popt.navsys, navflag, &rtk, &mut writer);
            }
            let r = qc2redis(site, qc, obss, redis_params.ip.clone(), redis_params.port.clone(), redis_params.auth.clone());
            if r.is_err() { outflag = false }
        }
    }

    if outflag { print!("\rProcessing {} at {}", site, time2str(obss.obs[0].time)) }
    io::stdout().flush().unwrap();

    // reset
    for i in 0..MAXSAT {
        qc.satinfo[i].code = [0; 6];
    }
}

/// Performs quality check on GNSS observation data and writes the results.
///
/// This function processes GNSS observation data and outputs the quality check results.
///
/// # Arguments
/// - `writer`: A mutable reference to a type that implements the `Write` trait, used for writing the output.
/// - `popt`: A reference to `PrcOpt`, which contains processing options.
/// - `detail_level`: An integer representing the detail level of the output.
/// - `obss`: A reference to `Obss`, which contains the observation data.
/// - `navs`: A reference to `Nav`, which contains the navigation data.
///
/// # Returns
/// - `Result<(), String>`: Returns `Ok(())` if the quality control process is successful, or an error message wrapped in `Err` if an error occurs.
pub fn qc_rnx<W: Write>(
    writer: &mut W,
    popt: &PrcOpt,
    detail_level: i32,
    obss: &Obss,
    navs: &Nav,
) -> Result<(), String> {
    let mut rtk = Rtk::init(&popt);
    let mut qc = QC::init();
    let mut nobs;
    let mut count = 0;
    let mut obs: Vec<Obs> = Vec::new();
    let mut nav_flag = true;

    if navs.n == 0 && navs.ng == 0 {
        nav_flag = false;
    }

    // set progressbar
    let pb = ProgressBar::new(obss.n as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}")
        .unwrap().progress_chars("#>-"));

    // store information of code
    for i in 0..obss.n {
        let mut sat_exist = false;
        let sys = satsys(obss.obs[i].sat);
        let syslog = sys.ilog2() as usize;
        for j in 0..qc.satlist.len() {
            if qc.satlist[j] == obss.obs[i].sat {
                sat_exist = true;
                break;
            }
        }
        if !sat_exist {
            qc.satlist.push(obss.obs[i].sat);
        }
        for j in 0..NFREOBS {
            if obss.obs[i].code[j] != 0 && qc.sysinfo[syslog].code[j] == 0 {
                qc.sysinfo[syslog].code[j] = obss.obs[i].code[j];
                qc.sysinfo[syslog].codenum += 1;
            }
            if qc.satinfo[obss.obs[i].sat - 1].code[j] == 0 {
                qc.satinfo[obss.obs[i].sat - 1].code[j] = obss.obs[i].code[j];
            }
        }
    }

    // qc obs by epoch
    nobs = inputobs(&obss, &mut obs);
    while nobs >= 0 {
        qc_obs(&mut qc, &mut rtk, &obs, &navs.clone());

        pb.inc(obs.len() as u64);
        count += obs.len();
        let message = format!("QC obs processing {:.2}%", count as f32 / obss.n as f32 * 100.0);
        pb.set_message(message.to_string());

        obs = Vec::new();
        nobs = inputobs(&obss, &mut obs);
    }

    if qc.ati.ti <= 0 {
        return Err(String::from("can not cal sample"));
    }

    qc_cal(&mut qc, popt.navsys);
    qc_out(&qc, popt.navsys, detail_level, nav_flag, writer);
    pb.finish_with_message(String::from("QC completed"));
    return Ok(());
}

/// Performs quality check (QC) on observation data and outputs the results.
///
/// This function reads observation data from input files, performs quality check processing,
/// and writes the results to a specified output file.
///
/// # Arguments
/// - `opt` - A mutable reference to the QcOpt struct containing processing options
/// - `ifile` - A reference to a vector of strings representing input file paths
/// - `ofile` - A string slice representing the output file path
pub fn obsqc(
    opt: &mut QcOpt,
    ifile: &Vec<String>,
    ofile: &str,
) {
    let mut obss: Obss = Obss::new();
    let mut navs = Nav::new();
    let mut sta = Sta::default();
    let mut nepoch = 0;
    let index = [0; MAXINFILE];

    // read rnx_file
    let read_res = readobsnav(
        opt.ts,
        opt.te,
        opt.ti,
        ifile.clone(),
        index,
        &mut obss,
        &mut navs,
        &mut sta,
        &mut nepoch,
    );

    match read_res {
        Ok(true) => {
            let file = File::create(ofile).expect("Fail to create output file!");
            let mut writer = BufWriter::new(file);

            //output files I/O
            for file in ifile{
                writeln!(writer, "inp file:{}", file).expect("error ifile");
            }
            writeln!(writer, "out file:{}", ofile).expect("error ofile");

            match qc_rnx(&mut writer, &opt.popt, opt.detail, &obss, &navs) {
                Ok(_) => {}
                Err(e) => eprintln!("Error is {}", e)
            }
        }
        Ok(false) => eprintln!("read fail"),
        Err(e) => eprintln!("read rnx error is {}", e)
    }
}

/// Performs quality check (QC) on RTCM data.
///
/// This function reads RTCM data from the input file, processes it to perform quality checks,
/// and writes the results to the output file.
///
/// # Arguments
///
/// - `opt` - A mutable reference to a `QcOpt` struct containing configuration options.
/// - `ifile` - A string slice representing the path to the input RTCM data file.
/// - `ofile` - A string slice representing the path to the output QC results file.
pub fn rtcmqc(
    opt: &mut QcOpt,
    ifile: &str,
    ofile: &str,
) {
    let rtcm_path = Path::new(ifile);
    let rtcmfile = File::open(rtcm_path).expect("Fail to open rtcm file!");
    let mut reader = BufReader::new(rtcmfile);
    let qcfile = File::create(ofile).expect("Fail to create qc file!");
    let mut writer = BufWriter::new(qcfile);

    let popt = opt.popt;
    let mut rtk = Rtk::init(&popt);
    let mut navs = Nav::new();
    let mut qc = QC::init();
    let mut loop_num = 0;
    let mut count = 0;

    let decoder;
    decoder = Decoder::new();

    let _ = scan_rtcm(&decoder, &ifile, &mut navs, &mut qc, popt.navsys, opt.tr, &mut loop_num);
    let nav_flag = if navs.n == 0 || navs.ng == 0 { false } else { true };

    let pb = ProgressBar::new(loop_num);

    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}")
        .unwrap().progress_chars("#>-"));

    let mut trt = opt.tr;

    loop {
        let mut obss = Obss::new();
        match decode_epoch(&decoder, &mut reader, &mut obss, None, trt) {
            Err(e) => match e {
                IoError => break,
                EndofFile => break,
                _ => continue
            }
            Ok(false) => continue,
            Ok(true) => {
                let mut obsforqc = Vec::new();
                for i in 0..obss.n {
                    if (satsys(obss.obs[i].sat) & popt.navsys) != 0 && popt.exsats[obss.obs[i].sat - 1] != 1 {
                        obsforqc.push(obss.obs[i].clone());
                    }
                }
                if obsforqc.len() == 0 {
                    continue;
                }

                qc_obs(&mut qc, &mut rtk, &obsforqc, &navs);


                trt = obsforqc[0].time;

                pb.inc(1);
                count += 1;
                let message = format!("QC obs processing: {:.2}%", count as f32 / loop_num as f32 * 100.0);
                pb.set_message(message.to_string());
            }
        }
    }

    qc_cal(&mut qc, popt.navsys);
    pb.set_message(String::from("Cal qc data"));

    // output files I/O
    writeln!(writer, "inp file:{}", ifile).expect("error ifile");
    writeln!(writer, "out file:{}", ofile).expect("error ofile");

    qc_out(&qc, popt.navsys, opt.detail, nav_flag, &mut writer);

    pb.set_message(String::from("Output qc data"));
    pb.finish_with_message(String::from("QC completed"));
}

/// Performs quality check (QC) on RTCM real-time data stream.
/// This function receives an RTCM data stream, decodes it, processes observation and navigation data,
/// and then performs quality check checks.
///
/// # Arguments
///
/// - `decoder` - A reference to the RTCM decoder used to decode the incoming data stream.
/// - `site` - A string slice that holds the site identifier.
/// - `data` - A vector of bytes containing the RTCM data to be processed.
/// - `epochtime` - A mutable reference to an i64 integer representing the epoch time.
/// - `obss` - A mutable reference to an `Obss` struct for storing observation data.
/// - `navs` - A mutable reference to a `Nav` struct for storing navigation data.
/// - `popt` - A reference to a `PrcOpt` struct containing processing options.
/// - `rtk` - A mutable reference to an `Rtk` struct for RTK processing.
/// - `qc` - A mutable reference to a `QC` struct for quality check.
/// - `supple_nav` - An optional reference to a `Nav` struct providing supplementary navigation data.
/// - `outtype` - Specifies the output type for the quality control export.
pub fn rt_rtcm2qc(
    decoder: &Decoder,
    site: &str,
    data: Vec<u8>,
    epochtime: &mut i64,
    obss: &mut Obss,
    navs: &mut Nav,
    popt: &PrcOpt,
    rtk: &mut Rtk,
    qc: &mut QC,
    supple_nav: &Option<Nav>,
    outtype: Outtype,
) {
    setup_panic_hook();

    // Direct replacement of stored data
    PANIC_DATA.with(|storage| {
        let mut data_ref = storage.borrow_mut();
        // Replace the contents of Vec<u8> in Option, storing new data each time.
        *data_ref = Some(data.clone());  // Directly overwrite previously stored Vec<u8>
    });

    let tr = GTime {
        time: timeget(),
        sec: 0.0,
    };
    update_rtnav(navs, tr);
    match decoder.decode_stream(&data) {
        None => {}
        Some(rtcm) => {
            if rtcm.eph.sat > 0 {
                let _ = navs.add_eph(rtcm.eph);
                if navs.eph.len() != navs.n {
                    eprintln!("nav eph long error")
                }
                uniqnav(navs);
            } else if rtcm.geph.sat > 0 {
                let _ = navs.add_geph(rtcm.geph);
                uniqnav(navs);
            } else if rtcm.obss.n > 0 {
                if *epochtime == 0 {
                    *epochtime = rtcm.time.time;
                }
                if *epochtime != rtcm.time.time {
                    let nav = if navs.n + navs.ng == 0 && supple_nav.is_some()
                    { supple_nav.clone().unwrap() } else { navs.clone() };
                    rt_qc_export(site, qc, rtk, popt, obss, &nav, &outtype);

                    *epochtime = 0;
                    *obss = Obss::new();
                }
                for i in 0..rtcm.obss.n {
                    let obs = rtcm.obss.obs[i].clone();
                    let sys = satsys(obs.sat);
                    let syslog = sys.ilog2() as usize;
                    if !qc.satlist.contains(&obs.sat) {
                        qc.satlist.push(obs.sat);
                    }
                    for j in 0..NFREOBS {
                        if obs.code[j] > 0 && qc.sysinfo[syslog].code[j] == 0 {
                            qc.sysinfo[syslog].code[j] = obs.code[j];
                            qc.sysinfo[syslog].codenum += 1;
                        }
                        if qc.satinfo[obs.sat - 1].code[j] == 0 &&
                            (obs.l[j] != 0.0 || obs.p[j] != 0.0) {
                            qc.satinfo[obs.sat - 1].code[j] = obs.code[j];
                        }
                    }
                    let _ = obss.add_obs_data(&rtcm.obss.obs[i]);
                }
                if rtcm.obsflag {
                    let nav = if navs.n + navs.ng == 0 && supple_nav.is_some()
                    { supple_nav.clone().unwrap() } else { navs.clone() };
                    rt_qc_export(site, qc, rtk, popt, obss, &nav, &outtype);
                    *epochtime = 0;
                    *obss = Obss::new();
                }
            }
        }
    }
}


