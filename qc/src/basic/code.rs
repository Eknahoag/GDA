use crate::basic::var::*;

use super::sat::{getprn, satsys};

/// GPS obs code to frequency
pub fn code2freq_gps(code: u8, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    match obs.chars().next() {
        Some('1') => {
            *freq = FREQ1;
            0
        }
        Some('2') => {
            *freq = FREQ2;
            1
        }
        Some('5') => {
            *freq = FREQ5;
            2
        }
        _ => -1,
    }
}

/// GLONASS obs code to frequency
pub fn code2freq_glo(code: u8, fcn: i32, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    if fcn < -7 || fcn > 6 {
        return -1;
    }

    match obs.chars().next() {
        Some('1') => {
            *freq = FREQ1_GLO + DFRQ1_GLO * fcn as f64;
            0
        }
        Some('2') => {
            *freq = FREQ2_GLO + DFRQ2_GLO * fcn as f64;
            1
        }
        Some('3') => {
            *freq = FREQ3_GLO;
            2
        }
        Some('4') => {
            *freq = FREQ1A_GLO;
            0
        }
        Some('6') => {
            *freq = FREQ2A_GLO;
            1
        }
        _ => -1,
    }
}

/// Galileo obs code to frequency
pub fn code2freq_gal(code: u8, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    match obs.chars().next() {
        Some('1') => {
            *freq = FREQ1;
            0
        }
        Some('7') => {
            *freq = FREQ7;
            1
        }
        Some('5') => {
            *freq = FREQ5;
            2
        }
        Some('6') => {
            *freq = FREQ6;
            3
        }
        Some('8') => {
            *freq = FREQ8;
            4
        }
        _ => -1,
    }
}

/// QZSS obs code to frequency
pub fn code2freq_qzs(code: u8, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    match obs.chars().next() {
        Some('1') => {
            *freq = FREQ1;
            0
        }
        Some('2') => {
            *freq = FREQ2;
            1
        }
        Some('5') => {
            *freq = FREQ5;
            2
        }
        Some('6') => {
            *freq = FREQ6;
            3
        }
        _ => -1,
    }
}

/// BeiDou obs code to frequency
pub fn code2freq_bds(code: u8, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    match obs.chars().next() {
        Some('1') => {
            *freq = FREQ1;
            0
        }
        Some('2') => {
            *freq = FREQ1_CMP;
            0
        }
        Some('7') => {
            *freq = FREQ2_CMP;
            1
        }
        Some('5') => {
            *freq = FREQ5;
            2
        }
        Some('6') => {
            *freq = FREQ3_CMP;
            3
        }
        Some('8') => {
            *freq = FREQ8;
            4
        }
        _ => -1,
    }
}

/// NavIC obs code to frequency
pub fn code2freq_irn(code: u8, freq: &mut f64) -> i32 {
    let obs = code2obs(code);

    match obs.chars().next() {
        Some('5') => {
            *freq = FREQ5;
            0
        }
        Some('9') => {
            *freq = FREQ9;
            1
        }
        _ => -1,
    }
}

/// obs type string to obs code
///
/// # Arguments
/// - `obs`: obs code string ("1C","1P","1Y",...)
///
/// # Returns
/// obs code (CODE_???)
///
/// # Notes
/// obs codes are based on RINEX 3.04
pub fn obs2code(obs: &str) -> u8 {
    for (i, &code) in OBS_CODES.iter().enumerate().skip(1) {
        if code == obs {
            return i as u8;
        }
    }
    CODE_NONE
}

/// obs code to obs code string
///
/// # Arguments
/// - 'code': obs code (CODE_???)
///
/// # Returns
/// obs code string ("1C","1P","1P",...)
///
/// # Notes
/// obs codes are based on RINEX 3.04
pub fn code2obs(code: u8) -> &'static str {
    if code <= CODE_NONE || code > MAXCODE {
        ""
    } else {
        OBS_CODES[code as usize]
    }
}

/// system and obs code to frequency index
///
/// # Arguments
/// - `sys`: system (SYS_???)
/// - `code`: obs code (CODE_???)
///
/// # Returns
/// frequency index (-1:error)
/// 0     1     2     3     4
/// --------------------------------------
/// GPS       L1    L2    L5     -     -
/// GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
/// Galileo   E1    E5b   E5a   E6   E5ab
/// QZSS      L1    L2    L5    L6     -
/// BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
/// NavIC     L5     S     -     -     -
pub fn code2idx(sys: usize, code: u8) -> i32 {
    let mut freq = 0.0;

    match sys {
        SYS_GPS => code2freq_gps(code, &mut freq),
        SYS_GLO => code2freq_glo(code, 0, &mut freq),
        SYS_GAL => code2freq_gal(code, &mut freq),
        SYS_QZS => code2freq_qzs(code, &mut freq),
        SYS_CMP => code2freq_bds(code, &mut freq),
        SYS_IRN => code2freq_irn(code, &mut freq),
        _ => -1,
    }
}

/// get code priority for multiple codes in a frequency
///
/// # Arguments
/// - `sys`: system (SYS_???)
/// - `code`: obs code (CODE_???)
/// - `opt`: code options (NULL:no option)
///
/// # Returns
/// priority (15:highest-1:lowest,0:error)
pub fn getcodepri(sys: usize, code: u8, opt: Option<&str>) -> i32 {
    let (i, optstr) = match sys {
        SYS_GPS => (0, "-GL"),
        SYS_GLO => (1, "-RL"),
        SYS_GAL => (2, "-EL"),
        SYS_QZS => (3, "-JL"),
        SYS_CMP => (4, "-CL"),
        SYS_IRN => (5, "-IL"),
        _ => return 0,
    };

    let j = code2idx(sys, code);
    if j < 0 {
        return 0;
    }

    let obs = code2obs(code);

    if let Some(opt) = opt {
        if let Some(pos) = opt.find(optstr) {
            let start = pos + optstr.len();
            if start + 1 < opt.len() {
                let str_slice = &opt[start..start + 2];
                if str_slice == &obs[0..2] {
                    return 15;
                }
            }
        }
    }

    // search code priority
    if let Some(p) = CODEPRIS[i][j as usize].find(obs.chars().nth(1).unwrap()) {
        return 14 - p as i32;
    }

    return 0;
}

/// system and obs code to frequency
///
/// # Arguments
/// - `sys`: system (SYS_???)
/// - `code`: obs code (CODE_???)
/// - `fcn`: frequency channel number for GLONASS
///
/// # Returns
/// carrier frequency (Hz) (0.0: error)
pub fn code2freq(sys: usize, code: u8, fcn: i32) -> f64 {
    let mut freq = 0.0;

    match sys {
        SYS_GPS => code2freq_gps(code, &mut freq),
        SYS_GLO => code2freq_glo(code, fcn, &mut freq),
        SYS_GAL => code2freq_gal(code, &mut freq),
        SYS_QZS => code2freq_qzs(code, &mut freq),
        SYS_CMP => code2freq_bds(code, &mut freq),
        SYS_IRN => code2freq_irn(code, &mut freq),
        _ => -1,
    };

    return freq;
}

/// satellite and obs code to frequency
///
/// # Arguments
/// - `sat`: satellite number (1..MAXSAT)
/// - `code`: obs code (CODE_???)
/// - `nav`: navigation data for GLONASS (NONE:not used)
pub fn sat2freq(sat: usize, code: u8, nav: Option<&Nav>) -> f64 {
    let mut fcn = 0;
    let sys = satsys(sat);
    let prn = getprn(sat);
    if sys == SYS_GLO {
        if nav.is_none() || nav.unwrap().ng == 0 {
            fcn = GLO_FCN[sat - 33];
        } else {
            let nav = nav.unwrap();
            let mut found = false;
            for i in 0..nav.ng {
                if nav.geph[i].sat == sat {
                    fcn = nav.geph[i].frq;
                    found = true;
                    break;
                }
            }
            if !found {
                if nav.glo_fcn[prn - 1] > 0 {
                    fcn = nav.glo_fcn[prn - 1] - 8;
                } else {
                    return 0.0;
                }
            }
        }
    }
    return code2freq(sys, code, fcn);
}

/// satellite to PRN
pub fn sat2code(sat:usize)->String{
    let prn=getprn(sat);
    return match satsys(sat) {
        SYS_GPS => format!("G{:02}", prn),
        SYS_GLO => format!("R{:02}", prn),
        SYS_GAL => format!("E{:02}", prn),
        SYS_QZS => format!("J{:02}", prn),
        SYS_CMP => format!("C{:02}", prn),
        SYS_IRN => format!("I{:02}", prn),
        _ => "".to_string()
    }
}