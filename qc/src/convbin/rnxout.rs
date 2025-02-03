use std::cmp::{max, min};
use std::collections::LinkedList;
use std::fs::{create_dir_all, File, remove_file};
use std::io::{self, BufRead, BufReader, BufWriter, Error, Seek, SeekFrom, Write};
use std::path::Path;
use std::string::String;
use indicatif::{ProgressBar, ProgressStyle};
use tokio::io::AsyncReadExt;
use tokio::sync::mpsc;
use crate::basic::code::{code2idx, code2obs, getcodepri, obs2code, sat2code};
use crate::basic::eph::{sisa_value, uravalue};
use crate::basic::func::norm;
use crate::basic::ntrip::conntrip;
use crate::basic::pos::{ecef2enu, ecef2pos};
use crate::basic::sat::{getprn, satid2no, satno2id, satsys};
use crate::basic::time::{gpst2bdt, gpst2utc, time2bdt, time2epoch, time2gpst, time2str, timediff, timestr_rnx};
use crate::basic::var::*;
use crate::convbin::decode::{decode_epoch, DecodeError, Decoder, read_single_rtcm};
use crate::convbin::decode::DecodeError::{EndofFile, IoError};
use crate::qc::qc_fun::cal_sample;
use crate::qc::qc_var::AutoTi;

static SYSCODES: &str = "GREJCI";

static VERCODE: [&str; 6] = [
    "00000000...0.0000000000000..........................................", // GPS
    "00...........0....0..........44.4..........222...................444", // GLO
    "0........0000..........0000000000...000.............................", // GAL
    "2.....22...22..222.....222......2422....................4444........", // QZS
    ".4...4...4.4.....1.......41114..1.....41111............444..44444...", // BDS
    ".........................3......................3333333.............", // IRN
];

fn save_slips(slips: &mut [[u8; NFREOBS]], data: &Vec<Obs>) {
    for i in 0..data.len() {
        for j in 0..NFREOBS {
            if (data[i].lli[j] & LLI_SLIP) != 0 { slips[data[i].sat - 1][j] = 1 }
        }
    }
}
#[allow(dead_code)]
fn rest_slips(slips: &mut [[u8; NFREOBS]], data: &mut Vec<Obs>) {
    for i in 0..data.len() {
        for j in 0..NFREOBS {
            if data[i].l[j] != 0.0 && slips[data[i].sat - 1][j] != 0 {
                data[i].lli[j] |= LLI_SLIP;
                slips[data[i].sat - 1][j] = 0;
            }
        }
    }
}

// output obs-types RINEX ver.2
fn outobstype_ver2<W: Write>(writer: &mut W, opt: &RnxOpt) -> Result<(), Error> {
    const LABEL: &str = "# / TYPES OF OBSERV";

    write!(writer, "{:6}", opt.nobs[0])?;

    for i in 0..opt.nobs[0] as usize {
        if i > 0 && i % 9 == 0 {
            write!(writer, "      ")?;
        }

        let tobs = opt.tobs[0][i].iter().collect::<String>();

        write!(writer, "{:6}", tobs.trim_end_matches('\0'))?;


        if i % 9 == 8 {
            write!(writer, "{:<20}\n", LABEL)?;
        }
    }

    if opt.nobs[0] == 0 || opt.nobs[0] % 9 != 0 {
        write!(writer, "{:width$}{:<20}\n", "", LABEL, width = (9 - opt.nobs[0] as usize % 9) * 6)?;
    }

    Ok(())
}

// output obs-types RINEX ver.3
fn outobstype_ver3<W: Write>(writer: &mut W, opt: &RnxOpt) -> Result<(), Error> {
    const LABEL: &str = "SYS / # / OBS TYPES";

    for i in 0..NUMSYS {
        if 1 << i == 0 || (1 << i & opt.navsys) == 0 || opt.nobs[i] == 0 {
            continue;
        }

        write!(writer, "{}  {:3}", SYSCODES.chars().nth(i).unwrap(), opt.nobs[i])?;

        for j in 0..opt.nobs[i] as usize {
            if j > 0 && j % 13 == 0 {
                write!(writer, "      ")?;
            }

            let mut tobs = opt.tobs[i][j].iter().collect::<String>();

            // BDS B2x -> 1x (3.02), 2x (other)
            if 1 << i == SYS_CMP && opt.rnxver == 302 && tobs.chars().nth(1) == Some('2') {
                tobs.replace_range(1..2, "1");
            }

            write!(writer, " {:3}", tobs.trim_end_matches('\0'))?;

            if j % 13 == 12 {
                write!(writer, "  {:<20}\n", LABEL)?;
            }
        }

        if opt.nobs[i] % 13 != 0 {
            write!(writer, "{:width$}  {:<20}\n", "", LABEL, width = (13 - opt.nobs[i] as usize % 13) * 4)?;
        }
    }

    Ok(())
}

// output RINEX phase shift
fn outrnx_phase_shift<W: Write>(writer: &mut W, opt: &RnxOpt) -> Result<(), Error> {
    const LABEL: &str = "SYS / PHASE SHIFT";
    const REF_CODE: [&'static [u8]; 6] = [
        &[CODE_L1C, CODE_L2P, CODE_L5I],                   // GPS
        &[CODE_L1C, CODE_L4A, CODE_L2C, CODE_L6A, CODE_L3I], // GLO
        &[CODE_L1B, CODE_L5I, CODE_L7I, CODE_L8I, CODE_L6B], // GAL
        &[CODE_L1C, CODE_L2S, CODE_L5I, CODE_L5D, CODE_L6S], // QZS
        &[CODE_L2I, CODE_L1D, CODE_L5D, CODE_L7I, CODE_L7D, CODE_L8D, CODE_L6I], // BDS
        &[CODE_L5A, CODE_L9A],                             // IRN
    ];

    for i in 0..NUMSYS {
        if 1 << i == 0 || (1 << i & opt.navsys) == 0 || opt.nobs[i] == 0 {
            continue;
        }

        for j in 0..opt.nobs[i] as usize {
            if opt.tobs[i][j][0] != 'L' {
                continue;
            }

            let mut obs = opt.tobs[i][j].iter().collect::<String>();

            let mut found_ref_code = false;
            for &ref_code in REF_CODE[i] {
                if obs2code(&obs[1..3]) == ref_code {
                    found_ref_code = true;
                    break;
                }
            }

            // BDS B2x -> 1x (3.02), 2x (other)
            if 1 << i == SYS_CMP && opt.rnxver == 302 && obs.chars().nth(1) == Some('2') {
                obs.replace_range(1..2, "1");
            }

            if found_ref_code {
                writeln!(writer, "{} {:3} {:54}{}", SYSCODES.chars().nth(i).unwrap(), obs.trim_end_matches('\0'), "", LABEL)?;
            } else {
                writeln!(writer, "{} {:3} {:8.5}{:46}{}", SYSCODES.chars().nth(i).unwrap(), obs.trim_end_matches('\0'), opt.shift[i][j], "", LABEL)?;
            }
        }
    }

    Ok(())
}

// output RINEX GLONASS slot/freq
fn outrnx_glo_fcn<W: Write>(writer: &mut W, opt: &RnxOpt, mut navs: Option<&Nav>) -> Result<(), Error> {
    const LABEL: &str = "GLONASS SLOT / FRQ #";
    const GLO_FCN: [i32; 24] = [
        1, -4, 5, 6, 1, -4, 5, 6,
        -2, -7, 0, -1, -2, -7, 0, -1,
        4, -3, 3, 2, 4, -3, 3, 2
    ];
    let mut prn = Vec::new();
    let mut fcn = Vec::new();

    if navs.is_some() && navs.unwrap().ng == 0 { navs = None }

    if opt.navsys & SYS_GLO != 0 {
        if let Some(nav) = navs {
            for i in 0..MAXPRNGLO {
                let sat = satid2no(format!("R{:02}", i + 1).as_str());
                if nav.geph.len() >= i + 1 && nav.geph[i].sat == sat {
                    prn.push(i + 1);
                    fcn.push(nav.geph[i].frq);
                } else if nav.glo_fcn[i] != 0 {
                    prn.push(i + 1);
                    fcn.push(nav.glo_fcn[i] - 8);
                }
            }
        } else {
            for i in 0..GLO_FCN.len() {
                prn.push(i + 1);
                fcn.push(GLO_FCN[i]);
            }
        }
    }

    let n = prn.len();
    for i in 0..if n <= 0 { 1 } else { (n - 1) / 8 + 1 } {
        if i == 0 {
            write!(writer, "{:3}", n)?;
        } else {
            write!(writer, "   ")?;
        }
        for k in 0..8.min(n - i * 8) {
            let j = i * 8 + k;
            write!(writer, " R{:02} {:2}", prn[j], fcn[j])?;
        }
        write!(writer, "{:width$} {}\n", "", LABEL, width = (8 - 8.min(n - i * 8)) * 7)?;
    }

    Ok(())
}

// output RINEX GLONASS code/phase/bias
fn outrnx_glo_bias<W: Write>(writer: &mut W, opt: &RnxOpt) -> Result<(), Error> {
    const LABEL: &str = "GLONASS COD/PHS/BIS";
    const TOBS: [&str; 4] = ["C1C", "C1P", "C2C", "C2P"];

    if opt.navsys & SYS_GLO != 0 {
        writeln!(
            writer,
            " {} {:8.3} {} {:8.3} {} {:8.3} {} {:8.3}{:8}{}",
            TOBS[0], opt.glo_cp_bias[0],
            TOBS[1], opt.glo_cp_bias[1],
            TOBS[2], opt.glo_cp_bias[2],
            TOBS[3], opt.glo_cp_bias[3],
            "", LABEL
        )?;
    } else {
        writeln!(writer, "{:60}{}", "", LABEL)?;
    }

    Ok(())
}

// output RINEX observation header
fn outrnxobsh<W: Write>(
    writer: &mut W,
    opt: &RnxOpt,
    nav: Option<&Nav>,
) -> Result<(), Error> {
    let sys = match opt.navsys {
        SYS_GPS => "G: GPS",
        SYS_GLO => "R: GLONASS",
        SYS_GAL => "E: Galileo",
        SYS_QZS => "J: QZSS",
        SYS_CMP => "C: BeiDou",
        SYS_IRN => "I: IRNSS",
        _ => "M: Mixed",
    };

    let date = timestr_rnx();
    let tsys = "GPS";

    writeln!(writer, "{:9.2}{:<11}{:<20}{:<20}{:<20}", opt.rnxver as f32 / 100.0, "", "OBSERVATION DATA", sys, "RINEX VERSION / TYPE")?;
    writeln!(writer, "{:<20.20}{:<20.20}{:<20.20}{:<20}", opt.prog, opt.runby, date, "PGM / RUN BY / DATE")?;

    for comment in opt.comment.iter() {
        if comment.is_empty() { continue; }
        writeln!(writer, "{:<60.60}{:<20}", comment, "COMMENT")?;
    }

    writeln!(writer, "{:<60.60}{:<20}", opt.marker, "MARKER NAME")?;
    writeln!(writer, "{:<20.20}{:<40.40}{:<20}", opt.markerno, "", "MARKER NUMBER")?;

    if opt.rnxver >= 300 {
        writeln!(writer, "{:<20.20}{:<40.40}{:<20}", opt.markertype, "", "MARKER TYPE")?;
    }

    writeln!(writer, "{:<20.20}{:<40.40}{:<20}", opt.name[0], opt.name[1], "OBSERVER / AGENCY")?;
    writeln!(writer, "{:<20.20}{:<20.20}{:<20.20}{:<20}", opt.rec[0], opt.rec[1], opt.rec[2], "REC # / TYPE / VERS")?;
    writeln!(writer, "{:<20.20}{:<20.20}{:<20.20}{:<20}", opt.ant[0], opt.ant[1], "", "ANT # / TYPE")?;

    let pos: Vec<f64> = opt.apppos.iter().map(|&x| if x.abs() < 1e8 { x } else { 0.0 }).collect();
    let del: Vec<f64> = opt.antdel.iter().map(|&x| if x.abs() < 1e8 { x } else { 0.0 }).collect();

    writeln!(writer, "{:14.4}{:14.4}{:14.4}{:<18}{:<20}", pos[0], pos[1], pos[2], "", "APPROX POSITION XYZ")?;
    writeln!(writer, "{:14.4}{:14.4}{:14.4}{:<18}{:<20}", del[0], del[1], del[2], "", "ANTENNA: DELTA H/E/N")?;

    if opt.rnxver <= 299 {
        writeln!(writer, "{:6}{:6}{:<48}{:<20}", 1, 1, "", "WAVELENGTH FACT L1/2")?;
        outobstype_ver2(writer, opt)?;
    } else {
        outobstype_ver3(writer, opt)?;
    }

    if opt.tint > 0 {
        writeln!(writer, "{:10.3}{:<50}{:<20}", opt.tint, "", "INTERVAL")?;
    }

    let mut ep = [0.0; 6];
    time2epoch(opt.ts, &mut ep);
    writeln!(writer, "  {:04.0}    {:02.0}    {:02.0}    {:02.0}    {:02.0}   {:010.7}     {:<12}{:<20}",
             ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], tsys, "TIME OF FIRST OBS")?;

    time2epoch(opt.te, &mut ep);
    writeln!(writer, "  {:04.0}    {:02.0}    {:02.0}    {:02.0}    {:02.0}   {:010.7}     {:<12}{:<20}",
             ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], tsys, "TIME OF LAST OBS")?;

    if opt.rnxver >= 301 {
        outrnx_phase_shift(writer, opt)?; // SYS / PHASE SHIFT
    }

    if opt.rnxver >= 302 {
        outrnx_glo_fcn(writer, opt, nav)?; // GLONASS SLOT / FRQ #
    }

    if opt.rnxver >= 302 {
        outrnx_glo_bias(writer, opt)?; // GLONASS COD/PHS/BIS
    }

    writeln!(writer, "{:<60.60}{:<20}", "", "END OF HEADER")?;
    Ok(())
}

// output data field in RINEX navigation data
fn outrnxnavf_n<W: Write>(writer: &mut W, value: f64, n: usize) -> Result<(), Error> {
    // 计算指数部分
    let e = if value.abs() < 1E-99 {
        0.0
    } else {
        value.abs().log10().floor() + 1.0
    };

    let coef = value.abs() / 10f64.powf(e - n as f64);

    write!(
        writer,
        " {}.{:0width$.0}{}{:03.0}",
        if value < 0.0 { "-" } else { " " },
        coef,
        "D",
        e,
        width = n
    )?;

    Ok(())
}

// output iono correction for a system
fn out_iono_sys<W: Write>(writer: &mut W, sys: &str, ion: &[f64], n: usize) -> Result<(), Error> {
    const LABEL1: [&str; 2] = ["ION ALPHA", "ION BETA"];
    const LABEL2: &str = "IONOSPHERIC CORR";

    if norm(ion) <= 0.0 {
        return Ok(());
    }

    for i in 0..(n + 3) / 4 {
        let mut str_buf = String::new();
        let label_suffix = if !sys.is_empty() && n >= 4 { 'A' as u8 + i as u8 } else { b' ' } as char;
        str_buf.push_str(sys);
        str_buf.push(label_suffix);

        write!(writer, "{:<4} ", if sys.is_empty() { "" } else { &str_buf })?;

        for j in 0..4.min(n - i * 4) {
            write!(writer, " ")?;
            outrnxnavf_n(writer, ion[i * 4 + j], 4)?;
        }

        write!(
            writer,
            "{:>width$}{}\n",
            "",
            if sys.is_empty() { LABEL1[i] } else { LABEL2 },
            width = if sys.is_empty() { 10 } else { 7 + 12 * (4 - 4.min(n - i * 4)) }
        )?;
    }

    Ok(())
}

// output iono corrections
fn out_iono<W: Write>(writer: &mut W, sys: usize, opt: &RnxOpt, nav: &Nav) -> Result<(), Error> {
    if opt.outiono == 0 {
        return Ok(());
    }

    if sys & opt.navsys & SYS_GPS != 0 {
        if opt.rnxver <= 211 {
            out_iono_sys(writer, "", &nav.ion_gps, 8)?;
        } else {
            out_iono_sys(writer, "GPS", &nav.ion_gps, 8)?;
        }
    }
    if sys & opt.navsys & SYS_GAL != 0 && opt.rnxver >= 212 {
        out_iono_sys(writer, "GAL", &nav.ion_gal, 3)?;
    }
    if sys & opt.navsys & SYS_QZS != 0 && opt.rnxver >= 302 {
        out_iono_sys(writer, "QZS", &nav.ion_qzs, 8)?;
    }
    if sys & opt.navsys & SYS_CMP != 0 && opt.rnxver >= 302 {
        out_iono_sys(writer, "BDS", &nav.ion_cmp, 8)?;
    }
    if sys & opt.navsys & SYS_IRN != 0 && opt.rnxver >= 303 {
        out_iono_sys(writer, "IRN", &nav.ion_irn, 8)?;
    }

    Ok(())
}

// output time system correction for a system
fn out_time_sys<W: Write>(writer: &mut W, sys: &str, utc: &[f64]) -> Result<(), Error> {
    const LABEL1: &str = "TIME SYSTEM CORR";
    const LABEL2: &str = "DELTA-UTC: A0,A1,T,W";

    if norm(utc) <= 0.0 {
        return Ok(());
    }

    if !sys.is_empty() {
        write!(writer, "{:<4} ", sys)?;
        outrnxnavf_n(writer, utc[0], 10)?;
        outrnxnavf_n(writer, utc[1], 9)?;
        write!(writer, "{:7.0}{:5.0}{:10}{}\n", utc[2], utc[3], "", LABEL1)?;
    } else {
        write!(writer, "   ")?;
        outrnxnavf_n(writer, utc[0], 12)?;
        outrnxnavf_n(writer, utc[1], 12)?;
        write!(writer, "{:9.0}{:9.0} {}\n", utc[2], utc[3], LABEL2)?;
    }

    Ok(())
}

// output time system corrections
fn out_time<W: Write>(writer: &mut W, sys: usize, opt: &RnxOpt, nav: &Nav) -> Result<(), Error> {
    let mut utc = [0.0; 8];

    if opt.outtime == 0 {
        return Ok(());
    }

    if sys & opt.navsys & SYS_GPS != 0 {
        if opt.rnxver <= 211 {
            out_time_sys(writer, "", &nav.utc_gps)?;
        } else {
            out_time_sys(writer, "GPUT", &nav.utc_gps)?;
        }
    }
    if sys & opt.navsys & SYS_GLO != 0 && opt.rnxver >= 212 {
        /* RINEX 2.12-3.02: tau_C, 3.03- : -tau_C */
        utc[0] = if opt.rnxver <= 302 { nav.utc_glo[0] } else { -nav.utc_glo[0] };
        out_time_sys(writer, "GLUT", &utc)?;
    }
    if sys & opt.navsys & SYS_GAL != 0 && opt.rnxver >= 212 {
        out_time_sys(writer, "GAUT", &nav.utc_gal)?;
    }
    if sys & opt.navsys & SYS_QZS != 0 && opt.rnxver >= 302 {
        out_time_sys(writer, "QZUT", &nav.utc_qzs)?;
    }
    if sys & opt.navsys & SYS_CMP != 0 && opt.rnxver >= 302 {
        out_time_sys(writer, "BDUT", &nav.utc_cmp)?;
    }
    if sys & opt.navsys & SYS_IRN != 0 && opt.rnxver >= 303 {
        out_time_sys(writer, "IRUT", &nav.utc_irn)?;
    }

    Ok(())
}

// output leap seconds
fn out_leaps<W: Write>(writer: &mut W, sys: usize, opt: &RnxOpt, nav: &Nav) -> Result<(), Error> {
    const LABEL: &str = "LEAP SECONDS";
    let leaps: &[f64];

    if opt.outleaps == 0 {
        return Ok(());
    }

    leaps = match sys {
        SYS_GAL => &nav.utc_gal[4..],
        SYS_QZS => &nav.utc_qzs[4..],
        SYS_CMP => &nav.utc_cmp[4..],
        SYS_IRN => &nav.utc_irn[4..],
        _ => &nav.utc_gps[4..],
    };

    if leaps[0] == 0.0 {
        return Ok(());
    }

    if opt.rnxver <= 300 {
        if sys == SYS_GPS {
            writeln!(writer, "{:6.0}{:54}{}", leaps[0], "", LABEL)?;
        }
    } else if norm(&leaps[1..]) <= 0.0 {
        writeln!(writer, "{:6.0}{:18}{:3}{:33}{}", leaps[0], "", if sys == SYS_CMP { "BDS" } else { "" }, "", LABEL)?;
    } else {
        writeln!(
            writer,
            "{:6.0}{:6.0}{:6.0}{:6.0}{:3}{:33}{}",
            leaps[0],
            leaps[3],
            leaps[1],
            leaps[2],
            if sys == SYS_CMP { "BDS" } else { "" },
            "",
            LABEL
        )?;
    }

    Ok(())
}

/// output navigation data header
fn outrnxnavh<W: Write>(writer: &mut W, opt: &RnxOpt, navs: Option<&Nav>) -> Result<(), Error> {
    let date = timestr_rnx();

    // RINEX version information and type
    if opt.rnxver <= 299 {
        writeln!(
            writer,
            "{:9.2}           {:<20}{:<20}{:<20}",
            opt.rnxver as f64 / 100.0,
            "N: GPS NAV DATA",
            "",
            "RINEX VERSION / TYPE"
        )?;
    } else {
        let sys = match opt.navsys {
            SYS_GPS => "G: GPS",
            SYS_GLO => "R: GLONASS",
            SYS_GAL => "E: Galileo",
            SYS_QZS => "J: QZSS",    // v.3.02
            SYS_CMP => "C: BeiDou",  // v.3.02
            SYS_IRN => "I: IRNSS",   // v.3.03
            _ if opt.sep_nav != 0 => "G: GPS",
            _ => "M: Mixed",
        };
        writeln!(
            writer,
            "{:9.2}           {:<20}{:<20}{:<20}",
            opt.rnxver as f64 / 100.0,
            "N: GNSS NAV DATA",
            sys,
            "RINEX VERSION / TYPE"
        )?;
    }

    writeln!(
        writer,
        "{:<20.20}{:<20.20}{:<20.20}{:<20}",
        opt.prog, opt.runby, date, "PGM / RUN BY / DATE"
    )?;

    for comment in opt.comment.iter() {
        if comment.is_empty() {
            continue;
        }
        writeln!(writer, "{:<60.60}{:<20}", comment, "COMMENT")?;
    }

    if let Some(nav) = navs {
        out_iono(writer, if opt.sep_nav != 0 { SYS_GPS } else { SYS_ALL }, opt, nav)?;
        out_time(writer, if opt.sep_nav != 0 { SYS_GPS } else { SYS_ALL }, opt, nav)?;
        out_leaps(writer, SYS_GPS, opt, nav)?;
    }

    writeln!(writer, "{:60}{:<20}", "", "END OF HEADER")?;

    Ok(())
}

// sort obs-types
fn sort_obstype(
    codes: &mut Vec<u8>,
    types: &mut [i32; 33],
    syslog: usize,
) {
    if codes.len() == 0 { return; }
    for i in 0..codes.len() - 1 {
        for j in i + 1..codes.len() {
            let idx1 = code2idx(1 << syslog, codes[i]);
            let idx2 = code2idx(1 << syslog, codes[j]);
            let pri1 = getcodepri(1 << syslog, codes[i], Some(""));
            let pri2 = getcodepri(1 << syslog, codes[j], Some(""));
            if idx1 < idx2 || (idx1 == idx2 && pri1 >= pri2) { continue; }
            let temp = codes[i];
            codes[i] = codes[j];
            codes[j] = temp;
            let temp = types[i];
            types[i] = types[j];
            types[j] = temp;
        }
    }
}

// convert RINEX obs-type ver.3 -> ver.2
fn convcode(rnxver: i32, sys: usize, type_vec: &mut Vec<char>) {
    let type_str: String = type_vec.iter().collect();

    if rnxver >= 212 && (sys == SYS_GPS || sys == SYS_QZS) && &type_str[1..] == "1C" {
        type_vec.truncate(1);
        type_vec.push('A');
    } else if rnxver >= 212 && (sys == SYS_GPS || sys == SYS_QZS) && (&type_str[1..] == "1S" || &type_str[1..] == "1L" || &type_str[1..] == "1X") {
        type_vec.truncate(1);
        type_vec.push('B');
    } else if rnxver >= 212 && (sys == SYS_GPS || sys == SYS_QZS) && (&type_str[1..] == "2S" || &type_str[1..] == "2L" || &type_str[1..] == "2X") {
        type_vec.truncate(1);
        type_vec.push('C');
    } else if rnxver >= 212 && sys == SYS_GLO && &type_str[1..] == "1C" {
        type_vec.truncate(1);
        type_vec.push('A');
    } else if rnxver >= 212 && sys == SYS_GLO && &type_str[1..] == "2C" {
        type_vec.truncate(1);
        type_vec.push('D');
    } else if sys == SYS_CMP && (&type_str[1..] == "2I" || &type_str[1..] == "2Q" || &type_str[1..] == "2X") {
        type_vec.truncate(1);
        type_vec.push('2');
    } else if type_str == "C1P" || &type_str == "C1W" || &type_str == "C1Y" || &type_str == "C1N" {
        type_vec.clear();
        type_vec.push('P');
        type_vec.push('1');
    } else if type_str == "C2P" || &type_str == "C2W" || &type_str == "C2Y" || &type_str == "C2N" || &type_str == "C2D" {
        type_vec.clear();
        type_vec.push('P');
        type_vec.push('2');
    } else {
        type_vec.truncate(2);
    }
}

// set obs-types in RINEX options
fn setopt_obstype(
    codes: &Vec<u8>,
    types: &[i32; 33],
    syslog: usize,
    opt: &mut RnxOpt,
) {
    let type_str = "CLDS";
    let mut type_buf: Vec<char>;

    opt.nobs[syslog] = 0;

    if (1 << syslog & opt.navsys) == 0 { return; }

    for i in 0..codes.len() {
        let id = code2obs(codes[i]);
        let idx = code2idx(1 << syslog, codes[i]);
        if id.is_empty() || idx < 0 { continue; }
        if (opt.freqtype & (1 << idx)) == 0 || opt.mask[syslog][codes[i] as usize - 1] == '0' { continue; }
        if opt.rnxver >= 300 {
            let ver = VERCODE[syslog].as_bytes()[(codes[i] - 1) as usize] as char;
            if ver < '0' || ver > ('0' as u8 + (opt.rnxver - 300) as u8) as char {
                continue;
            }
        }

        for j in 0..4 {
            if (opt.obstype & (1 << j)) == 0 { continue; }
            if (types[i] & (1 << j)) == 0 { continue; }
            type_buf = Vec::new();
            if let Some(ch) = type_str.chars().nth(j) {
                type_buf.push(ch);
            }
            for ch in id.chars() {
                type_buf.push(ch);
            }
            if type_buf.get(0) == Some(&'C') && type_buf.get(2) == Some(&'N') {
                return;
            }

            if opt.rnxver <= 299 {
                /* ver.3 -> ver.2 */
                convcode(opt.rnxver, 1 << syslog, &mut type_buf);
                /* check duplicated obs-type */
                let mut flag = true;
                for k in 0..opt.nobs[0] {
                    let tobs = opt.tobs[0][k as usize].iter().collect::<String>();
                    if tobs.trim_end_matches('\0') == type_buf.iter().collect::<String>() {
                        flag = false;
                        break;
                    }
                }
                if flag && opt.nobs[0] < MAXOBSTYPE as i32 {
                    for (idx, &ch) in type_buf.iter().enumerate() {
                        if idx < 4 {
                            opt.tobs[0][opt.nobs[0] as usize][idx] = ch;
                        }
                    }
                    opt.nobs[0] += 1;
                }
            } else if opt.nobs[syslog] < MAXOBSTYPE as i32 {
                for (idx, &ch) in type_buf.iter().enumerate() {
                    if idx < 4 {
                        opt.tobs[syslog][opt.nobs[syslog] as usize][idx] = ch;
                    }
                }
                opt.nobs[syslog] += 1;
            }
        }
    }
}

// search obsservattion data index
fn obsindex(
    rnxver: i32,
    sys: usize,
    code: &[u8],
    tobs: &[char],
    mask: &[char]) -> i32 {
    for i in 0..(NFREQ + NEXOBS) {
        // signal mask
        if code[i] == 0 || mask[(code[i] - 1) as usize] == '0' {
            continue;
        }

        if rnxver <= 299 {
            // ver.2
            if (tobs[0] == 'C' && tobs[1] == '1') && (sys == SYS_GPS || sys == SYS_GLO || sys == SYS_QZS || sys == SYS_CMP) {
                if code[i] == CODE_L1C {
                    return i as i32;
                }
            } else if tobs[0] == 'P' && tobs[1] == '1' {
                if code[i] == CODE_L1P || code[i] == CODE_L1W || code[i] == CODE_L1Y || code[i] == CODE_L1N {
                    return i as i32;
                }
            } else if (tobs[0] == 'C' && tobs[1] == '2') && (sys == SYS_GPS || sys == SYS_QZS) {
                if code[i] == CODE_L2S || code[i] == CODE_L2L || code[i] == CODE_L2X {
                    return i as i32;
                }
            } else if (tobs[0] == 'C' && tobs[1] == '2') && (sys == SYS_GLO) {
                if code[i] == CODE_L2C {
                    return i as i32;
                }
            } else if tobs[0] == 'P' && tobs[1] == '2' {
                if code[i] == CODE_L2P || code[i] == CODE_L2W || code[i] == CODE_L2Y || code[i] == CODE_L2N || code[i] == CODE_L2D {
                    return i as i32;
                }
            } else if rnxver >= 212 && tobs[1] == 'A' {
                // L1C/A
                if code[i] == CODE_L1C {
                    return i as i32;
                }
            } else if rnxver >= 212 && tobs[1] == 'B' {
                // L1C
                if code[i] == CODE_L1S || code[i] == CODE_L1L || code[i] == CODE_L1X {
                    return i as i32;
                }
            } else if rnxver >= 212 && tobs[1] == 'C' {
                // L2C
                if code[i] == CODE_L2S || code[i] == CODE_L2L || code[i] == CODE_L2X {
                    return i as i32;
                }
            } else if rnxver >= 212 && tobs[1] == 'D' && sys == SYS_GLO {
                // GLO L2C/A
                if code[i] == CODE_L2C {
                    return i as i32;
                }
            } else if tobs[1] == '2' && sys == SYS_CMP {
                // BDS B1
                if code[i] == CODE_L2I || code[i] == CODE_L2Q || code[i] == CODE_L2X {
                    return i as i32;
                }
            } else {
                let id = code2obs(code[i]);
                if id.chars().nth(0) == Some(tobs[1]) {
                    return i as i32;
                }
            }
        } else {
            // ver.3
            let id = code2obs(code[i]);
            let tob: Vec<char> = tobs[1..3].iter().copied().collect();
            if id == &tob.into_iter().collect::<String>() {
                return i as i32;
            }
        }
    }
    -1
}

// output observation data field
fn outrnxobsf<W: Write>(
    writer: &mut W,
    obs: f64,
    lli: i32,
) -> Result<(), Error> {
    if obs == 0.0 || obs <= -1E9 || obs >= 1E9 { write!(writer, "              ")? } else { write!(writer, "{:14.3}", obs)? }

    if lli < 0 || (lli & (LLI_SLIP | LLI_HALFC | LLI_BOCTRK)) == 0 { write!(writer, "  ")? } else { write!(writer, "{:1.1} ", lli & (LLI_SLIP | LLI_HALFC | LLI_BOCTRK))? }

    Ok(())
}

// output navigation data field
fn outrnxnavf<W: Write>(writer: &mut W, value: f64) -> Result<(), Error> {
    outrnxnavf_n(writer, value, 12)
}

// update station list
fn update_stas(stas_list: &mut LinkedList<Stas>, staid: i32, time: GTime) {
    if let Some(front) = stas_list.front_mut() {
        if front.staid == staid {
            front.te = time;
            return;
        }
    }

    let new_sta = Stas {
        staid,
        ts: time,
        te: time,
        sta: Sta::default(),
    };
    stas_list.push_front(new_sta);
}

// update station info in station list
fn update_stainf(stas_list: &mut LinkedList<Stas>, sta: &Sta, staid: i32) {
    if let Some(front) = stas_list.front_mut() {
        if front.staid == staid {
            front.sta = sta.clone()
        }
    }
}

// set station ID list to RINEX options comments
fn setopt_sta_list(stas_list: &LinkedList<Stas>, opt: &mut RnxOpt) {
    let n = stas_list.len();

    if n <= 1 {
        return;
    }

    let mut i = 0;
    while i < MAXCOMMENT && !opt.comment[i].is_empty() {
        i += 1;
    }

    if i < MAXCOMMENT {
        opt.comment[i] = format!("{:5}  {:22}  {:22}", "STAID", "TIME OF FIRST OBS", "TIME OF LAST OBS");
        i += 1;
    }

    for (index, stas) in stas_list.iter().enumerate().rev() {
        if i + index >= MAXCOMMENT {
            continue;
        }
        let s1 = time2str(stas.ts);
        let s2 = time2str(stas.te);
        opt.comment[i + index] = format!(" {:04}  {}  {}", stas.staid, s1, s2);
    }
}

// set station info in RINEX options
fn setopt_sta<'a>(stas_list: &'a LinkedList<Stas>, mut sta: &'a Sta, opt: &mut RnxOpt) {
    let mut s = None;
    let mut iter = stas_list.iter().peekable();
    while let Some(stas) = iter.next() {
        if iter.peek().is_none() || (opt.ts.time != 0 && timediff(stas.te, opt.ts) < 0.0) {
            s = Some(&stas.sta);
            break;
        }
    }

    if s.is_some() {
        setopt_sta_list(stas_list, opt);
        sta = s.unwrap();
    }

    /* marker name and number */
    if opt.marker.is_empty() && opt.markerno.is_empty() {
        opt.marker = sta.name.clone();
        opt.markerno = sta.marker.clone();
    }

    /* receiver and antenna info */
    if opt.rec[0].is_empty() && opt.rec[1].is_empty() && opt.rec[2].is_empty() {
        opt.rec[0] = sta.recsno.clone();
        opt.rec[1] = sta.rectype.clone();
        opt.rec[2] = sta.recver.clone();
    }
    if opt.ant[0].is_empty() && opt.ant[1].is_empty() && opt.ant[2].is_empty() {
        opt.ant[0] = sta.antsno.clone();
        opt.ant[1] = sta.antdes.clone();
        if sta.antsetup == 0 {
            opt.ant[2] = format!("{}", sta.antsetup).to_string();
        } else { opt.ant[2] = "\0".to_string(); }
    }

    /* antenna approx position */
    if opt.autopos == 0 && norm(&sta.pos) > 0.0 {
        opt.apppos = sta.pos;
    }

    /* antenna delta */
    if norm(&opt.antdel) <= 0.0 && norm(&sta.del) > 0.0 {
        if sta.deltype == 0 {
            opt.antdel[0] = sta.del[2];
            opt.antdel[1] = sta.del[0];
            opt.antdel[2] = sta.del[1];
        } else if norm(&sta.pos) > 0.0 {
            let pos = &mut [];
            ecef2pos(&sta.pos, pos);
            let enu = ecef2enu(pos, &sta.del);
            opt.antdel[0] = enu[2];
            opt.antdel[1] = enu[0];
            opt.antdel[2] = enu[1];
        }
    } else {
        opt.antdel[0] = sta.hgt;
        opt.antdel[1] = 0.0;
        opt.antdel[2] = 0.0;
    }
}

// set phase shift in RINEX options (RINEX 3.04 A23)
fn setopt_phshift(opt: &mut RnxOpt) {
    for i in 0..NUMSYS {
        for j in 0..opt.nobs[i] as usize {
            if opt.tobs[i][j][0] != 'L' {
                continue;
            }
            let code = obs2code(&opt.tobs[i][j][1..3].iter().collect::<String>());

            opt.shift[i][j] = match 1 << i {
                SYS_GPS => {
                    if matches!(code, CODE_L1S | CODE_L1L | CODE_L1X | CODE_L1P | CODE_L1W | CODE_L1N) {
                        0.25 // +1/4 cyc
                    } else if matches!(code, CODE_L2C | CODE_L2S | CODE_L2L | CODE_L2X | CODE_L5Q) {
                        -0.25 // -1/4 cyc
                    } else {
                        continue;
                    }
                }
                SYS_GLO => {
                    if matches!(code, CODE_L1P | CODE_L2P | CODE_L3Q) {
                        0.25 // +1/4 cyc
                    } else {
                        continue;
                    }
                }
                SYS_GAL => {
                    if code == CODE_L1C {
                        0.5 // +1/2 cyc
                    } else if matches!(code, CODE_L5Q | CODE_L7Q | CODE_L8Q) {
                        -0.25 // -1/4 cyc
                    } else if code == CODE_L6C {
                        -0.5 // -1/2 cyc
                    } else {
                        continue;
                    }
                }
                SYS_QZS => {
                    if matches!(code, CODE_L1S | CODE_L1L | CODE_L1X) {
                        0.25 // +1/4 cyc
                    } else if matches!(code, CODE_L5Q | CODE_L5P) {
                        -0.25 // -1/4 cyc
                    } else {
                        continue;
                    }
                }
                SYS_CMP => {
                    if matches!(code, CODE_L2P | CODE_L7Q | CODE_L6Q) {
                        -0.25 // -1/4 cyc
                    } else if matches!(code, CODE_L1P | CODE_L5P | CODE_L7P) {
                        0.25 // +1/4 cyc
                    } else {
                        continue;
                    }
                }
                _ => continue,
            };
        }
    }
}

// scan rtcm file
fn scan_rtcm(
    decoder: &Decoder,
    rtcmfile: &str,
    navs: &mut Option<&mut Nav>,
    opt: &mut RnxOpt,
    tr: GTime,
) -> Result<usize, DecodeError> {
    let path = Path::new(rtcmfile);
    let file = File::open(path).map_err(|_| IoError)?;
    let mut reader = BufReader::new(file);

    let mut epoch_time = 0;
    let mut obss = Obss::new();
    let mut deocde_type = 0;
    let mut scan_flag = false;
    let mut staid;
    let mut stas_list: LinkedList<Stas> = LinkedList::new();
    let mut sta = Sta::default();
    let mut te = 0;
    let mut ts = 0;
    let mut time = GTime::default();
    let mut types = [[0; 33]; NUMSYS];
    let mut codes: [Vec<u8>; NUMSYS] = std::array::from_fn(|_| Vec::new());
    let mut n = 0;
    let mut ati = AutoTi::new();
    let mut trt = tr;

    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap()
        .tick_strings(&["-", "\\", "|", "/"]));

    loop {
        let single_rtcm;
        let pos = reader.seek(SeekFrom::Current(0)).map_err(|_| IoError)?;

        match read_single_rtcm(&mut reader) {
            Ok(value) => {
                if value.len() <= 0 {
                    continue;
                }
                single_rtcm = value;
            }
            Err(e) => match e {
                IoError => break,
                EndofFile => break,
                _ => continue
            }
        };

        match decoder.decode_rtcm(trt, &single_rtcm) {
            None => continue,
            Some(rtcm) => {
                if rtcm.obss.n > 0 {
                    deocde_type = 1;
                    if ts == 0 { ts = rtcm.time.time }
                    if te == 0 { te = rtcm.time.time }
                    if epoch_time == 0 { epoch_time = rtcm.time.time }
                    // epoch time change, previous calendar element obs not properly terminated outputs
                    if epoch_time != rtcm.time.time {
                        reader.seek(SeekFrom::Start(pos)).map_err(|_| IoError)?;
                        scan_flag = true;
                        ts = min(rtcm.time.time, ts.clone());
                        te = max(rtcm.time.time, te.clone());
                        time = rtcm.time;
                        epoch_time = 0;
                    }
                    if !scan_flag {
                        for i in 0..rtcm.obss.n {
                            let _ = obss.add_obs_data(&rtcm.obss.obs[i]);
                        }
                    }
                    // Determine if the current calendar element obs is finished
                    if rtcm.obsflag {
                        scan_flag = true;
                        ts = min(rtcm.time.time, ts.clone());
                        te = max(rtcm.time.time, te.clone());
                        time = rtcm.time;
                        epoch_time = 0;
                    }
                    trt=rtcm.time;
                } else if rtcm.ephsat != 0 {
                    deocde_type = 2;
                    scan_flag = true;
                    ts = min(rtcm.time.time, ts.clone());
                    te = max(rtcm.time.time, te.clone());
                    time = rtcm.time;
                    if let Some(nav) = navs {
                        if satsys(rtcm.ephsat) == SYS_GLO {
                            nav.add_geph(rtcm.geph);
                        } else { nav.add_eph(rtcm.eph); }
                    }
                    trt=rtcm.time;
                } else if !rtcm.sta.name.is_empty() {
                    deocde_type = 5;
                    scan_flag = true;
                    sta = rtcm.sta;
                }
                staid = rtcm.staid;
            }
        }

        if scan_flag {
            if opt.ts.time != 0 && time.time != 0 && timediff(time, opt.ts) < 0.0 { continue; }
            if opt.te.time != 0 && time.time != 0 && timediff(time, opt.te) > 0.0 { break; }

            if deocde_type == 1 {
                for i in 0..obss.obs.len() {
                    let obs = obss.obs[i].clone();
                    let sys = satsys(obs.sat);
                    if (sys & opt.navsys) == 0 { continue; }
                    let syslog = sys.ilog2() as usize;
                    for j in 0..NFREOBS {
                        if obs.code[j] == 0 { continue; }
                        let mut flag = true;
                        for (idx, data) in codes[syslog].iter().enumerate() {
                            if *data == obs.code[j] {
                                flag = false;
                                if obs.p[j] != 0.0 { types[syslog][idx] |= 1 }
                                if obs.l[j] != 0.0 { types[syslog][idx] |= 2 }
                                if obs.d[j] != 0.0 { types[syslog][idx] |= 4 }
                                if obs.snr[j] != 0 { types[syslog][idx] |= 8 }
                            }
                        }
                        if flag && codes[syslog].len() < 32 {
                            codes[syslog].push(obs.code[j])
                        }
                    }
                }
                /* update station list */
                update_stas(&mut stas_list, staid, time);

                cal_sample(&mut ati, obss.obs[0].time);
                trt = obss.obs[0].time;
                n += 1;
                let message = format!("Scanning rtcm file: epoch:{}", n);
                pb.set_message(message.to_string());
                pb.tick();

                obss = Obss::new();
            } else if deocde_type == 5 {
                update_stainf(&mut stas_list, &sta, staid);
            }
            scan_flag = false;
        }
    }

    /* sort and set obs-types in RINEX options */
    for i in 0..NUMSYS {
        sort_obstype(&mut codes[i], &mut types[i], i);
        setopt_obstype(&codes[i], &types[i], i, opt);
    }

    /* set station info in RINEX options */
    setopt_sta(&stas_list, &sta, opt);

    /* set phase shifts in RINEX options */
    if opt.phshift {
        setopt_phshift(opt);
    }

    /* set GLONASS FCN*/
    if let Some(nav) = navs {
        for geph in nav.clone().geph {
            let prn = getprn(geph.sat);
            nav.glo_fcn[prn - 1] = geph.frq + 8;
        }
    }

    if opt.ts.time == 0 && ts != 0 { opt.ts.time = ts }
    if opt.te.time == 0 && te != 0 { opt.te.time = te }
    if opt.tint == 0 { opt.tint = ati.ti }

    pb.finish_with_message(String::from("Finish scan"));

    Ok(n as usize)
}

// output RINEX observation data body
fn outrnxobsb<W: Write>(
    writer: &mut W,
    opt: &RnxOpt,
    obs: &Vec<Obs>,
) -> Result<(), Error> {
    let mut sats = Vec::new();
    let mut s = Vec::new();
    let mut ind = [0; MAXOBS];
    let mut ep = [0.0; 6];
    let mut mask: [char; 64];
    let mut m;
    time2epoch(obs[0].time, &mut ep);

    let mut ns = 0;
    for i in 0..obs.len() {
        if ns >= MAXOBS { break; }
        let sys = satsys(obs[i].sat);
        let id = satno2id(obs[i].sat);
        if (sys & opt.navsys) == 0 || opt.exsats[obs[i].sat - 1] != 0 { continue; }
        if id.is_empty() { continue; }
        sats.push(id);
        s.push(sys.ilog2() as usize);
        if opt.nobs[if opt.rnxver <= 299 { 0 } else { s[ns] }] == 0 { continue; }
        ind[ns] = i;
        ns += 1;
    }

    if ns <= 0 { return Ok(()); }

    if opt.rnxver <= 299 {
        write!(writer, " {:02} {:02.0} {:02.0} {:02.0} {:02.0} {:010.7}  {}{:3}",
               ep[0] as i32 % 100, ep[1], ep[2], ep[3], ep[4], ep[5], 0, ns)?;
        for i in 0..ns {
            if i > 0 && i % 12 == 0 { write!(writer, "\n{:>32}", "")?; }
            write!(writer, "{:3}", sats[i])?;
        }
    } else {
        writeln!(writer, "> {:04.0} {:02.0} {:02.0} {:02.0} {:02.0} {:010.7}  {}{:3}{:>21}",
                 ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], 0, ns, "")?;
    }

    for i in 0..ns {
        let sys = satsys(obs[ind[i]].sat);

        if opt.rnxver <= 299 {
            m = 0;
            mask = opt.mask[s[i]];
        } else {
            write!(writer, "{:3}", sats[i])?;
            m = s[i];
            mask = opt.mask[s[i]];
        }

        let mut k;
        for j in 0..opt.nobs[m] as usize {
            if opt.rnxver <= 299 && j % 5 == 0 { writeln!(writer, "")?; }

            k = obsindex(opt.rnxver, sys, &obs[ind[i]].code, &opt.tobs[m][j], &mask);
            if k < 0 {
                outrnxobsf(writer, 0.0, -1)?;
                continue;
            }
            /* phase shift (cyc) */
            let dl = if obs[ind[i]].l[k as usize] != 0.0 { opt.shift[m][j] } else { 0.0 };
            /* output field */
            match opt.tobs[m][j][0] {
                'C' => outrnxobsf(writer, obs[ind[i]].p[k as usize], -1)?,
                'P' => outrnxobsf(writer, obs[ind[i]].p[k as usize], -1)?,
                'L' => outrnxobsf(writer, obs[ind[i]].l[k as usize] + dl, obs[ind[i]].lli[k as usize])?,
                'D' => outrnxobsf(writer, obs[ind[i]].d[k as usize] as f64, -1)?,
                'S' => outrnxobsf(writer, obs[ind[i]].snr[k as usize] as f64 / 1000.0, -1)?,
                _ => ()
            }
        }
        if opt.rnxver >= 300 && writeln!(writer, "").is_err() { return Ok(()); }
    }

    if opt.rnxver >= 300 { return Ok(()); }
    writeln!(writer,"")?;
    Ok(())
}

// output glonass navigation data body
fn outrnxgnavb<W: Write>(writer: &mut W, opt: &RnxOpt, geph: &Geph) -> Result<(), Error> {
    let prn = getprn(geph.sat);
    let sep;

    // Check if satellite system is GLO and is part of selected navigation systems
    if satsys(geph.sat) & opt.navsys != SYS_GLO {
        return Ok(());
    }

    // Convert gpst to utc and calculate tof
    let mut tof = time2gpst(gpst2utc(geph.tof), None);
    if opt.rnxver <= 299 {
        tof = tof % 86400.0; // Ver.2: tod in UTC
    }

    let toe = gpst2utc(geph.toe); // Convert GPST to UTC
    let mut ep = [0.0; 6];
    time2epoch(toe, &mut ep); // Get time components as [year, month, day, hour, minute, second]

    if opt.rnxver <= 299 {
        // RINEX version 2
        write!(writer, "{:2} {:02} {:02.0} {:02.0} {:02.0} {:02.0} {:04.1}",
               prn,
               (ep[0] as i32 % 100), ep[1], ep[2], ep[3], ep[4], ep[5])?;
        sep = "   ";
    } else {
        // RINEX version 3
        let code = sat2code(geph.sat);
        if code.is_empty() {
            return Ok(());
        }
        write!(writer, "{:<-3} {:04.0} {:02.0} {:02.0} {:02.0} {:02.0} {:02.0}",
               code, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5])?;
        sep = "    ";
    }

    // Output geph elements
    outrnxnavf(writer, -geph.taun)?;
    outrnxnavf(writer, geph.gamn)?;
    outrnxnavf(writer, tof)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, geph.pos[0] / 1E3)?;
    outrnxnavf(writer, geph.vel[0] / 1E3)?;
    outrnxnavf(writer, geph.acc[0] / 1E3)?;
    outrnxnavf(writer, geph.svh as f64)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, geph.pos[1] / 1E3)?;
    outrnxnavf(writer, geph.vel[1] / 1E3)?;
    outrnxnavf(writer, geph.acc[1] / 1E3)?;
    outrnxnavf(writer, geph.frq as f64)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, geph.pos[2] / 1E3)?;
    outrnxnavf(writer, geph.vel[2] / 1E3)?;
    outrnxnavf(writer, geph.acc[2] / 1E3)?;

    // Use dtaun or age based on a compile-time condition
    #[cfg(feature = "use_dtaun")]
    outrnxnavf(writer, geph.dtaun)?;
    #[cfg(not(feature = "use_dtaun"))]
    outrnxnavf(writer, geph.age as f64)?;

    // Final newline
    writeln!(writer)?;

    Ok(())
}

// output navigation data body
fn outrnxnavb<W: Write>(writer: &mut W, opt: &RnxOpt, eph: &Eph) -> Result<(), Error> {
    let sys = satsys(eph.sat);
    let prn = getprn(eph.sat);
    let sep;
    let mut ep = [0.0; 6]; // Time components [year, month, day, hour, minute, second]
    let ttr;
    let mut week = 0;

    // Check satellite system and navigation system
    if sys == 0 || (sys & opt.navsys) == 0 {
        return Ok(());
    }

    // Time conversion (toc to ep)
    if sys != SYS_CMP {
        time2epoch(eph.toc, &mut ep);
    } else {
        time2epoch(gpst2bdt(eph.toc), &mut ep); // Convert GPST to BDT
    }

    // Output satellite info based on RINEX version
    if (opt.rnxver >= 300 && sys == SYS_GPS) ||
        (opt.rnxver >= 212 && sys == SYS_GAL) ||
        (opt.rnxver >= 302 && sys == SYS_QZS) ||
        (opt.rnxver >= 302 && sys == SYS_CMP) ||
        (opt.rnxver >= 303 && sys == SYS_IRN) {
        let code = sat2code(eph.sat);
        if code.is_empty() {
            return Ok(());
        }
        write!(writer, "{:<3} {:04.0} {:02.0} {:02.0} {:02.0} {:02.0} {:02.0}",
               code, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5])?;
        sep = "    ";
    } else if opt.rnxver <= 299 && sys == SYS_GPS {
        write!(writer, "{:2} {:02} {:02.0} {:02.0} {:02.0} {:02.0} {:04.1}",
               prn, (ep[0] as i32 % 100), ep[1], ep[2], ep[3], ep[4], ep[5])?;
        sep = "   ";
    } else {
        return Ok(());
    }

    // Output other fields
    outrnxnavf(writer, eph.f0)?;
    outrnxnavf(writer, eph.f1)?;
    outrnxnavf(writer, eph.f2)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, eph.iode as f64)?; // GPS/QZS: IODE, GAL: IODnav, BDS: AODE
    outrnxnavf(writer, eph.crs)?;
    outrnxnavf(writer, eph.deln)?;
    outrnxnavf(writer, eph.m0 as f64)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, eph.cuc)?;
    outrnxnavf(writer, eph.e)?;
    outrnxnavf(writer, eph.cus)?;
    outrnxnavf(writer, eph.a.sqrt())?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, eph.toes)?;
    outrnxnavf(writer, eph.cic)?;
    outrnxnavf(writer, eph.omg0)?;
    outrnxnavf(writer, eph.cis)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, eph.i0)?;
    outrnxnavf(writer, eph.crc)?;
    outrnxnavf(writer, eph.omg)?;
    outrnxnavf(writer, eph.omgd)?;
    write!(writer, "\n{}", sep)?;

    outrnxnavf(writer, eph.idot)?;
    outrnxnavf(writer, eph.code as f64)?;
    outrnxnavf(writer, eph.week as f64)?; // GPS/QZS: GPS week, GAL: GAL week, BDS: BDT week
    if sys == SYS_GPS || sys == SYS_QZS {
        outrnxnavf(writer, eph.flag as f64)?;
    } else {
        outrnxnavf(writer, 0.0)?; // spare
    }
    write!(writer, "\n{}", sep)?;

    // Handle SISA and URA values based on system
    if sys == SYS_GAL {
        outrnxnavf(writer, sisa_value(eph.sva))?;
    } else {
        outrnxnavf(writer, uravalue(eph.sva as usize))?;
    }
    outrnxnavf(writer, eph.svh as f64)?;
    outrnxnavf(writer, eph.tgd[0])?; // GPS/QZS: TGD, GAL: BGD E5a/E1, BDS: TGD1 B1/B3
    if sys == SYS_GAL || sys == SYS_CMP {
        outrnxnavf(writer, eph.tgd[1])?; // GAL: BGD E5b/E1, BDS: TGD2 B2/B3
    } else if sys == SYS_GPS || sys == SYS_QZS {
        outrnxnavf(writer, eph.iodc as f64)?; // GPS/QZS: IODC
    } else {
        outrnxnavf(writer, 0.0)?; // spare
    }
    write!(writer, "\n{}", sep)?;

    // Time and week handling
    if sys != SYS_CMP {
        ttr = time2gpst(eph.ttr, Some(&mut week));
    } else {
        ttr = time2bdt(gpst2bdt(eph.ttr), &mut week); // GPST -> BDT
    }
    outrnxnavf(writer, ttr + (week - eph.week) as f64 * 604800.0)?;

    // Additional fields based on system
    if sys == SYS_GPS {
        outrnxnavf(writer, eph.fit)?;
    } else if sys == SYS_QZS {
        outrnxnavf(writer, if eph.fit > 2.0 { 1.0 } else { 0.0 })?;
    } else if sys == SYS_CMP {
        outrnxnavf(writer, eph.iodc as f64)?; // AODC
    } else {
        outrnxnavf(writer, 0.0)?; // spare
    }
    writeln!(writer)?;

    Ok(())
}

/// decode rtcm3 file to rinex
pub fn rtcm2rnx(
    rtcmfile: &str,
    obsfile: Option<String>,
    navfile: Option<String>,
    qc_opt: QcOpt,
) -> Result<bool, DecodeError> {
    if obsfile.is_none() && navfile.is_none() || rtcmfile.is_empty() { return Ok(false); }

    /* open rtcm file */
    let path = Path::new(rtcmfile);
    let file = File::open(path).map_err(|_| IoError)?;
    let mut reader = BufReader::new(file);

    // default option
    let mut trt = qc_opt.tr;
    let mut opt = RnxOpt::new();
    opt.obstype = 15; // All types: pseudorange/carrier/Doppler/signal-to-noise ratio
    opt.trtcm = qc_opt.tr;
    opt.rnxver = qc_opt.rnxver;
    opt.freqtype = 31; // five frequencies
    if opt.rnxver <= 210 { opt.freqtype &= 0x3 }
    opt.navsys = qc_opt.navsys;
    let sys_gr = SYS_GPS | SYS_GLO;
    if opt.rnxver <= 210 {
        opt.navsys &= sys_gr
    } else if opt.rnxver <= 211 {
        opt.navsys &= sys_gr | SYS_GAL
    } else if opt.rnxver <= 212 {
        opt.navsys &= sys_gr | SYS_GAL | SYS_CMP
    } else if opt.rnxver <= 300 {
        opt.navsys &= sys_gr | SYS_GAL
    } else if opt.rnxver <= 301 {
        opt.navsys &= sys_gr | SYS_GAL | SYS_CMP
    } else if opt.rnxver <= 302 {
        opt.navsys &= sys_gr | SYS_GAL | SYS_CMP | SYS_QZS
    }

    for i in 0..NUMSYS { for j in 0..64 { opt.mask[i][j] = '1' } }
    /* set rtcm3 comment for rnx */
    opt.comment[0] = "format: RTCM3".to_string();
    opt.comment[1] = format!("inp: {:.56}", rtcmfile);
    opt.prog = "GDA".to_string();
    opt.runby = "Wuhan Shouming Technology".to_string();
    //opt.phshift=true;
    /* get decoder*/
    let decoder = Decoder::new();
    /* scan rtcm and get nav */
    let mut navs = Nav::new();
    let n = scan_rtcm(&decoder, rtcmfile, &mut Some(&mut navs), &mut opt, trt)?;
    /* processing count */
    let mut count = 0;
    /* restore slips */
    let mut slips = [[0; NFREOBS]; MAXSAT];
    let mut first_ep = true;
    let mut obs_flag = false;

    /*open output file */
    let mut obs_writer: Option<BufWriter<File>> = None;
    if let Some(filename) = obsfile.clone() {
        let path = Path::new(&filename);
        if let Some(parent) = path.parent() {
            create_dir_all(parent).expect("Fail to create path");
        }
        let file = File::create(&filename).expect("Fail to create output obs file!");
        obs_writer = Some(BufWriter::new(file));
        if let Some(writer) = obs_writer.as_mut() {
            outrnxobsh(writer, &opt, Some(&navs)).expect("Fail to output header")
        }
    }
    if let Some(filename) = navfile.clone() {
        let path = Path::new(&filename);
        if let Some(parent) = path.parent() {
            create_dir_all(parent).expect("Fail to create path");
        }
        let file = File::create(&filename).expect("Fail to create output nav file!");
        let mut nav_writer = BufWriter::new(file);
        if navs.n + navs.ng <= 0 {
            remove_file(filename).expect("Failed to delete file")
        } else {
            outrnxnavh(&mut nav_writer, &opt, Some(&navs)).expect("Fail to output header");
            for eph in navs.eph {
                let _ = outrnxnavb(&mut nav_writer, &opt, &eph);
            }
            for geph in navs.geph {
                let _ = outrnxgnavb(&mut nav_writer, &opt, &geph);
            }
        }
    }

    let pb = ProgressBar::new(n as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}")
        .unwrap().progress_chars("#>-"));

    loop {
        if let Some(writer) = obs_writer.as_mut() {
            let mut obss = Obss::new();
            match decode_epoch(&decoder, &mut reader, &mut obss, None, trt) {
                Err(e) => match e {
                    IoError => break,
                    EndofFile => break,
                    _ => continue
                }
                Ok(true) => {
                    /* set slip */
                    save_slips(&mut slips, &obss.obs);
                    rest_slips(&mut slips, &mut obss.obs);
                    if first_ep {
                        for i in 0..obss.obs.len() {
                            for j in 0..NFREOBS {
                                if obss.obs[i].l[j] != 0.0 { obss.obs[i].lli[j] |= LLI_SLIP }
                            }
                        }
                        first_ep = false;
                    }
                    /* output obs */
                    let _ = outrnxobsb(writer, &opt, &obss.obs);
                    obs_flag = true;
                    /* update tr */
                    trt = obss.obs[0].time;
                    /* update */
                    pb.inc(1);
                    count += 1;
                    let message = format!("Decoding obs: {:.2}%", count as f32 / n as f32 * 100.0);
                    pb.set_message(message.to_string());
                }
                Ok(false) => continue
            }
        }
    }

    unsafe {
        (decoder.free_rtcm)(decoder.rtcm_ptr);
        libc::free(decoder.rtcm_ptr as *mut libc::c_void);
    }

    if !obs_flag && obsfile.is_some() { remove_file(obsfile.unwrap()).expect("Failed to delete file") }

    pb.finish_with_message("Finish decode");

    Ok(true)
}

/// Receive nrtip real-time data stream and save as rtcm file
pub async fn ntrip2rtcm(
    host: &str,
    port: &str,
    mountpoint: &str,
    username: &str,
    password: &str,
    outfile: &str,
) {
    let connect = conntrip(host, port, mountpoint, username, password).await;
    if connect.is_none() { return; }
    let mut stream = connect.unwrap();

    if let Some(parent) = Path::new(outfile).parent() {
        if !parent.exists() {
            create_dir_all(parent).expect("error in creating dir")
        }
    }

    let file = File::create(&outfile).expect("Fail to create output obs file!");
    let mut writer = BufWriter::new(file);

    // Create a channel for communication
    let (tx, mut rx) = mpsc::channel::<String>(100);
    // Start a task to listen for user input
    tokio::spawn(async move {
        let stdin = io::stdin();
        let reader = io::BufReader::new(stdin);
        let mut lines = reader.lines();
        while let Some(Ok(input)) = lines.next() {
            if tx.send(input.clone()).await.is_err() {
                break;
            }
            if input.trim() == "stop" {
                return;
            }
        }
    });

    // Processing NTRIP data streams
    let mut buf = vec![0; 2048];
    loop {
        tokio::select! {
            _ = async {
                let n=stream.read(&mut buf).await.expect("error in stream read");
                if n>0
                {
                    //println!("Received {} data in hex: {}", n,encode(&buf[0..n]));
                    writer.write_all(&buf[..n]).expect("Fail to store rtcm");
                    writer.flush().expect("Failed to flush writer");
                    writer.get_ref().sync_all().expect("Failed to sync data");
                    print!("\rStoring data at {}",timestr_rnx());
                    io::stdout().flush().expect("Flush failed");
                }
            } => {},

            // Processing user input
            Some(input) = rx.recv() => {
                match input.trim() {
                    "stop" => {
                        println!("Stopping...");
                        return;
                    },
                    _ => {
                        println!("Unknown command: {}", input.trim());
                    }
                }
            },
        }
    }
}

