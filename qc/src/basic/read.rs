use crate::basic::code::*;
use crate::basic::eph::*;
use crate::basic::sat::*;
use crate::basic::time::*;
use crate::basic::var::*;
use std::{
    cmp::Ordering,
    fs::File,
    io::{BufRead, BufReader, Lines},
    path::Path,
    str::{self, FromStr},
};
use indicatif::{ProgressBar, ProgressStyle};
/// satellite system codes
const SYSCODES: &str = "GREJCI";
/// observation type codes
const OBSCODES: &str = "CLDS";

/// read RINEX observation data
pub fn readrnxobs<R: BufRead>(
    lines: &mut Lines<R>,
    ts: GTime,
    te: GTime,
    tint: f64,
    opt: Option<&str>,
    rcv: usize,
    ver: f64,
    tsys: &mut i32,
    tobs: &mut [[[char; 4]; MAXOBSTYPE]; NUMSYS],
    obss: &mut Obss,
    sta: &mut Sta,
) -> Result<bool, String> {
    let mut data: Vec<Obs> = Vec::new();
    let mut slips = [[0u8; NFREQ + NEXOBS]; MAXSAT];
    let mut flag = 0;
    let mut stat = false;

    // 创建一个不显示具体进度的滚动进度条
    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner()
        .tick_chars("/|\\- ")
        .template("{spinner:.green} {msg}")
        .expect("Failed to set progress style"));
    pb.set_message("Reading obs...");

    if rcv > MAXRCV {
        return Ok(false);
    }

    while let Ok(n) = readrnxobsb(lines, opt, ver, tsys, tobs, &mut flag, &mut data, sta) {
        pb.tick();
        for i in 0..n as usize {
            if *tsys == TSYS_UTC {
                data[i].time = utc2gpst(data[i].time);
            }
            saveslips(&mut slips, &data[i]);
        }

        if n > 0 && !screent(data[0].time, ts, te, tint) {
            continue;
        }

        for i in 0..n as usize {
            restslips(&mut slips, &mut data[i]);
            data[i].rcv = rcv as u8;
            stat = obss.add_obs_data(&data[i].clone());
            if !stat {
                break;
            }
        }
        data = Vec::new();
    }
    pb.finish_with_message("Finish reading obs");
    return Ok(stat);
}

/// read RINEX header
pub fn read_rnxh<R: BufRead>(
    lines: &mut Lines<R>,
    ver: &mut f64,
    type_: &mut char,
    sys: &mut usize,
    tsys: &mut i32,
    tobs: &mut [[[char; 4]; MAXOBSTYPE]; NUMSYS],
    nav: &mut Nav,
    sta: &mut Sta,
) -> Result<bool, String> {
    let mut i = 0;

    *ver = 2.10;
    *type_ = ' ';
    *sys = SYS_GPS; // Assuming SYS_GPS is predefined somewhere

    while let Some(Ok(buff)) = lines.next() {
        let label = &buff[60..];

        if buff.len() <= 60 {
            continue;
        } else if label.contains("RINEX VERSION / TYPE") {
            *ver = buff[0..9].trim().parse().unwrap_or(0.0);
            *type_ = buff.chars().nth(20).unwrap_or(' ');

            match buff.chars().nth(40).unwrap_or(' ') {
                ' ' | 'G' => {
                    *sys = SYS_GPS;
                    *tsys = TSYS_GPS;
                }
                'R' => {
                    *sys = SYS_GLO;
                    *tsys = TSYS_UTC;
                }
                'E' => {
                    *sys = SYS_GAL;
                    *tsys = TSYS_GAL;
                }
                'J' => {
                    *sys = SYS_QZS;
                    *tsys = TSYS_QZS;
                }
                'C' => {
                    *sys = SYS_CMP;
                    *tsys = TSYS_CMP;
                }
                'I' => {
                    *sys = SYS_IRN;
                    *tsys = TSYS_IRN;
                }
                'M' => {
                    *sys = SYS_NONE;
                    *tsys = TSYS_GPS;
                }
                _ => {
                    return Err(String::from(format!(
                        "not supported satellite system: {}",
                        buff.chars().nth(40).unwrap_or(' ')
                    )));
                }
            }
            continue;
        } else if label.contains("PGM / RUN BY / DATE") || label.contains("COMMENT") {
            continue;
        }

        if *type_ == 'O' {
            decode_obsh(
                lines,
                &mut buff.clone(),
                *ver,
                tsys,
                tobs,
                &mut Some(nav),
                sta,
            );
        } else if *type_ == 'N' {
            decode_navh(buff.clone(), nav);
        }

        if label.contains("END OF HEADER") {
            return Ok(true);
        }

        if i >= MAXPOSHEAD && *type_ == ' ' {
            break; // No RINEX file
        }
        i += 1;
    }

    Ok(true)
}

/// read RINEX observation data body
pub fn readrnxobsb<R: BufRead>(
    lines: &mut Lines<R>,
    opt: Option<&str>,
    ver: f64,
    tsys: &mut i32,
    tobs: &mut [[[char; 4]; MAXOBSTYPE]; 6],
    flag: &mut i32,
    data: &mut Vec<Obs>,
    sta: &mut Sta,
) -> Result<i32, String> {
    let mut time: GTime = Default::default();
    let mut index = vec![Sigidx::default(); NUMSYS];
    let mut buff;
    let mut sats = vec![0; MAXOBS];
    let mut mask = 0;
    let mut nsat: i32 = 0;
    if let Some(value) = opt {
        mask = set_sysmask(value);
    }

    for i in 0..6 {
        if NUMSYS >= i + 1 {
            set_index(1 << i, opt, tobs[i], &mut index[i]);
        }
    }

    let mut n = 0;
    let mut i = 0;
    while let Some(Ok(line)) = lines.next() {
        buff = line;
        if i == 0 {
            if let Ok(x) = decode_obsepoch(lines, &mut buff, ver, &mut time, flag, &mut sats) {
                nsat = x;
                if nsat <= 0 {
                    continue;
                }
            }
        } else if (*flag <= 2 || *flag == 6) && n < MAXOBS {
            let mut obs_data = Obs {
                time: time.clone(),
                sat: sats[i - 1],
                rcv: 0,
                snr: [0; NFREQ + NEXOBS],
                lli: [0; NFREQ + NEXOBS],
                code: [0; NFREQ + NEXOBS],
                l: [0.0; NFREQ + NEXOBS],
                p: [0.0; NFREQ + NEXOBS],
                d: [0.0; NFREQ + NEXOBS],
            };
            // except SBS
            if buff.chars().next().unwrap_or(' ') == 'S' {
                nsat -= 1;
                if i > nsat as usize {
                    return Ok(n as i32);
                }
                continue;
            }
            if decode_obsdata(lines, &mut buff, ver, mask, &index, &mut obs_data)? {
                data.push(obs_data.clone());
                n += 1;
            }
        } else if *flag == 3 || *flag == 4 {
            decode_obsh(lines, &mut buff, ver, tsys, tobs, &mut None, sta);
        }
        if i >= nsat as usize {
            return Ok(n as i32);
        }
        i += 1;
    }
    Err(String::from("no obs"))
}

/// decode RINEX observation data file header
pub fn decode_obsh<R: BufRead>(
    lines: &mut Lines<R>,
    buff: &mut String,
    ver: f64,
    tsys: &mut i32,
    tobs: &mut [[[char; 4]; MAXOBSTYPE]; 6],
    nav: &mut Option<&mut Nav>,
    sta: &mut Sta,
) {
    const FRQ_CODES: &str = "1256789";
    const DEF_CODES: [&str; 7] = [
        "CWX    ", // GPS: L125____
        "CCXX X ", // GLO: L1234_6_
        "C XXXX ", // GAL: L1_5678_
        "CXXX   ", // QZS: L1256___
        "C X    ", // SBS: L1_5____
        "XIXIIX ", // BDS: L125678_
        "  A   A", // IRN: L__5___9
    ];

    let label = &buff[60..];

    if label.contains("MARKER NAME") {
        sta.name = buff[..60].trim().to_string();
    } else if label.contains("MARKER NUMBER") {
        sta.marker = buff[..20].trim().to_string();
    } else if label.contains("REC # / TYPE / VERS") {
        sta.recsno = buff[..20].trim().to_string();
        sta.rectype = buff[20..40].trim().to_string();
        sta.recver = buff[40..60].trim().to_string();
    } else if label.contains("ANT # / TYPE") {
        sta.antsno = buff[..20].trim().to_string();
        sta.antdes = buff[20..40].trim().to_string();
    } else if label.contains("APPROX POSITION XYZ") {
        for (i, chunk) in buff[..42].as_bytes().chunks(14).enumerate() {
            sta.pos[i] = str::from_utf8(chunk).unwrap().trim().parse().unwrap_or(0.0);
        }
    } else if label.contains("ANTENNA: DELTA H/E/N") {
        let mut del = [0.0; 3];
        for (i, chunk) in buff[..42].as_bytes().chunks(14).enumerate() {
            del[i] = str::from_utf8(chunk).unwrap().trim().parse().unwrap_or(0.0);
        }
        sta.del[2] = del[0]; // h
        sta.del[0] = del[1]; // e
        sta.del[1] = del[2]; // n
    } else if label.contains("SYS / # / OBS TYPES") {
        if let Some(p) = SYSCODES.find(buff.chars().next().unwrap()) {
            let i = p;
            let n: usize = buff[3..6].trim().parse().unwrap_or(0);
            let mut nt = 0;
            let mut k = 7;
            for _j in 0..n {
                if k > 58 {
                    if let Some(Ok(line)) = lines.next() {
                        *buff = line;
                        k = 7;
                    } else {
                        break;
                    }
                }
                if nt < MAXOBSTYPE - 1 {
                    tobs[i][nt] = [
                        buff.chars().nth(k).unwrap(),
                        buff.chars().nth(k + 1).unwrap(),
                        buff.chars().nth(k + 2).unwrap(),
                        '\0',
                    ];
                }
                k += 4;
                nt += 1;
            }
            tobs[i][nt][0] = '\0';

            if i == 5 && (ver - 3.02).abs() < 1e-3 {
                for obs in &mut tobs[i][..nt] {
                    if obs[1] == '1' {
                        obs[1] = '2';
                    }
                }
            }
            for obs in &mut tobs[i][..nt] {
                if !obs[2..3].is_empty() {
                    continue;
                }
                if let Some(p) = FRQ_CODES.find(obs[1]) {
                    obs[2] = DEF_CODES[i].chars().nth(p).unwrap();
                    println!(
                        "set default for unknown code: sys={} code={}",
                        buff.chars().nth(0).unwrap(),
                        obs.iter().collect::<String>()
                    );
                }
            }
        } else {
            let sys_name = buff.chars().next().unwrap();
            if sys_name != 'S' {
                println!("invalid system code: sys={}", sys_name);
                return;
            }
        }
    } else if label.contains("# / TYPES OF OBSERV") {
        let n: usize = buff[0..6].trim().parse().unwrap_or(0);
        let mut nt = 0;
        let mut j = 10;
        for _ in 0..n {
            if j > 58 {
                if let Some(Ok(line)) = lines.next() {
                    *buff = line;
                    j = 10;
                } else {
                    break;
                }
            }
            if nt >= MAXOBSTYPE - 1 {
                continue;
            }
            if ver <= 2.99 {
                let str = &buff[j..j + 2];
                convcode(ver, SYS_GPS, str, &mut tobs[0][nt]);
                convcode(ver, SYS_GLO, str, &mut tobs[1][nt]);
                convcode(ver, SYS_GAL, str, &mut tobs[2][nt]);
                convcode(ver, SYS_QZS, str, &mut tobs[3][nt]);
                convcode(ver, SYS_CMP, str, &mut tobs[5][nt]);
            }
            nt += 1;
            j += 6;
        }
        tobs[0][nt][0] = '\0';
    } else if label.contains("TIME OF FIRST OBS") {
        *tsys = match &buff[48..51] {
            "GPS" => 0,
            "GLO" => 1,
            "GAL" => 2,
            "QZS" => 3,
            "BDT" => 4,
            "IRN" => 5,
            _ => *tsys,
        };
    } else if label.contains("GLONASS SLOT / FRQ #") {
        for i in 0..8 {
            if &buff[4 + i * 7..5 + i * 7] != "R" {
                continue;
            }
            let prn: usize = buff[5 + i * 7..7 + i * 7].trim().parse().unwrap_or(0);
            let fcn: i32 = buff[8 + i * 7..10 + i * 7].trim().parse().unwrap_or(0);
            if prn >= 1 && prn <= 27 && fcn >= -7 && fcn <= 6 {
                if let Some(nav_ref) = nav {
                    nav_ref.glo_fcn[prn - 1] = fcn + 8;
                }
            }
        }
    } else if label.contains("GLONASS COD/PHS/BIS") {
        sta.glo_cp_bias[0] = buff[5..13].trim().parse().unwrap_or(0.0);
        sta.glo_cp_bias[1] = buff[18..26].trim().parse().unwrap_or(0.0);
        sta.glo_cp_bias[2] = buff[31..39].trim().parse().unwrap_or(0.0);
        sta.glo_cp_bias[3] = buff[44..52].trim().parse().unwrap_or(0.0);
    } else if label.contains("LEAP SECONDS") {
        if let Some(ref mut nav) = nav {
            nav.utc_gps[4] = buff[0..6].trim().parse().unwrap_or(0.0);
            nav.utc_gps[7] = buff[6..12].trim().parse().unwrap_or(0.0);
            nav.utc_gps[5] = buff[12..18].trim().parse().unwrap_or(0.0);
            nav.utc_gps[6] = buff[18..24].trim().parse().unwrap_or(0.0);
        }
    }
}

/// decode observation epoch
pub fn decode_obsepoch<R: BufRead>(
    lines: &mut Lines<R>,
    buff: &mut String,
    ver: f64,
    time: &mut GTime,
    flag: &mut i32,
    sats: &mut [usize],
) -> Result<i32, String> {
    let n: i32;

    if ver <= 2.99 {
        // Version 2.x processing
        n = buff[29..32].trim().parse().unwrap_or(0);
        if n <= 0 {
            return Err(String::from("0 sat in this epoch"));
        }

        // epoch flag: 3:new site,4:header info,5:external event
        *flag = buff[28..29].trim().parse().unwrap_or(0);
        if (3..=5).contains(flag) {
            return Ok(n);
        }

        if str2time(&buff[0..26], time).is_err() {
            return Err(format!("rinex obs invalid epoch: epoch={}", &buff[0..26]));
        }

        let mut j = 32;
        for i in 0..n as usize {
            if j >= 68 {
                buff.clear();
                if let Some(Ok(line)) = lines.next() {
                    *buff = line;
                } else {
                    break;
                }
                j = 32;
            }
            if i < MAXOBS {
                let satid = &buff[j..j + 3];
                sats[i] = satid2no(satid);
            }
            j += 3;
        }
    } else {
        // Version 3.x processing
        n = buff[32..35].trim().parse().unwrap_or(0);
        if n <= 0 {
            return Err(String::from("0 sat in this epoch"));
        }

        *flag = buff[31..32].trim().parse().unwrap_or(0);
        if (3..=5).contains(flag) {
            return Ok(n);
        }

        if buff.chars().nth(0) != Some('>') || str2time(&buff[1..29], time).is_err() {
            return Err(format!("rinex obs invalid epoch: epoch={}", &buff[1..30]));
        }
    }

    //return Err(format!("decode_obsepoch: time={} flag={}", time2str(*time), flag));

    Ok(n)
}

/// decode observation data
pub fn decode_obsdata<R: BufRead>(
    lines: &mut Lines<R>,
    buff: &mut String,
    ver: f64,
    mask: usize,
    index: &Vec<Sigidx>,
    obs: &mut Obs,
) -> Result<bool, String> {
    let mut val = [0.0; MAXOBSTYPE];
    let mut lli = [0u8; MAXOBSTYPE];
    let mut satid = String::new();
    let stat = 1;
    let mut p = [0; MAXOBSTYPE];
    let mut k = [0; 16];
    let mut l = [0; 16];
    let mut n = 0;
    let mut m = 0;

    //,println!("decode_obsdata: ver={:.2}", ver);

    if ver > 2.99 {
        satid = buff.chars().take(3).collect();
        obs.sat = satid2no(&satid);
    }
    if obs.sat == 0 {
        return Err(format!("decode_obsdata: unsupported sat sat={}", satid));
    } else if (satsys(obs.sat) & mask) == 0 {
        return Ok(false);
    }

    let ind = match satsys(obs.sat) {
        SYS_GLO => &index[1],
        SYS_GAL => &index[2],
        SYS_QZS => &index[3],
        SYS_CMP => &index[4],
        SYS_IRN => &index[5],
        _ => &index[0],
    };

    let mut j = if ver <= 2.99 { 0 } else { 3 };
    for i in 0..ind.n {
        if ver <= 2.99 && j >= 80 {
            if let Some(Ok(line)) = lines.next() {
                *buff = line;
                j = 0;
            } else {
                break;
            }
        }
        if stat != 0 && j < buff.len() {
            if j + 14 > buff.len() { break; }
            val[i] = f64::from_str(&buff[j..j + 14].trim()).unwrap_or(0.0) + ind.shift[i];
            if j + 15 > buff.len() { break; }
            lli[i] = u8::from_str(&buff[j + 14..j + 15].trim()).unwrap_or(0) & 3;
        }
        j += 16;
    }
    if stat == 0 {
        return Ok(false);
    }

    for i in 0..NFREQ + NEXOBS {
        obs.p[i] = 0.0;
        obs.l[i] = 0.0;
        obs.d[i] = 0.0;
        obs.snr[i] = 0;
        obs.lli[i] = 0;
        obs.code[i] = 0;
    }

    for i in 0..ind.n {
        p[i] = if ver <= 2.11 { ind.idx[i] } else { ind.pos[i] };

        if ind.type_[i] == 0 && p[i] == 0 {
            k[n] = i;
            n += 1;
        }
        if ind.type_[i] == 0 && p[i] == 1 {
            l[m] = i;
            m += 1;
        }
    }

    if ver <= 2.11 {
        if n >= 2 {
            if val[k[0]] == 0.0 && val[k[1]] == 0.0 {
                p[k[0]] = -1;
                p[k[1]] = -1;
            } else if val[k[0]] != 0.0 && val[k[1]] == 0.0 {
                p[k[0]] = 0;
                p[k[1]] = -1;
            } else if val[k[0]] == 0.0 && val[k[1]] != 0.0 {
                p[k[0]] = -1;
                p[k[1]] = 0;
            } else if ind.pri[k[1]] > ind.pri[k[0]] {
                p[k[1]] = 0;
                p[k[0]] = if NEXOBS < 1 { -1 } else { NFREQ as i32 };
            } else {
                p[k[0]] = 0;
                p[k[1]] = if NEXOBS < 1 { -1 } else { NFREQ as i32 };
            }
        }
        if m >= 2 {
            if val[l[0]] == 0.0 && val[l[1]] == 0.0 {
                p[l[0]] = -1;
                p[l[1]] = -1;
            } else if val[l[0]] != 0.0 && val[l[1]] == 0.0 {
                p[l[0]] = 1;
                p[l[1]] = -1;
            } else if val[l[0]] == 0.0 && val[l[1]] != 0.0 {
                p[l[0]] = -1;
                p[l[1]] = 1;
            } else if ind.pri[l[1]] > ind.pri[l[0]] {
                p[l[1]] = 1;
                p[l[0]] = if NEXOBS < 2 { -1 } else { (NFREQ + 1) as i32 };
            } else {
                p[l[0]] = 1;
                p[l[1]] = if NEXOBS < 2 { -1 } else { (NFREQ + 1) as i32 };
            }
        }
    }

    for i in 0..ind.n {
        if p[i] < 0 || val[i] == 0.0 {
            continue;
        }
        match ind.type_[i] {
            0 => {
                obs.p[p[i] as usize] = val[i];
                obs.code[p[i] as usize] = ind.code[i];
            }
            1 => {
                obs.l[p[i] as usize] = val[i];
                obs.lli[p[i] as usize] = lli[i] as i32;
            }
            2 => {
                obs.d[p[i] as usize] = val[i];
            }
            3 => {
                obs.snr[p[i] as usize] = (val[i] / SNR_UNIT + 0.5) as i32;
            }
            _ => {}
        }
    }
    Ok(true)
}

/// convert RINEX obs-type ver.2 -> ver.3
pub fn convcode(ver: f64, sys: usize, str: &str, type_: &mut [char; 4]) {
    let mut tob: String = String::from("   ");

    if str == "P1" { // ver.2.11 GPS L1PY,GLO L2P
        if sys == SYS_GPS {
            tob = format!("{}1W", 'C');
        } else if sys == SYS_GLO {
            tob = format!("{}1P", 'C');
        } else if sys == SYS_CMP {
            tob = format!("{}2I", 'C');
        }
    } else if str == "P2" { // ver.2.11 GPS L2PY,GLO L2P
        if sys == SYS_GPS {
            tob = format!("{}2W", 'C');
        } else if sys == SYS_GLO {
            tob = format!("{}2P", 'C');
        } else if sys == SYS_CMP {
            tob = format!("{}7I", 'C');
        }
    } else if str == "C1" { // ver.2.11 GPS L1C,GLO L1C/A
        if ver < 2.12 {
            if sys == SYS_GPS {
                tob = format!("{}1C", 'C');
            } else if sys == SYS_GLO {
                tob = format!("{}1C", 'C');
            } else if sys == SYS_GAL {
                tob = format!("{}1X", 'C'); // ver.2.12
            } else if sys == SYS_QZS {
                tob = format!("{}1C", 'C');
            } else if sys == SYS_CMP {
                tob = format!("{}2I", 'C');
            }
        }
    } else if str == "C2" {
        if sys == SYS_GPS {
            if ver >= 2.12 {
                tob = format!("{}2W", 'C'); // L2P(Y)
            } else {
                tob = format!("{}2X", 'C'); // L2C
            }
        } else if sys == SYS_GLO {
            tob = format!("{}2C", 'C');
        } else if sys == SYS_QZS {
            tob = format!("{}2X", 'C');
        } else if sys == SYS_CMP {
            tob = format!("{}7I", 'C'); // ver.2.12 B1_2
        }
    } else if ver >= 2.12 && str.chars().nth(1) == Some('A') { // ver.2.12 L1C/A
        if sys == SYS_GPS || sys == SYS_GLO || sys == SYS_QZS {
            tob = format!("{}1C", str.chars().nth(0).unwrap());
        }
    } else if ver >= 2.12 && str.chars().nth(1) == Some('B') { // ver.2.12 GPS L1C
        if sys == SYS_GPS || sys == SYS_QZS {
            tob = format!("{}1X", str.chars().nth(0).unwrap());
        }
    } else if ver >= 2.12 && str.chars().nth(1) == Some('C') { // ver.2.12 GPS L2C
        if sys == SYS_GPS || sys == SYS_QZS {
            tob = format!("{}2X", str.chars().nth(0).unwrap());
        }
    } else if ver >= 2.12 && str.chars().nth(1) == Some('D') { // ver.2.12 GLO L2C/A
        if sys == SYS_GLO {
            tob = format!("{}2C", str.chars().nth(0).unwrap());
        }
    } else if ver >= 2.12 && str.chars().nth(1) == Some('1') { // ver.2.12 GPS L1PY,GLO L1P
        if sys == SYS_GPS {
            tob = format!("{}1W", str.chars().nth(0).unwrap());
        } else if sys == SYS_GLO {
            tob = format!("{}1P", str.chars().nth(0).unwrap());
        } else if sys == SYS_GAL {
            tob = format!("{}1X", str.chars().nth(0).unwrap());
        } else if sys == SYS_CMP {
            tob = format!("{}2X", str.chars().nth(0).unwrap());
        }
    } else if ver < 2.12 && str.chars().nth(1) == Some('1') {
        if sys == SYS_GPS || sys == SYS_GLO || sys == SYS_QZS {
            tob = format!("{}1C", str.chars().nth(0).unwrap());
        } else if sys == SYS_GAL {
            tob = format!("{}1X", str.chars().nth(0).unwrap());
        } else if sys == SYS_CMP {
            tob = format!("{}2I", str.chars().nth(0).unwrap());
        }
    } else if str.chars().nth(1) == Some('2') {
        if sys == SYS_GPS {
            tob = format!("{}2W", str.chars().nth(0).unwrap());
        } else if sys == SYS_GLO {
            tob = format!("{}2P", str.chars().nth(0).unwrap());
        } else if sys == SYS_QZS {
            tob = format!("{}2X", str.chars().nth(0).unwrap());
        } else if sys == SYS_CMP {
            tob = format!("{}7I", str.chars().nth(0).unwrap()); // ver.2.12 B1_2
        }
    } else if str.chars().nth(1) == Some('5') {
        if sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS {
            tob = format!("{}5X", str.chars().nth(0).unwrap());
        }
    } else if str.chars().nth(1) == Some('6') {
        if sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP {
            tob = format!("{}6X", str.chars().nth(0).unwrap());
        }
    } else if str.chars().nth(1) == Some('7') {
        if sys == SYS_GAL || sys == SYS_CMP {
            tob = format!("{}7X", str.chars().nth(0).unwrap());
        }
    } else if str.chars().nth(1) == Some('8') {
        if sys == SYS_GAL {
            tob = format!("{}8X", str.chars().nth(0).unwrap());
        }
    }

    let chars: Vec<char> = tob.chars().collect();
    for (i, c) in chars.iter().enumerate() {
        type_[i] = *c;
    }
}

/// decode RINEX NAV header
fn decode_navh(buff: String, nav: &mut Nav) {
    let label = &buff[60..];

    match label {
        l if l.starts_with("ION ALPHA") => { // opt ver.2
            for i in 0..4 {
                nav.ion_gps[i] = buff[2 + i * 12..2 + (i + 1) * 12]
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            }
        }
        l if l.starts_with("ION BETA") => { // opt ver.2
            for i in 0..4 {
                nav.ion_gps[4 + i] = buff[2 + i * 12..2 + (i + 1) * 12]
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            }
        }
        l if l.starts_with("DELTA-UTC: A0,A1,T,W") => { // opt ver.2
            for i in 0..2 {
                nav.utc_gps[i] = buff[3 + i * 19..3 + (i + 1) * 19]
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            }
            for i in 2..4 {
                nav.utc_gps[i] = buff[3 + i * 9..3 + (i + 1) * 9]
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            }
        }
        l if l.starts_with("IONOSPHERIC CORR") => match &buff[..4] {  // opt ver.3
            "GPSA" | "GPSB" => {
                let offset = if buff.starts_with("GPSA") { 0 } else { 4 };
                for i in 0..4 {
                    nav.ion_gps[offset + i] = buff[5 + i * 12..5 + (i + 1) * 12]
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                }
            }
            "GAL " => {
                for i in 0..4 {
                    nav.ion_gal[i] = buff[5 + i * 12..5 + (i + 1) * 12]
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                }
            }
            "QZSA" | "QZSB" => { // v.3.02
                let offset = if buff.starts_with("QZSA") { 0 } else { 4 };
                for i in 0..4 {
                    nav.ion_qzs[offset + i] = buff[5 + i * 12..5 + (i + 1) * 12]
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                }
            }
            "BDSA" | "BDSB" => { // v.3.02
                let offset = if buff.starts_with("BDSA") { 0 } else { 4 };
                for i in 0..4 {
                    nav.ion_cmp[offset + i] = buff[5 + i * 12..5 + (i + 1) * 12]
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                }
            }
            "IRNA" | "IRNB" => { // v.3.03
                let offset = if buff.starts_with("IRNA") { 0 } else { 4 };
                for i in 0..4 {
                    nav.ion_irn[offset + i] = buff[5 + i * 12..5 + (i + 1) * 12]
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                }
            }
            _ => {}
        },
        l if l.starts_with("TIME SYSTEM CORR") => { // opt ver.3
            match &buff[..4] {
                "GPUT" => {
                    nav.utc_gps[0] = buff[5..22].trim().parse().unwrap_or(0.0);
                    nav.utc_gps[1] = buff[22..38].trim().parse().unwrap_or(0.0);
                    nav.utc_gps[2] = buff[38..45].trim().parse().unwrap_or(0.0);
                    nav.utc_gps[3] = buff[45..50].trim().parse().unwrap_or(0.0);
                }
                "GLUT" => { // tau_C
                    nav.utc_glo[0] = -buff[5..22].trim().parse().unwrap_or(0.0);
                    // tau_C
                }
                "GLGP" => { // tau_GPS
                    nav.utc_glo[1] = buff[5..22].trim().parse().unwrap_or(0.0); // tau_GPS
                }
                "GAUT" => { // v.3.02
                    nav.utc_gal[0] = buff[5..22].trim().parse().unwrap_or(0.0);
                    nav.utc_gal[1] = buff[22..38].trim().parse().unwrap_or(0.0);
                    nav.utc_gal[2] = buff[38..45].trim().parse().unwrap_or(0.0);
                    nav.utc_gal[3] = buff[45..50].trim().parse().unwrap_or(0.0);
                }
                "QZUT" => { // v.3.02
                    nav.utc_qzs[0] = buff[5..22].trim().parse().unwrap_or(0.0);
                    nav.utc_qzs[1] = buff[22..38].trim().parse().unwrap_or(0.0);
                    nav.utc_qzs[2] = buff[38..45].trim().parse().unwrap_or(0.0);
                    nav.utc_qzs[3] = buff[45..50].trim().parse().unwrap_or(0.0);
                }
                "BDUT" => { // v.3.02
                    nav.utc_cmp[0] = buff[5..22].trim().parse().unwrap_or(0.0);
                    nav.utc_cmp[1] = buff[22..38].trim().parse().unwrap_or(0.0);
                    nav.utc_cmp[2] = buff[38..45].trim().parse().unwrap_or(0.0);
                    nav.utc_cmp[3] = buff[45..50].trim().parse().unwrap_or(0.0);
                }
                "IRUT" => { // v.3.03
                    nav.utc_irn[0] = buff[5..22].trim().parse().unwrap_or(0.0);
                    nav.utc_irn[1] = buff[22..38].trim().parse().unwrap_or(0.0);
                    nav.utc_irn[2] = buff[38..45].trim().parse().unwrap_or(0.0);
                    nav.utc_irn[3] = buff[45..50].trim().parse().unwrap_or(0.0);
                    nav.utc_irn[8] = 0.0;
                }
                _ => {}
            }
        }
        l if l.starts_with("LEAP SECONDS") => {
            nav.utc_gps[4] = buff[0..6].trim().parse().unwrap_or(0.0);
            nav.utc_gps[7] = buff[6..12].trim().parse().unwrap_or(0.0);
            nav.utc_gps[5] = buff[12..18].trim().parse().unwrap_or(0.0);
            nav.utc_gps[6] = buff[18..24].trim().parse().unwrap_or(0.0);
        }
        _ => {}
    }
}

// set system mask
pub fn set_sysmask(opt: &str) -> usize {
    if let Some(pos) = opt.find("-SYS=") {
        let mut mask = SYS_NONE;
        let mut chars = opt[pos + 5..].chars();

        while let Some(c) = chars.next() {
            match c {
                'G' => mask |= SYS_GPS,
                'R' => mask |= SYS_GLO,
                'E' => mask |= SYS_GAL,
                'J' => mask |= SYS_QZS,
                'C' => mask |= SYS_CMP,
                'I' => mask |= SYS_IRN,
                ' ' => break,
                _ => {}
            }
        }
        mask
    } else {
        SYS_ALL
    }
}

// set signal index
pub fn set_index(sys: usize, opt: Option<&str>, tobs: [[char; 4]; MAXOBSTYPE], ind: &mut Sigidx) {
    let mut n = 0;

    for (i, obs_str) in tobs.iter().enumerate() {
        let obs_string: String = obs_str.iter().collect();
        let obs_trimmed: &str = obs_string.trim_end_matches('\0');
        let obs = match obs_trimmed {
            "" => continue, // 如果转换后是空字符串，则跳过
            valid_str => valid_str,
        };
        if obs.is_empty() {
            break;
        }

        ind.code[i] = obs2code(&obs[1..]);
        ind.type_[i] = OBSCODES
            .find(obs.chars().next().unwrap())
            .map_or(0, |p| p as i32);
        ind.idx[i] = code2idx(sys, ind.code[i]);
        ind.pri[i] = getcodepri(sys, ind.code[i], opt);
        ind.pos[i] = -1;
        n += 1;
    }

    let optstr = match sys {
        SYS_GPS => "-GL",
        SYS_GLO => "-RL",
        SYS_GAL => "-EL",
        SYS_QZS => "-JL",
        SYS_CMP => "-CL",
        SYS_IRN => "-IL",
        _ => "",
    };

    if let Some(opt) = opt {
        let opt_parts: Vec<&str> = opt.split('-').collect();
        for p in opt_parts.iter() {
            if p.is_empty() {
                continue;
            }

            let shift_str;
            if p.starts_with(optstr) {
                let remaining = &p[optstr.len()..];
                let parts: Vec<&str> = remaining.split('=').collect();
                if parts.len() == 2 {
                    let str_part = parts[0];
                    let shift: f64 = match parts[1].parse() {
                        Ok(val) => val,
                        Err(_) => continue,
                    };
                    shift_str = str_part.to_string();
                    for i in 0..n {
                        if code2obs(ind.code[i]) == shift_str {
                            ind.shift[i] = shift;
                        }
                    }
                }
            }
        }
    }

    for i in 0..NFREQ {
        let mut k: i32 = -1;
        for j in 0..n {
            if ind.idx[j] == i as i32
                && ind.pri[j] != 0
                && (k < 0 || ind.pri[j] > ind.pri[k as usize])
            {
                k = j as i32;
            }
        }
        if k < 0 {
            continue;
        }

        for j in 0..n {
            if ind.code[j] == ind.code[k as usize] {
                ind.pos[j] = i as i32;
            }
        }
    }

    for i in 0..NEXOBS {
        let mut j: usize = 0;

        while j < n {
            if ind.code[j] != 0 && ind.pri[j] != 0 && ind.pos[j] < 0 {
                break;
            }
            j += 1;
        }
        if j >= n {
            break;
        }

        for k in 0..n {
            if ind.code[k] == ind.code[j] {
                ind.pos[k] = NFREQ as i32 + i as i32;
            }
        }
    }

    for i in 0..n {
        if ind.code[i] == 0 || ind.pri[i] == 0 || ind.pos[i] >= 0 {
            continue;
        }
    }

    ind.n = n;
}

// save cycle slips
pub fn saveslips(slips: &mut [[u8; NFREQ + NEXOBS]], data: &Obs) {
    for i in 0..NFREQ + NEXOBS {
        if data.lli[i] & 1 != 0 {
            slips[data.sat - 1][i] |= LLI_SLIP as u8;
        }
    }
}

// restore cycle slips
pub fn restslips(slips: &mut [[u8; NFREQ + NEXOBS]], data: &mut Obs) {
    for i in 0..NFREQ + NEXOBS {
        if slips[data.sat - 1][i] & 1 != 0 {
            data.lli[i] |= LLI_SLIP;
        }
        slips[data.sat - 1][i] = 0;
    }
}

/// read RINEX navigation data body
pub fn readrnxnavb<R: BufRead>(
    lines: &mut Lines<R>,
    opt: Option<&str>,
    ver: f64,
    mut sys: usize,
    type_: &mut i32,
    eph: &mut Eph,
    geph: &mut Geph,
) -> Result<bool, String> {
    let mut toc = GTime::default();
    let mut data = [0.0; 64];
    let mut i = 0;
    let mut prn: i32;
    let mut sat = 0;
    let mut sp = 3;
    let mut id;
    let mut mask = 0;
    if let Some(value) = opt {
        mask = set_sysmask(value);
    }

    while let Some(Ok(line)) = lines.next() {
        let buff = line;
        if i == 0 {
            // if !SYSCODES.contains(buff.chars().nth(0).unwrap()) { continue }
            if ver >= 3.0 || sys == SYS_GAL || sys == SYS_QZS {
                id = buff[0..3].to_string();
                sat = satid2no(&id);
                sp = 4;
                if ver >= 3.0 {
                    sys = satsys(sat);
                    if sys == 0 {
                        sys = match id.chars().nth(0).unwrap() {
                            'R' => SYS_GLO,
                            'S' => {
                                for _ in 0..3 {
                                    if let Some(result) = lines.next() {
                                        result.map_err(|e| format!("Error reading line: {}", e))?;
                                    }
                                }
                                continue;
                            }
                            _ => SYS_GPS,
                        };
                    }
                }
            } else {
                prn = i32::from_str(&buff[0..2].trim()).map_err(|e| e.to_string())?;
                sat = if sys == SYS_GLO {
                    satno(SYS_GLO, prn as usize)
                } else if (93..=97).contains(&prn) {
                    satno(SYS_QZS, (prn + 100) as usize)
                } else {
                    satno(SYS_GPS, prn as usize)
                };
            }

            str2time(&buff[sp..sp + 19], &mut toc).map_err(|e| e.to_string())?;

            for j in 0..3 {
                let value = &buff[sp + (j + 1) * 19..sp + (j + 2) * 19].trim();
                data[i] = value
                    .replace('D', "E")
                    .parse::<f64>()
                    .map_err(|e| e.to_string())?;
                i += 1;
            }
        } else {
            // Decoding data fields
            for j in 0..4 {
                if sp + (j + 1) * 19 > buff.len() {
                    i += 1;
                    continue;
                }
                let value = &buff[sp + j * 19..sp + (j + 1) * 19].trim();
                if !value.is_empty() {
                    data[i] = value
                        .replace('D', "E")
                        .parse::<f64>()
                        .map_err(|e| e.to_string())?;
                }
                i += 1;
            }

            // decode the ephemeris
            if sys == SYS_GLO && i >= 15 {
                if (mask & sys) == 0 {
                    return Err(String::from("geph failed"));
                }
                *type_ = 1;
                return decode_geph(ver, sat, toc, &data, geph);
            } else if i >= 31 {
                if (mask & sys) == 0 {
                    return Err(String::from("eph failed"));
                }
                *type_ = 0;
                return decode_eph(sat, toc, &data, eph);
            }
        }
    }
    Ok(false)
}

/// read RINEX navigation data
pub fn readrnxnav<R: BufRead>(
    lines: &mut Lines<R>,
    opt: Option<&str>,
    ver: f64,
    sys: usize,
    nav: &mut Nav,
) -> Result<bool, String> {
    let mut eph = Default::default();
    let mut geph = Default::default();
    let mut type_ = 0;

    // Create a scrolling progress bar
    let pb = ProgressBar::new_spinner();
    // Setting the Progress Bar Style
    pb.set_style(ProgressStyle::default_spinner()
        .tick_chars("/|\\- ")
        .template("{spinner:.green} {msg}")
        .expect("Failed to set progress style"));
    pb.set_message("Reading nav...");

    while let Some(stat) = match readrnxnavb(lines, opt, ver, sys, &mut type_, &mut eph, &mut geph) {
        Ok(s) => Some(s),
        Err(e) => return Err(e),
    } {
        pb.tick();
        if stat {
            let add_stat = match type_ {
                1 => nav.add_geph(geph.clone()),
                _ => nav.add_eph(eph.clone()),
            };
            if !add_stat {
                return Ok(false);
            }
        } else {
            break;
        }
    }
    pb.finish_with_message("Finish reading nav");
    return Ok(nav.eph.len() > 0 || nav.geph.len() > 0);
}

/// read RINEX file
pub fn readrnxfile(
    file: &String,
    ts: GTime,
    te: GTime,
    tint: f64,
    opt: Option<&str>,
    flag: i32,
    index: usize,
    obss: &mut Obss,
    nav: &mut Nav,
    sta: &mut Sta,
) -> Result<bool, String> {
    let path = Path::new(&file);
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut lines = BufReader::new(file).lines();
    let mut ver = 0.0;
    let mut sys = 0;
    let mut tsys = TSYS_GPS;
    let mut tobs = [[['\0'; 4]; MAXOBSTYPE]; NUMSYS];
    let mut type_ = '\0';
    if !matches!(
        read_rnxh(&mut lines, &mut ver, &mut type_, &mut sys, &mut tsys, &mut tobs, nav, sta),
        Ok(true)
    ) {
        return Ok(false);
    }
    if (flag == 0 && type_ == 'C') || (flag != 0 && type_ != 'C') {
        return Ok(false);
    }
    match type_ {
        'O' => {
            return readrnxobs(
                &mut lines, ts, te, tint, opt, index, ver, &mut tsys, &mut tobs, obss, sta,
            )
        }
        'N' => return readrnxnav(&mut lines, opt, ver, sys, nav),
        'G' => return readrnxnav(&mut lines, opt, ver, SYS_GLO, nav),
        'J' => return readrnxnav(&mut lines, opt, ver, SYS_QZS, nav),
        'L' => return readrnxnav(&mut lines, opt, ver, SYS_GAL, nav),
        _ => Err(String::from("file type error")),
    }
}

/// read obs and nav data
pub fn readobsnav(
    ts: GTime,
    te: GTime,
    tint: f64,
    ifile: Vec<String>,
    index: [i32; MAXINFILE],
    //prcopt
    obss: &mut Obss,
    nav: &mut Nav,
    sta: &mut Sta,
    nepoch: &mut usize,
) -> Result<bool, String> {
    let n = ifile.len();
    let mut ind = 0;
    let mut nobs = 0;
    let mut rcv = 1;
    let opt = Some("");

    for i in 0..n {
        if index[i] != ind {
            if obss.n > nobs {
                rcv += 1;
            }
            ind = index[i];
            nobs = obss.n;
        }
        readrnxfile(&ifile[i], ts, te, tint, opt, 0, rcv, obss, nav, sta)?;
    }
    *nepoch = sortobs(obss);
    uniqnav(nav);
    Ok(true)
}

//  compare observation data
fn cmpobs(a: &Obs, b: &Obs) -> Ordering {
    let tt = timediff(a.time, b.time);
    if tt.abs() > DTTOL {
        if tt < 0.0 {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    } else if a.rcv != b.rcv {
        a.rcv.cmp(&b.rcv)
    } else {
        a.sat.cmp(&b.sat)
    }
}

/// sort and unique observation data
pub fn sortobs(obss: &mut Obss) -> usize {
    if obss.n <= 0 {
        return 0;
    }

    obss.obs.sort_by(cmpobs);

    // 删除重复数据
    let mut j = 0;
    for i in 0..obss.n {
        if obss.obs[i].sat != obss.obs[j].sat
            || obss.obs[i].rcv != obss.obs[j].rcv
            || timediff(obss.obs[i].time, obss.obs[j].time) != 0.0
        {
            j += 1;
            obss.obs[j] = obss.obs[i].clone();
        }
    }
    obss.n = j + 1;
    obss.obs.truncate(obss.n); // Make sure the length of vec is consistent with obs.n

    let mut n = 0;
    let mut i = 0;
    while i < obss.n {
        let mut j = i + 1;
        while j < obss.n {
            if timediff(obss.obs[j].time, obss.obs[i].time) > DTTOL {
                break;
            }
            j += 1;
        }
        i = j;
        n += 1;
    }
    n
}

// compare ephemeris
fn cmpeph(a: &Eph, b: &Eph) -> Ordering {
    if a.ttr.time != b.ttr.time {
        return a.ttr.time.partial_cmp(&b.ttr.time).unwrap();
    }
    if a.toe.time != b.toe.time {
        return a.toe.time.partial_cmp(&b.toe.time).unwrap();
    }
    if a.sat == b.sat && satsys(a.sat) == SYS_GAL {
        let bit_a = (a.code & (1 << 9)) != 0;
        let bit_b = (b.code & (1 << 9)) != 0;
        if bit_a != bit_b {
            return bit_b.cmp(&bit_a); // bit_b comes first, so bit_b returns Greater if it is true.
        }
    }
    a.sat.cmp(&b.sat)
}

// sort and unique ephemeris
fn uniqeph(nav: &mut Nav) {
    if nav.n <= 0 {
        return;
    }

    nav.eph.sort_by(cmpeph);

    // Deletion of duplicate data
    let mut j = 0;
    for i in 1..nav.n {
        if nav.eph[i].sat != nav.eph[j].sat || nav.eph[i].iode != nav.eph[j].iode {
            j += 1;
            nav.eph[j] = nav.eph[i].clone();
        }
    }
    nav.n = j + 1;
    nav.eph.truncate(nav.n);

    // 确Ensure that the capacity and length of the Vec matches
    nav.eph.shrink_to_fit();
}

// compare glonass ephemeris
fn cmpgeph(a: &Geph, b: &Geph) -> Ordering {
    if a.tof.time != b.tof.time {
        return a.tof.time.partial_cmp(&b.tof.time).unwrap();
    }
    if a.toe.time != b.toe.time {
        return a.toe.time.partial_cmp(&b.toe.time).unwrap();
    }
    a.sat.cmp(&b.sat)
}

// sort and unique glonass ephemeris
fn uniqgeph(nav: &mut Nav) {
    if nav.ng <= 0 {
        return;
    }

    nav.geph.sort_by(cmpgeph);

    // 删除重复数据
    let mut j = 0;
    for i in 1..nav.ng {
        if nav.geph[i].sat != nav.geph[j].sat
            || nav.geph[i].toe.time != nav.geph[j].toe.time
            || nav.geph[i].svh != nav.geph[j].svh
        {
            j += 1;
            nav.geph[j] = nav.geph[i].clone();
        }
    }
    nav.ng = j + 1;
    nav.geph.truncate(nav.ng);

    // 确保 Vec 的容量和长度匹配
    nav.geph.shrink_to_fit();
}

/// sort and unique ephemeris
pub fn uniqnav(nav: &mut Nav) {
    uniqeph(nav);
    uniqgeph(nav);
}
