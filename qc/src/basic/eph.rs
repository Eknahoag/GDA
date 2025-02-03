use nalgebra::DVector;
use crate::basic::func::sqr;
use crate::basic::sat::*;
use crate::basic::time::*;
use crate::basic::var::*;
use nalgebra::DMatrix;
//use tracy_client::span;

/// RAa values
const URA_EPH: [f64; 16] = [
    2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0,
    6144.0, 0.0,
];

/// max time difference to GPS Toe (s)
const MAXDTOE: f64 = 7200.0;
/// max time difference to QZS Toe (s)
const MAXDTOE_QZS: f64 = 7200.0;
/// max time difference to GAL Toe (s)
const MAXDTOE_GAL: f64 = 7200.0;
/// max time difference to BDS Toe (s)
const MAXDTOE_CMP: f64 = 7200.0;
/// max time difference to GLO Toe (s)
const MAXDTOE_GLO: f64 = 7200.0;
/// max time difference to IRN Toe (s)
const MAXDTOE_IRN: f64 = 7200.0;
/// error of broadcast clock (m)
const STD_BRDCCLK: f64 = 30.0;
/// error of galileo ephemeris for NAPA (m)
const STD_GAL_NAPA: f64 = 500.0;
/// error of glonass ephemeris (m)
const ERREPH_GLO: f64 = 5.0;
/// integration step glonass ephemeris (s)
const TSTEP: f64 = 60.0;
/// relative tolerance for Kepler equation
const RTOL_KEPLER: f64 = 1E-13;
/// max number of iteration of Kelpler
const MAX_ITER_KEPLER: usize = 30;
/// radius of earth (m)
const RE_GLO: f64 = 6378136.0;
/// gravitational constant
const MU_GPS: f64 = 3.9860050E14;
/// gravitational constant
const MU_GLO: f64 = 3.9860044E14;
/// gravitational constant
const MU_GAL: f64 = 3.986004418E14;
/// gravitational constant
const MU_CMP: f64 = 3.986004418E14;
/// 2nd zonal harmonic of geopot
const J2_GLO: f64 = 1.0826257E-3;
/// earth angular velocity (rad/s)
const OMGE_GLO: f64 = 7.292115E-5;
/// earth angular velocity (rad/s)
const OMGE_GAL: f64 = 7.2921151467E-5;
/// earth angular velocity (rad/s)
const OMGE_CMP: f64 = 7.292115E-5;
/// cos(-5.0 deg)
const COS_5: f64 = 0.9961946980917456;
/// sin(-5.0 deg)
const SIN_5: f64 = -0.0871557427476582;

#[allow(unused_mut)]
/// URA value (m) to URA index
fn uraindex(value: f64) -> i32 {
    let mut i = 0;
    while i < 15 {
        if URA_EPH[i] >= value {
            break;
        }
        i += 1;
    }
    return i as i32;
}

/// URA index to URA nominal value (m)
pub fn uravalue(sva: usize) -> f64 {
    const URA_NOMINAL: [f64; 16] = [
        2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0,
        2048.0, 4096.0, 8192.0,
    ];
    if sva < 15 {
        URA_NOMINAL[sva]
    } else {
        8192.0
    }
}

/// Galileo SISA index to SISA nominal value (m)
pub fn sisa_value(sisa:i32)->f64{
    let sisa=sisa.clone() as f64;
    if sisa<=49.0 {
        return sisa *0.01
    }
    if sisa<=74.0{
        return 0.5+(sisa- 50.0)*0.02
    }
    if sisa<=99.0{
        return 1.0+(sisa- 75.0)*0.04
    }
    if sisa<=125.0{
        return 2.0+(sisa-100.0)*0.16
    }
    -1.0
}

/// Galileo SISA value (m) to SISA index
fn sisa_index(value: f64) -> i32 {
    if value < 0.0 || value > 6.0 {
        return 255;
    } else if value <= 0.5 {
        return (value / 0.01) as i32;
    } else if value <= 1.0 {
        return ((value - 0.5) / 0.02) as i32 + 50;
    } else if value <= 2.0 {
        return ((value - 1.0) / 0.04) as i32 + 75;
    }
    return ((value - 2.0) / 0.16) as i32 + 100;
}

/// decode ephemeris
pub fn decode_eph(sat: usize, toc: GTime, data: &[f64], eph: &mut Eph) -> Result<bool, String> {
    *eph = Eph::default();

    let sys = satsys(sat);

    if sys & (SYS_GPS | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN) == 0 {
        // return Err(format!("ephemeris error: invalid satellite sat={}", sat));
        return Ok(false)
    }

    eph.sat = sat;
    eph.toc = toc;

    eph.f0 = data[0];
    eph.f1 = data[1];
    eph.f2 = data[2];

    eph.a = data[10].powi(2);
    eph.e = data[8];
    eph.i0 = data[15];
    eph.omg0 = data[13];
    eph.omg = data[17];
    eph.m0 = data[6];
    eph.deln = data[5];
    eph.omgd = data[18];
    eph.idot = data[19];
    eph.crc = data[16];
    eph.crs = data[4];
    eph.cuc = data[7];
    eph.cus = data[9];
    eph.cic = data[12];
    eph.cis = data[14];

    if sys == SYS_GPS || sys == SYS_QZS {
        eph.iode = data[3] as i32;
        eph.iodc = data[26] as i32;
        eph.toes = data[11];
        eph.week = data[21] as i32;
        eph.toe = adjweek(gpst2time(eph.week, data[11]), toc);
        eph.ttr = adjweek(gpst2time(eph.week, data[27]), toc);

        eph.code = data[20] as i32;
        eph.svh = data[24] as i32;
        eph.sva = uraindex(data[23]);
        eph.flag = data[22] as i32;

        eph.tgd[0] = data[25];
        eph.tgd[1] = data[26];

        if sys == SYS_GPS {
            eph.fit = data[28];
        } else {
            eph.fit = if data[28] == 0.0 { 1.0 } else { 2.0 };
        }
    } else if sys == SYS_GAL {
        eph.iode = data[3] as i32;
        eph.toes = data[11];
        eph.week = data[21] as i32;
        eph.toe = adjweek(gpst2time(eph.week, data[11]), toc);
        eph.ttr = adjweek(gpst2time(eph.week, data[27]), toc);

        eph.code = data[20] as i32;
        eph.svh = data[24] as i32;
        eph.sva = sisa_index(data[23]);

        eph.tgd[0] = data[25];
        eph.tgd[1] = data[26];
    } else if sys == SYS_CMP {
        eph.toc = bdt2gpst(eph.toc);
        eph.iode = data[3] as i32;
        eph.iodc = data[28] as i32;
        eph.toes = data[11];
        eph.week = data[21] as i32;
        eph.toe = bdt2gpst(bdt2time(eph.week, data[11]));
        eph.ttr = bdt2gpst(bdt2time(eph.week, data[27]));
        eph.toe = adjweek(eph.toe, toc);
        eph.ttr = adjweek(eph.ttr, toc);

        eph.svh = data[24] as i32;
        eph.sva = uraindex(data[23]);

        eph.tgd[0] = data[25];
        eph.tgd[1] = data[26];
    } else if sys == SYS_IRN {
        eph.iode = data[3] as i32;
        eph.toes = data[11];
        eph.week = data[21] as i32;
        eph.toe = adjweek(gpst2time(eph.week, data[11]), toc);
        eph.ttr = adjweek(gpst2time(eph.week, data[27]), toc);
        eph.svh = data[24] as i32;
        eph.sva = uraindex(data[23]);
        eph.tgd[0] = data[25];
    }

    Ok(true)
}

/// decode GLONASS ephemeris
pub fn decode_geph(
    ver: f64,
    sat: usize,
    mut toc: GTime,
    data: &[f64],
    geph: &mut Geph,
) -> Result<bool, String> {
    *geph = Geph::default();

    if satsys(sat) != SYS_GLO {
        return Err(format!(
            "glonass ephemeris error: invalid satellite sat={}",
            sat
        ));
    }

    geph.sat = sat;

    /* Toc rounded by 15 min in utc */
    let mut week = 0;
    let tow = time2gpst(toc, Some(&mut week));
    toc = gpst2time(week, ((tow + 450.0) / 900.0).floor() * 900.0);
    let dow = (tow / 86400.0).floor() as i32;

    /* time of frame in UTC */
    let tod = if ver <= 2.99 {
        data[2]
    } else {
        data[2] % 86400.0
    }; /* Tod (v.2), Tow (v.3) in UTC */
    let mut tof = gpst2time(week, tod + (dow * 86400) as f64);
    tof = adjday(tof, toc);

    geph.toe = utc2gpst(toc); /* Toc (GPST) */
    geph.tof = utc2gpst(tof); /* Tof (GPST) */

    /* IODE = Tb (7bit), Tb =index of UTC+3H within current day */
    geph.iode = ((tow + 10800.0) % 86400.0 / 900.0 + 0.5).floor() as i32;

    geph.taun = -data[0]; /* -taun */
    geph.gamn = data[1]; /* +gamman */

    geph.pos[0] = data[3] * 1E3;
    geph.pos[1] = data[7] * 1E3;
    geph.pos[2] = data[11] * 1E3;
    geph.vel[0] = data[4] * 1E3;
    geph.vel[1] = data[8] * 1E3;
    geph.vel[2] = data[12] * 1E3;
    geph.acc[0] = data[5] * 1E3;
    geph.acc[1] = data[9] * 1E3;
    geph.acc[2] = data[13] * 1E3;

    geph.svh = data[6] as i32;
    geph.frq = data[10] as i32;

    //geph.dtaun = data[14];
    geph.age = data[14] as i32;
    /* some receiver output >128 for minus frequency number */
    if geph.frq > 128 {
        geph.frq -= 256;
    }

    Ok(true)
}

/// broadcast ephemeris to satellite clock bias
fn eph2clk(time: GTime, eph: &Eph) -> f64 {
    let mut t = timediff(time, eph.toc);
    let ts = t;

    for _ in 0..2 {
        t = ts - (eph.f0 + eph.f1 * t + eph.f2 * t * t);
    }

    eph.f0 + eph.f1 * t + eph.f2 * t * t
}

/// broadcast ephemeris to satellite position and clock bias
fn geph2clk(time: GTime, geph: &Geph) -> f64 {
    let mut t = timediff(time, geph.toe);
    let ts = t;

    for _ in 0..2 {
        t = ts - (-geph.taun + geph.gamn * t);
    }

    -geph.taun + geph.gamn * t
}

/// select ephemeris
fn seleph(time: GTime, sat: usize, iode: i32, nav: &Nav) -> Option<&Eph> {
    //let _z=span!("seleph");
    let tmax;
    let mut tmin;
    let mut j = -1;
    let sys;

    sys = satsys(sat);
    match sys {
        SYS_GPS => {
            tmax = MAXDTOE + 1.0;
        }
        SYS_GAL => {
            tmax = MAXDTOE_GAL;
        }
        SYS_QZS => {
            tmax = MAXDTOE_QZS + 1.0;
        }
        SYS_CMP => {
            tmax = MAXDTOE_CMP + 1.0;
        }
        SYS_IRN => {
            tmax = MAXDTOE_IRN + 1.0;
        }
        _ => {
            tmax = MAXDTOE + 1.0;
        }
    }
    tmin = tmax + 1.0;

    for (i, eph) in nav.eph.iter().enumerate() {
        if eph.sat != sat {
            continue;
        }
        if iode >= 0 && eph.iode != iode {
            continue;
        }
        if sys == SYS_GAL {
            //let sel = 0;
            //let sel = getseleph(SYS_GAL);
            // if sel == 0 && (eph.code & (1 << 9)) == 0 {
            //     continue;
            // } // I/NAV
            // if sel == 1 && (eph.code & (1 << 8)) == 0 {
            //     continue;
            // } // F/NAV
            if timediff(eph.toe, time) >= 0.0 {
                continue;
            } // AOD <= 0
        }
        let t = timediff(eph.toe, time).abs();
        if t > tmax {
            continue;
        }
        if iode >= 0 {
            return Some(eph);
        }
        if t <= tmin {
            j = i as i32;
            tmin = t;
        }
    }
    if iode >= 0 || j < 0 {
        return None;
    }
    Some(&nav.eph[j as usize])
}

/// select glonass ephemeris
fn selgeph(time: GTime, sat: usize, iode: i32, nav: &Nav) -> Option<&Geph> {
    let tmax = MAXDTOE_GLO;
    let mut tmin = tmax + 1.0;
    let mut j = -1;

    for (i, geph) in nav.geph.iter().enumerate() {
        if geph.sat != sat {
            continue;
        }
        if iode >= 0 && geph.iode != iode {
            continue;
        }
        let t = timediff(geph.toe, time).abs();
        if t > tmax {
            continue;
        }
        if iode >= 0 {
            return Some(geph);
        }
        if t <= tmin {
            j = i as i32;
            tmin = t;
        }
    }

    if iode >= 0 || j < 0 {
        return None;
    }
    Some(&nav.geph[j as usize])
}

/// satellite clock with broadcast ephemeris
pub fn ephclk(time: GTime, teph: GTime, sat: usize, nav: &Nav, dts: &mut f64) -> bool {
    let sys = satsys(sat);
    match sys {
        SYS_GPS | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN => {
            if let Some(eph) = seleph(teph, sat, -1, nav) {
                *dts = eph2clk(time, eph);
            } else {
                return false;
            }
        }
        SYS_GLO => {
            if let Some(geph) = selgeph(teph, sat, -1, nav) {
                *dts = geph2clk(time, geph);
            } else {
                return false;
            }
        }
        _ => return false,
    }

    true
}

/// variance by ura ephemeris
fn var_uraeph(sys: usize, ura: i32) -> f64 {
    if sys == SYS_GAL {
        // Galileo SISA
        if ura <= 49 {
            sqr(ura as f64 * 0.01)
        } else if ura <= 74 {
            sqr(0.5 + (ura as f64 - 50.0) * 0.02)
        } else if ura <= 99 {
            sqr(1.0 + (ura as f64 - 75.0) * 0.04)
        } else if ura <= 125 {
            sqr(2.0 + (ura as f64 - 100.0) * 0.16)
        } else {
            sqr(STD_GAL_NAPA)
        }
    } else {
        // GPS URA
        if ura < 0 || ura > 14 {
            sqr(6144.0)
        } else {
            sqr(URA_EPH[ura as usize])
        }
    }
}

// broadcast ephemeris to satellite position and clock bias
#[allow(non_snake_case)]
fn eph2pos(time: GTime, eph: &Eph, rs: &mut [f64], dts: &mut f64, var: &mut f64) {
    let (tk, M, mut E, sinE, cosE, mut u, mut r, mut i, O, x, y, sinO, cosO, cosi, sin2u, cos2u, mu, omge): (
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    );
    let (xg, yg, zg, sino, coso);
    let mut iter: usize = 0;
    let sys;

    if eph.a <= 0.0 {
        rs[0] = 0.0;
        rs[1] = 0.0;
        rs[2] = 0.0;
        *dts = 0.0;
        *var = 0.0;
        return;
    }
    tk = timediff(time, eph.toe);

    let prn = getprn(eph.sat);
    sys = satsys(eph.sat);
    match sys {
        SYS_GAL => {
            mu = MU_GAL;
            omge = OMGE_GAL;
        }
        SYS_CMP => {
            mu = MU_CMP;
            omge = OMGE_CMP;
        }
        _ => {
            mu = MU_GPS;
            omge = OMGE;
        }
    }

    M = eph.m0 + (tk * (eph.deln + (mu / (eph.a * eph.a * eph.a)).sqrt()));

    E = M;

    for n in 0..MAX_ITER_KEPLER {
        iter = n;
        let de = (M - (E - eph.e * E.sin())) / (1.0 - eph.e * E.cos());
        E += de;
        if de.abs() <= RTOL_KEPLER {
            break;
        }
    }

    if iter == MAX_ITER_KEPLER {
        return;
    }

    sinE = E.sin();
    cosE = E.cos();

    u = ((1.0 - eph.e * eph.e).sqrt() * sinE).atan2(cosE - eph.e) + eph.omg;
    r = eph.a * (1.0 - eph.e * cosE);
    i = eph.i0 + eph.idot * tk;
    sin2u = (2.0 * u).sin();
    cos2u = (2.0 * u).cos();
    u += eph.cus * sin2u + eph.cuc * cos2u;
    r += eph.crs * sin2u + eph.crc * cos2u;
    i += eph.cis * sin2u + eph.cic * cos2u;
    cosi = i.cos();
    x = r * u.cos();
    y = r * u.sin();

    if sys == SYS_CMP && (prn <= 5 || prn >= 59) {
        O = eph.omg0 + eph.omgd * tk - omge * eph.toes;
        sinO = O.sin();
        cosO = O.cos();
        xg = x * cosO - y * cosi * sinO;
        yg = x * sinO + y * cosi * cosO;
        zg = y * i.sin();
        sino = (omge * tk).sin();
        coso = (omge * tk).cos();
        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
        rs[2] = -yg * SIN_5 + zg * COS_5;
    } else {
        O = eph.omg0 + (eph.omgd - omge) * tk - omge * eph.toes;
        sinO = O.sin();
        cosO = O.cos();
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * i.sin();
    }

    *dts = eph.f0 + eph.f1 * tk + eph.f2 * tk * tk
        - 2.0 * (mu * eph.a).sqrt() * eph.e * sinE / sqr(CLIGHT);
    *var = var_uraeph(sys, eph.sva);
}

// glonass orbit differential equations
fn deq(x: &[f64], xdot: &mut [f64], acc: &[f64]) {
    let x3 = DVector::from_row_slice(&x[0..3]);
    let r2 = x3.dot(&x3);
    let r3 = r2 * r2.sqrt();
    let omg2 = sqr(OMGE_GLO);

    if r2 <= 0.0 {
        for i in 0..6 {
            xdot[i] = 0.0;
        }
        return;
    }

    let a = 1.5 * J2_GLO * MU_GLO * sqr(RE_GLO) / r2 / r3; // 3/2*J2*mu*Ae^2/r^5
    let b = 5.0 * x[2] * x[2] / r2; // 5*z^2/r^2
    let c = -MU_GLO / r3 - a * (1.0 - b); // -mu/r^3-a(1-b)

    xdot[0] = x[3];
    xdot[1] = x[4];
    xdot[2] = x[5];
    xdot[3] = (c + omg2) * x[0] + 2.0 * OMGE_GLO * x[4] + acc[0];
    xdot[4] = (c + omg2) * x[1] - 2.0 * OMGE_GLO * x[3] + acc[1];
    xdot[5] = (c - 2.0 * a) * x[2] + acc[2];
}

// glonass position and velocity by numerical integration
fn glorbit(t: f64, x: &mut [f64], acc: &[f64]) {
    let mut k1 = vec![0.0; 6];
    let mut k2 = vec![0.0; 6];
    let mut k3 = vec![0.0; 6];
    let mut k4 = vec![0.0; 6];
    let mut w = vec![0.0; 6];

    deq(x, &mut k1, acc);
    for i in 0..6 {
        w[i] = x[i] + k1[i] * t / 2.0;
    }
    deq(&w, &mut k2, acc);
    for i in 0..6 {
        w[i] = x[i] + k2[i] * t / 2.0;
    }
    deq(&w, &mut k3, acc);
    for i in 0..6 {
        w[i] = x[i] + k3[i] * t;
    }
    deq(&w, &mut k4, acc);
    for i in 0..6 {
        x[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * t / 6.0;
    }
}

// glonass ephemeris to satellite position and clock bias
fn geph2pos(time: GTime, geph: &Geph, rs: &mut [f64], dts: &mut f64, var: &mut f64) {
    let mut t = timediff(time, geph.toe);
    let mut tt;
    let mut x = [0.0; 6];

    *dts = -geph.taun + geph.gamn * t;

    for i in 0..3 {
        x[i] = geph.pos[i];
        x[i + 3] = geph.vel[i];
    }

    tt = if t < 0.0 { -TSTEP } else { TSTEP };

    while t.abs() > 1E-9 {
        if t.abs() < TSTEP {
            tt = t;
        }
        glorbit(tt, &mut x, &geph.acc);
        t -= tt;
    }

    for i in 0..3 {
        rs[i] = x[i];
    }

    *var = sqr(ERREPH_GLO);
}

/// satellite position and clock by broadcast ephemeris
fn ephpos(
    time: GTime,
    teph: GTime,
    sat: usize,
    nav: &Nav,
    iode: i32,
    rs: &mut [f64],
    dts: &mut [f64],
    var: &mut f64,
    svh: &mut i32,
) -> bool {
    let mut rst = [0.0; 3];
    let mut dtst = [0.0; 1];
    let tt = 1E-3;
    let sys;
    sys = satsys(sat);

    *svh = -1;

    if sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP || sys == SYS_IRN {
        if let Some(eph) = seleph(teph, sat, iode, nav) {
            eph2pos(time, eph, rs, &mut dts[0], var);
            let time_tt = timeadd(time, tt);
            eph2pos(time_tt, eph, &mut rst, &mut dtst[0], var);
            *svh = eph.svh;
        } else {
            return false;
        }
    } else if sys == SYS_GLO {
        if let Some(geph) = selgeph(teph, sat, iode, nav) {
            geph2pos(time, geph, rs, &mut dts[0], var);
            let time_tt = timeadd(time, tt);
            geph2pos(time_tt, geph, &mut rst, &mut dtst[0], var);
            *svh = geph.svh;
        } else {
            return false;
        }
    } else {
        return false;
    }

    // Satellite velocity and clock drift by differential approx
    for i in 0..3 {
        rs[i + 3] = (rst[i] - rs[i]) / tt;
    }
    dts[1] = (dtst[0] - dts[0]) / tt;

    true
}


pub fn satpos(
    time: GTime,
    teph: GTime,
    sat: usize,
    nav: &mut Nav,
    rs: &mut [f64],
    dts: &mut [f64],
    var: &mut f64,
    svh: &mut i32,
) -> bool {
    return ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
}

/// satellite positions and clocks
pub fn satposs(
    teph: GTime,
    obs: &Vec<Obs>,
    n: usize,
    nav: &mut Nav,
    rs: &mut DMatrix<f64>,
    dts: &mut DMatrix<f64>,
    var: &mut DVector<f64>,
    svh: &mut [i32],
) {
    let mut time ;
    let mut dt: f64 = 0.0;
    let mut pr: f64;

    for i in 0..n.min(2 * MAXOBS) {
        for j in 0..6 {
            rs[(j, i)] = 0.0;
        }
        for j in 0..2 {
            dts[(j, i)] = 0.0;
        }

        var[i] = 0.0;
        svh[i] = 0;
        // search any pseudorange
        pr = 0.0;
        for j in 0..NFREQ {
            pr = obs[i].p[j];
            if pr != 0.0 {
                break;
            }
        }

        if pr == 0.0 {
            continue;
        }

        // transmission time by satellite clock
        time = timeadd(obs[i].time, -pr / CLIGHT);

        // satellite clock bias by broadcast ephemeris
        if !ephclk(time, teph, obs[i].sat, nav, &mut dt) {
            // trace(3, &format!("no broadcast clock {} sat={}", time_str(time[i], 3), obs[i].sat), 0);
            continue;
        }
        time = timeadd(time, -dt);
        // satellite position and clock at transmission time
        if !satpos(
            time,
            teph,
            obs[i].sat,
            nav,
            &mut rs.column_mut(i).as_mut_slice(),
            &mut dts.column_mut(i).as_mut_slice(),
            &mut var[i],
            &mut svh[i],
        ) {
            continue;
        }

        // if no precise clock available, use broadcast clock instead
        if dts[(0, i)] == 0.0 {
            if !ephclk(time, teph, obs[i].sat, nav, &mut dts[(0, i)]) {
                continue;
            }
            dts[(1, i)] = 0.0;
            var[i] = sqr(STD_BRDCCLK);
        }
    }
}
