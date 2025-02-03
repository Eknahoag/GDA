use crate::basic::code::sat2freq;
use crate::basic::eph::satposs;
use crate::basic::func::sqr;
use crate::basic::sat::*;
use crate::basic::time::{time2gpst, timeadd, timediff};
use crate::basic::var::*;
use nalgebra::{DMatrix, DVector, Matrix3, RealField, Vector2, Vector3, Vector4};
use std::ops::Sub;
use num_traits::Float;

/// rad to deg
const R2D: f64 = 180.0 / PI;
/// deg to rad
const D2R: f64 = PI / 180.0;
/// earth semimajor axis (WGS84) (m)
const RE_WGS84: f64 = 6378137.0;
/// max number of iteration for point pos
const MAXITR: usize = 10;
/// broadcast ionosphere model error factor
const ERR_BRDCI: f64 = 0.5;
#[allow(dead_code)]
/// ionospheric delay Std (m)
const ERR_ION: f64 = 5.0;
#[allow(dead_code)]
/// tropspheric delay Std (m)
const ERR_TROP: f64 = 3.0;
/// Saastamoinen model error Std (m)
const ERR_SAAS: f64 = 0.3;
/// code bias error Std (m)
const ERR_CBIAS: f64 = 0.3;
/// relative humidity for Saastamoinen model
const REL_HUMI: f64 = 0.7;
/// estimated parameters
const NX: usize = 8;
/// earth flattening (WGS84)
const FE_WGS84: f64 = 1.0 / 298.257223563;
/// error factor: GPS
const EFACT_GPS: f64 = 1.0;
/// error factor: GLONASS
const EFACT_GLO: f64 = 1.5;
/// min elevation for measurement error (rad)
const MIN_EL: f64 = 5.0 * D2R;
// chi-sqr(n) (alpha=0.001)
const CHISQR: [f64; 100] = [
    10.8, 13.8, 16.3, 18.5, 20.5, 22.5, 24.3, 26.1, 27.9, 29.6, 31.3, 32.9, 34.5, 36.1, 37.7, 39.3,
    40.8, 42.3, 43.8, 45.3, 46.8, 48.3, 49.7, 51.2, 52.6, 54.1, 55.5, 56.9, 58.3, 59.7, 61.1, 62.5,
    63.9, 65.2, 66.6, 68.0, 69.3, 70.7, 72.1, 73.4, 74.7, 76.0, 77.3, 78.6, 80.0, 81.3, 82.6, 84.0,
    85.4, 86.7, 88.0, 89.3, 90.6, 91.9, 93.3, 94.7, 96.0, 97.4, 98.7, 100.0, 101.0, 102.0, 103.0,
    104.0, 105.0, 107.0, 108.0, 109.0, 110.0, 112.0, 113.0, 114.0, 115.0, 116.0, 118.0, 119.0,
    120.0, 122.0, 123.0, 125.0, 126.0, 127.0, 128.0, 129.0, 131.0, 132.0, 133.0, 134.0, 135.0,
    137.0, 138.0, 139.0, 140.0, 142.0, 143.0, 144.0, 145.0, 147.0, 148.0, 149.0,
];
static mut IOB: usize = 0;

/// test excluded satellite
pub fn satexclude(sat: usize, var: f64, svh: &mut i32, opt: Option<&PrcOpt>) -> bool {
    let sys = satsys(sat);

    if *svh < 0 {
        return true; // ephemeris unavailable
    }

    if let Some(opt) = opt {
        if opt.exsats[sat - 1] == 1 {
            return true; // excluded satellite
        }
        if opt.exsats[sat - 1] == 2 {
            return false; // included satellite
        }
        if (sys & opt.navsys) == 0 {
            return true; // unselected sat sys
        }
    }

    if sys == SYS_QZS {
        *svh &= 0xFE; // mask QZSS LEX health
    }

    if *svh != 0 {
        return true;
    }

    if var > MAX_VAR_EPH {
        return true;
    }

    return false;
}

/// compute geometric distance and receiver-to-satellite unit vector
pub fn geodist(rs: &Vector3<f64>, rr: &[f64], e: &mut [f64; 3]) -> f64 {
    if rs.clone().norm() < RE_WGS84 {
        return -1.0;
    }

    for i in 0..3 {
        e[i] = rs[i] - rr[i];
    }

    let r = Vector3::from(e.clone()).norm();

    for i in 0..3 {
        e[i] /= r;
    }

    r + OMGE * (rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT
}

/// compute azimuth and elevation angles
pub fn satazel(pos: &[f64], e: &[f64], azel: Option<&mut [f64]>) -> f64 {
    let mut az = 0.0;
    let mut el = std::f64::consts::PI / 2.0;
    let enu;

    if pos[2] > -RE_WGS84 {
        enu = ecef2enu(&pos[0..2], e);
        let enu_xy = Vector2::new(enu[0], enu[1]);
        az = if enu_xy.dot(&enu_xy) < 1E-12 {
            0.0
        } else {
            enu[0].atan2(enu[1])
        };
        if az < 0.0 {
            az += 2.0 * PI;
        }
        el = enu[2].asin();
    }

    if let Some(azel) = azel {
        azel[0] = az;
        azel[1] = el;
    }

    return el;
}

/// compute ecef to local coordinate transfromation matrix
pub fn xyz2enu<T: RealField + Copy>(pos: &[T]) -> Matrix3<T> {
    let sinp = pos[0].sin();
    let cosp = pos[0].cos();
    let sinl = pos[1].sin();
    let cosl = pos[1].cos();

    Matrix3::new(
        -sinl,
        cosl,
        T::zero(),
        -sinp * cosl,
        -sinp * sinl,
        cosp,
        cosp * cosl,
        cosp * sinl,
        sinp,
    )
}

/// transform ecef vector to local tangental coordinate
pub fn ecef2enu<T: RealField + Copy>(pos: &[T], r: &[T]) -> Vector3<T> {
    let e = xyz2enu(pos);
    let r_vec = Vector3::new(r[0], r[1], r[2]);
    e * r_vec
}

/// transform local tangental coordinate vector to ecef
pub fn ecef2pos<T: RealField + Float + Sub<Output=T>>(r: &[T], pos: &mut [T]) {
    let e2 = T::from(FE_WGS84 * (2.0 - FE_WGS84)).unwrap();
    let r2 = Vector2::new(r[0], r[1]).dot(&Vector2::new(r[0], r[1]));
    let mut z = r[2];
    let mut zk = T::zero();
    let mut v = T::from(RE_WGS84).unwrap();
    let mut sinp: T;

    while Float::abs(z - zk) >= T::from(1E-4).unwrap() {
        zk = z;
        sinp = z / Float::sqrt(r2 + z * z);
        v = T::from(RE_WGS84).unwrap() / Float::sqrt(T::from(1.0).unwrap() - e2 * sinp * sinp);
        z = r[2] + v * e2 * sinp;
    }
    pos[0] = if r2 > T::from(1E-12).unwrap() {
        Float::atan(z / Float::sqrt(r2))
    } else {
        if r[2] > T::zero() {
            T::from(PI / 2.0).unwrap()
        } else {
            T::from(-PI / 2.0).unwrap()
        }
    };
    pos[1] = if r2 > T::from(1E-12).unwrap() {
        Float::atan2(r[1], r[0])
    } else {
        T::zero()
    };
    pos[2] = Float::sqrt(r2 + z * z) - v;
}

/// test snr mask
pub fn snrmask(obs: &Obs, azel: &Vector2<f64>, opt: PrcOpt) -> bool {
    if testsnr(0, 0, azel[1], obs.snr[0] as f64 * SNR_UNIT, &opt.snrmask) {
        return false;
    }
    if opt.ionoopt == IONOOPT_IFLC {
        if testsnr(0, 1, azel[1], obs.snr[0] as f64 * SNR_UNIT, &opt.snrmask) {
            return false;
        }
    }
    true
}

// test snr mask
fn testsnr(base: usize, idx: usize, el: f64, snr: f64, mask: &SnrMask) -> bool {
    if mask.ena[base] == 0 || idx >= NFREQ {
        return false;
    }

    let mut a = (el * R2D + 5.0) / 10.0;
    let i = a.floor() as usize;
    a -= i as f64;

    let minsnr = if i < 1 {
        mask.mask[idx][0]
    } else if i > 8 {
        mask.mask[idx][8]
    } else {
        (1.0 - a) * mask.mask[idx][i - 1] + a * mask.mask[idx][i]
    };

    return snr < minsnr;
}

// compute ionospheric delay by broadcast ionosphere model (klobuchar model)
fn ionmodel(t: GTime, ion: &mut [f64; 8], pos: &[f64], azel: &Vector2<f64>) -> f64 {
    let ion_default = [
        0.1118E-07,
        -0.7451E-08,
        -0.5961E-07,
        0.1192E-06,
        0.1167E+06,
        -0.2294E+06,
        -0.1311E+06,
        0.1049E+07,
    ];

    let psi;
    let mut phi;
    let lam;
    let amp;
    let per;
    let x;
    let f;
    let mut week = 0;

    if pos[2] < -1E3 || azel[1] <= 0.0 {
        return 0.0;
    }

    if ion.iter().all(|&v| v == 0.0) {
        *ion = ion_default;
    }

    // Earth centered angle (semicircle)
    psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;

    // Subionospheric latitude/longitude (semi-circle)
    phi = pos[0] / PI + psi * azel[0].cos();
    if phi > 0.416 {
        phi = 0.416;
    } else if phi < -0.416 {
        phi = -0.416;
    }
    lam = pos[1] / PI + psi * azel[0].sin() / (phi * PI).cos();

    // Geomagnetic latitude (semi-circle)
    phi += 0.064 * ((lam - 1.617) * PI).cos();

    // Local time (s)
    let mut tt = 43200.0 * lam + time2gpst(t, Some(&mut week));
    tt = tt - (tt / 86400.0).floor() * 86400.0; // 0 <= tt < 86400

    // Slant factor
    f = 1.0 + 16.0 * (0.53 - azel[1] / PI).powi(3);

    // Ionospheric delay
    amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]));
    per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]));
    let amp = if amp < 0.0 { 0.0 } else { amp };
    let per = if per < 72000.0 { 72000.0 } else { per };
    x = 2.0 * PI * (tt - 50400.0) / per;

    CLIGHT
        * f
        * if x.abs() < 1.57 {
        5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0))
    } else {
        5E-9
    }
}

// compute ionospheric correction
fn ionocorr(
    time: GTime,
    nav: &mut Nav,
    pos: &[f64],
    azel: &Vector2<f64>,
    ion: &mut f64,
    var: &mut f64,
) -> bool {
    *ion = ionmodel(time, &mut nav.ion_gps, pos, azel);
    *var = sqr(*ion * ERR_BRDCI);
    true
}

// compute tropospheric delay by standard atmosphere and saastamoinen model
fn tropmodel(pos: &[f64], azel: &Vector2<f64>, humi: f64) -> f64 {
    const TEMP0: f64 = 15.0; // temperature at sea level
    let hgt: f64;
    let pres: f64;
    let temp: f64;
    let e: f64;
    let z: f64;
    let trph: f64;
    let trpw: f64;

    if pos[2] < -100.0 || pos[2] > 1E4 || azel[1] <= 0.0 {
        return 0.0;
    }

    // standard atmosphere
    hgt = if pos[2] < 0.0 { 0.0 } else { pos[2] };
    pres = 1013.25 * (1.0 - 2.2557E-5 * hgt).powf(5.2568);
    temp = TEMP0 - 6.5E-3 * hgt + 273.16;
    e = 6.108 * humi * ((17.15 * temp - 4684.0) / (temp - 38.45)).exp();

    // Saastamoinen model
    z = PI / 2.0 - azel[1];
    trph =
        0.0022768 * pres / (1.0 - 0.00266 * (2.0 * pos[0]).cos() - 0.00028 * hgt / 1E3) / z.cos();
    trpw = 0.002277 * (1255.0 / temp + 0.05) * e / z.cos();

    trph + trpw
}

// compute tropospheric correction
fn tropcorr(pos: &[f64; 3], azel: &Vector2<f64>, trp: &mut f64, var: &mut f64) -> bool {
    *trp = tropmodel(pos, azel, REL_HUMI);
    *var = sqr(ERR_SAAS / (azel[1].sin() + 0.1));
    true
}

// get group delay parameter (m)
fn gettgd(sat: usize, nav: &Nav, r#type: usize) -> f64 {
    let sys = satsys(sat);

    if sys == SYS_GLO {
        for i in 0..nav.ng {
            if nav.geph[i].sat == sat {
                return -nav.geph[i].dtaun * CLIGHT;
            }
        }
        0.0
    } else {
        for i in 0..nav.n {
            if nav.eph[i].sat == sat {
                return nav.eph[i].tgd[r#type] * CLIGHT;
            }
        }
        0.0
    }
}

/// psendorange with code bias correction
fn prange(obs: &Obs, nav: &Nav, opt: &PrcOpt, var: &mut f64) -> f64 {
    let mut p1 = obs.p[0];
    let mut p2 = obs.p[1];
    let b1;
    let b2;
    let gamma;
    let sat = obs.sat;
    let sys = satsys(sat);
    *var = 0.0;

    if p1 == 0.0 || (opt.ionoopt == IONOOPT_IFLC && p2 == 0.0) {
        return 0.0;
    }

    if sys == SYS_GPS || sys == SYS_GLO {
        if obs.code[0] == CODE_L1C {
            p1 += nav.cbias[sat - 1][1];
        }
        if obs.code[1] == CODE_L2C {
            p2 += nav.cbias[sat - 1][2];
        }
    }

    if opt.ionoopt == IONOOPT_IFLC {
        // dual-frequency
        if sys == SYS_GPS || sys == SYS_QZS {
            gamma = sqr(FREQ1 / FREQ2);
            return (p2 - gamma * p1) / (1.0 - gamma);
        } else if sys == SYS_GLO {
            gamma = sqr(FREQ1_GLO / FREQ2_GLO);
            return (p2 - gamma * p1) / (1.0 - gamma);
        } else if sys == SYS_GAL {
            gamma = sqr(FREQ1 / FREQ7);
            // if getseleph(SYS_GAL) {
            //     p2 -= gettgd(sat, nav, 0) - gettgd(sat, nav, 1);
            // }
            return (p2 - gamma * p1) / (1.0 - gamma);
        } else if sys == SYS_CMP {
            gamma = sqr(if obs.code[0] == CODE_L2I {
                FREQ1_CMP
            } else {
                FREQ1
            } / FREQ2_CMP);
            if obs.code[0] == CODE_L2I {
                b1 = gettgd(sat, nav, 0);
            } else if obs.code[0] == CODE_L1P {
                b1 = gettgd(sat, nav, 2);
            } else {
                b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4);
            }
            b2 = gettgd(sat, nav, 1);
            return ((p2 - gamma * p1) - (b2 - gamma * b1)) / (1.0 - gamma);
        } else if sys == SYS_IRN {
            gamma = sqr(FREQ5 / FREQ9);
            return (p2 - gamma * p1) / (1.0 - gamma);
        }
    } else {
        // single-freq (L1/E1/B1)
        *var = sqr(ERR_CBIAS);

        if sys == SYS_GPS || sys == SYS_QZS {
            b1 = gettgd(sat, nav, 0);
            return p1 - b1;
        } else if sys == SYS_GLO {
            gamma = sqr(FREQ1_GLO / FREQ2_GLO);
            b1 = gettgd(sat, nav, 0);
            return p1 - b1 / (gamma - 1.0);
        } else if sys == SYS_GAL {
            b1 = gettgd(sat, nav, 1);
            //if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0);
            return p1 - b1;
        } else if sys == SYS_CMP {
            b1 = if obs.code[0] == CODE_L2I {
                gettgd(sat, nav, 0)
            } else if obs.code[0] == CODE_L1P {
                gettgd(sat, nav, 2)
            } else {
                gettgd(sat, nav, 2) + gettgd(sat, nav, 4)
            };
            return p1 - b1;
        } else if sys == SYS_IRN {
            gamma = sqr(FREQ9 / FREQ5);
            b1 = gettgd(sat, nav, 0);
            return p1 - gamma * b1;
        }
    }
    p1
}

/// pseudorange residuals
fn rescode(
    iter: usize,
    obs: &Vec<Obs>,
    n: usize,
    rs: &DMatrix<f64>,
    dts: &DMatrix<f64>,
    vare: &DVector<f64>,
    svh: &mut [i32],
    nav: &mut Nav,
    x: &DVector<f64>,
    opt: &PrcOpt,
    v: &mut DVector<f64>,
    h: &mut DMatrix<f64>,
    var: &mut [f64],
    azel: &mut DMatrix<f64>,
    vsat: &mut [i32],
    resp: &mut DVector<f64>,
    ns: &mut i32,
) -> usize {
    let mut time;
    let mut r;
    let mut freq;
    let mut rr = [0.0; 3];
    let mut pos = [0.0; 3];
    let dtr;
    let mut e = [0.0; 3];
    let mut p;
    let mut dion = 0.0;
    let mut dtrp = 0.0;
    let mut vmeas = 0.0;
    let mut vion = 0.0;
    let mut vtrp = 0.0;
    let mut nv = 0;
    let mut sat;
    let mut sys;
    let mut mask = [0; NX - 3];

    for i in 0..3 {
        rr[i] = x[i];
    }
    dtr = x[3];

    ecef2pos(&rr, &mut pos);

    for i in 0..n.min(MAXOBS) {
        vsat[i] = 0;
        azel[(0, i)] = 0.0;
        azel[(1, i)] = 0.0;
        resp[i] = 0.0;
        time = obs[i].time;
        sat = obs[i].sat;
        sys = satsys(sat);
        if sys == 0 {
            continue;
        }

        if i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat {
            continue;
        }

        if satexclude(sat, vare[i], &mut svh[i], Some(opt)) {
            continue;
        }
        r = geodist(
            &Vector3::new(rs[(0, i)], rs[(1, i)], rs[(2, i)]),
            &rr,
            &mut e,
        );
        if r <= 0.0 {
            continue;
        }

        if iter > 0 {
            if satazel(
                &pos,
                &e,
                Some(&mut azel.column_mut(i).as_mut_slice()),
            ) < opt.elmin
            {
                continue;
            }

            if !snrmask(&obs[i], &Vector2::new(azel[(0, i)], azel[(1, i)]), *opt) {
                continue;
            }

            if !ionocorr(
                time,
                nav,
                &pos,
                &Vector2::new(azel[(0, i)], azel[(1, i)]),
                &mut dion,
                &mut vion,
            ) {
                continue;
            }
            freq = sat2freq(sat, obs[i].code[0], Some(nav));
            if freq == 0.0 {
                continue;
            }
            dion *= (FREQ1 / freq).powi(2);
            vion *= (FREQ1 / freq).powi(2);

            if !tropcorr(
                &pos,
                &mut Vector2::new(azel[(0, i)], azel[(1, i)]),
                &mut dtrp,
                &mut vtrp,
            ) {
                continue;
            }
        }
        p = prange(&obs[i], nav, opt, &mut vmeas);
        if p == 0.0 {
            continue;
        }

        v[nv] = p - (r + dtr - CLIGHT * dts[(0, i)] + dion + dtrp);

        for j in 0..NX {
            h[(j, nv)] = if j < 3 {
                -e[j]
            } else {
                if j == 3 {
                    1.0
                } else {
                    0.0
                }
            };
        }

        match sys {
            SYS_GLO => {
                v[nv] -= x[4];
                h[(4, nv)] = 1.0;
                mask[1] = 1;
            }
            SYS_GAL => {
                v[nv] -= x[5];
                h[(5, nv)] = 1.0;
                mask[2] = 1;
            }
            SYS_CMP => {
                v[nv] -= x[6];
                h[(6, nv)] = 1.0;
                mask[3] = 1;
            }
            SYS_IRN => {
                v[nv] -= x[7];
                h[(7, nv)] = 1.0;
                mask[4] = 1;
            }
            _ => mask[0] = 1,
        }

        vsat[i] = 1;
        resp[i] = v[nv];
        *ns += 1;

        var[nv] = varerr(opt, azel[(1, i)], sys) + vare[i] + vmeas + vion + vtrp;
        nv += 1;
    }

    for i in 0..NX - 3 {
        if mask[i] != 0 {
            continue;
        }
        v[nv] = 0.0;
        for j in 0..NX {
            h[(j, nv)] = if j == i + 3 { 1.0 } else { 0.0 };
        }
        var[nv] = 0.01;
        nv += 1;
    }

    nv
}

// pseudorange measurement error variance
fn varerr(opt: &PrcOpt, el: f64, sys: usize) -> f64 {
    let fact = match sys {
        SYS_GLO => EFACT_GLO,
        _ => EFACT_GPS,
    };

    let el = if el < MIN_EL { MIN_EL } else { el };
    let mut varr = opt.err[0].powi(2) * (opt.err[1].powi(2) + opt.err[2].powi(2) / el.sin());
    varr = if opt.ionoopt == IONOOPT_IFLC {
        varr * 3.0_f64.powi(2)
    } else {
        varr
    };

    fact.powi(2) * varr
}

// least square estimation by solving normal equation (x=(A*A')^-1*A*y
fn lsq(a: &DMatrix<f64>, y: &DVector<f64>, x: &mut DVector<f64>, q: &mut DMatrix<f64>) -> bool {
    let ay = a * y;
    *q = a * a.transpose();

    if let Some(inv_q) = q.clone().try_inverse() {
        // x = Q^-1 * Ay
        *x = inv_q * ay;
        true
    } else {
        false
    }
}

/// compute DOPs (dilution of precision)
pub fn dops(ns: usize, azel: &DMatrix<f64>, elmin: f64, dop: &mut [f64; 4]) {
    let mut h = DMatrix::zeros(MAXSAT, 4);
    let mut n = 0;

    for i in 0..4 {
        dop[i] = 0.0;
    }

    for i in 0..ns.min(MAXSAT) {
        if azel[(1, i)] < elmin || azel[(1, i)] <= 0.0 {
            continue;
        }

        let cosel = azel[(1, i)].cos();
        let sinel = azel[(1, i)].sin();

        h[(n, 0)] = cosel * azel[i * 2].sin();
        h[(n, 1)] = cosel * azel[i * 2].cos();
        h[(n, 2)] = sinel;
        h[(n, 3)] = 1.0;

        n += 1;
    }

    if n < 4 {
        return;
    }

    let h = h.rows(0,n).into_owned(); // 只保留有效行
    let q = h.transpose() * &h;

    if let Some(inv_q) = q.clone().try_inverse() {
        dop[0] = (inv_q[(0, 0)] + inv_q[(1, 1)] + inv_q[(2, 2)] + inv_q[(3, 3)]).sqrt(); // GDOP
        dop[1] = (inv_q[(0, 0)] + inv_q[(1, 1)] + inv_q[(2, 2)]).sqrt(); // PDOP
        dop[2] = (inv_q[(0, 0)] + inv_q[(1, 1)]).sqrt(); // HDOP
        dop[3] = inv_q[(2, 2)].sqrt(); // VDOP
    }
}

// validate solution
fn valsol(
    azel: &DMatrix<f64>,
    vsat: &[i32],
    n: usize,
    opt: &PrcOpt,
    v: &DVector<f64>,
    nv: usize,
    nx: usize,
) -> Result<bool, String> {
    let mut azels = DMatrix::zeros(2, MAXOBS);
    let mut dop = [0.0; 4];

    // Chi-square validation of residuals
    let vv = v.dot(v);
    if nv > nx && vv > CHISQR[nv - nx - 1] {
        return Err(String::from("chi-square error"));
    }

    // Large GDOP check
    let mut ns = 0;
    for i in 0..n {
        if vsat[i] == 0 {
            continue;
        }
        azels[(0, ns)] = azel[(0, i)].clone();
        azels[(1, ns)] = azel[(1, i)].clone();
        ns += 1;
    }

    dops(ns, &azels, opt.elmin, &mut dop);
    if dop[0] <= 0.0 || dop[0] > opt.maxgdop {
        return Err(String::from("gdop error"));
    }
    Ok(true)
}

/// estimate receiver position
pub fn estpos(
    obs: &Vec<Obs>,
    n: usize,
    rs: &DMatrix<f64>,
    dts: &DMatrix<f64>,
    vare: &DVector<f64>,
    svh: &mut [i32],
    nav: &mut Nav,
    opt: &PrcOpt,
    sol: &mut SolOpt,
    azel: &mut DMatrix<f64>,
    vsat: &mut [i32],
    resp: &mut DVector<f64>,
) -> Result<bool, String> {
    let mut x = DVector::zeros(NX);
    let mut dx = DVector::zeros(NX);
    let mut q = DMatrix::zeros(NX, NX);
    let mut v = DVector::zeros(n + 4);
    let mut h = DMatrix::zeros(NX, n + 4);
    let mut var = vec![0.0; n + 4];
    let mut iter = 0;

    for i in 0..3 {
        x[i] = sol.rr[i];
    }

    for i in 0..MAXITR {
        iter = i;
        let mut ns = 0;
        let nv = rescode(
            i, obs, n, rs, dts, vare, svh, nav, &x, opt, &mut v, &mut h, &mut var, azel, vsat,
            resp, &mut ns,
        );
        // for j in 0..nv{println!("{:20.8}",v[j])}
        if nv < NX {
            return Err(String::from("lack of valid sats ns=%d"));
        }

        // Weighted by Std
        for j in 0..nv {
            let sig = var[j].sqrt();
            v[j] /= sig;
            for k in 0..NX {
                h[(k, j)] /= sig;
            }
        }
        // for j in 0..(n+4)*NX{println!("{:20.8}",h[j])}

        // Least square estimation
        if !lsq(&h.columns(0, nv).into_owned(), &v.rows(0, nv).into_owned(), &mut dx, &mut q) {
            return Err(String::from("lsq error"));
        }

        for j in 0..NX {
            x[j] += dx[j];
        }

        if dx.norm() < 1E-4 {
            sol.type_ = 0;
            sol.time = timeadd(obs[0].time, -x[3] / CLIGHT);
            sol.dtr[0] = x[3] / CLIGHT;
            sol.dtr[1] = x[4] / CLIGHT;
            sol.dtr[2] = x[5] / CLIGHT;
            sol.dtr[3] = x[6] / CLIGHT;
            sol.dtr[4] = x[7] / CLIGHT;
            for j in 0..3 {
                sol.rr[j] = x[j];
            }
            for j in 0..3 {
                sol.qr[j] = q[(j, j)];
            }
            sol.qr[3] = q[(0, 1)]; // cov xy
            sol.qr[4] = q[(1, 2)]; // cov yz
            sol.qr[5] = q[(2, 0)]; // cov zx
            sol.ns = ns as u8;
            sol.age = 0.0;
            sol.ratio = 0.0;

            // Validate solution
            if let Ok(true) = valsol(azel, vsat, n, opt, &v.rows(0, nv).into_owned(), nv, NX) {
                sol.stat = if opt.sateph == 2 { 1 } else { 0 }; // SOLQ_SBAS: 1, SOLQ_SINGLE: 0
                return Ok(true);
            }
            return Ok(false);
        }
    }

    if iter >= MAXITR {
        return Err(String::from("over MAXITR"));
    }

    Err(String::from("fail iter in MAXITR"))
}

#[allow(non_snake_case)]
// range rate residuals
fn resdop(
    obs: &[Obs],
    n: usize,
    rs: &mut DMatrix<f64>,
    dts: &mut DMatrix<f64>,
    nav: &Nav,
    rr: &[f64],
    x: &Vector4<f64>,
    azel: &mut DMatrix<f64>,
    vsat: &[i32],
    err: f64,
    v: &mut DVector<f64>,
    H: &mut DMatrix<f64>,
) -> usize {
    let mut pos = [0.0; 3];
    let mut nv = 0;

    ecef2pos(rr, &mut pos);
    let E = xyz2enu(&pos);

    for i in 0..n.min(MAXOBS) {
        let freq = sat2freq(obs[i].sat, obs[i].code[0], Some(nav));

        if obs[i].d[0] == 0.0
            || freq == 0.0
            || vsat[i] == 0
            || Vector3::new(rs[(3, i)], rs[(4, i)], rs[(5, i)]).norm() <= 0.0
        {
            continue;
        }

        // LOS (line-of-sight) vector in ECEF
        let cosel = azel[(1, i)].cos();
        let a = Vector3::new(
            azel[(0, i)].sin() * cosel,
            azel[(0, i)].cos() * cosel,
            azel[(1, i)].sin(),
        );
        let e = E.transpose() * a;

        // Satellite velocity relative to receiver in ECEF
        let vs = Vector3::new(rs[(3, i)] - x[0], rs[(4, i)] - x[1], rs[(5, i)] - x[2]);

        // Range rate with earth rotation correction
        let rate = vs.dot(&e)
            + OMGE / CLIGHT
            * (rs[(4, i)] * rr[0] + rs[(1, i)] * x[0]
            - rs[(3, i)] * rr[1]
            - rs[(0, i)] * x[1]);

        // Std of range rate error (m/s)
        let sig = if err <= 0.0 { 1.0 } else { err * CLIGHT / freq };

        // Range rate residual (m/s)
        v[nv] = (-obs[i].d[0] * CLIGHT / freq - (rate + x[3] - CLIGHT * dts[(1, i)])) / sig;

        // Design matrix
        for j in 0..4 {
            H[j + nv * 4] = (if j < 3 { -e[j] } else { 1.0 }) / sig;
        }
        nv += 1;
    }
    nv
}

#[allow(non_snake_case)]
// estimate receiver velocity
pub fn estvel(
    obs: &[Obs],
    n: usize,
    rs: &mut DMatrix<f64>,
    dts: &mut DMatrix<f64>,
    nav: &Nav,
    opt: &PrcOpt,
    sol: &mut SolOpt,
    azel: &mut DMatrix<f64>,
    vsat: &[i32],
) {
    let mut x = Vector4::zeros();
    let mut dx = DVector::zeros(4);
    let mut Q = DMatrix::zeros(4, 4);
    let err = opt.err[4]; // Doppler error (Hz)
    let mut v = DVector::zeros(n);
    let mut H = DMatrix::zeros(4, n);

    for _ in 0..MAXITR {
        // Range rate residuals (m/s)
        let nv = resdop(
            obs,
            n,
            rs,
            dts,
            nav,
            &sol.rr[0..3],
            &x,
            azel,
            vsat,
            err,
            &mut v,
            &mut H,
        );
        if nv < 4 {
            break;
        }

        // Least square estimation
        if !lsq(&H, &v, &mut dx, &mut Q) {
            break;
        }

        x += dx.clone();

        if dx.clone().norm() < 1E-6 {
            sol.rr[3] = x[0];
            sol.rr[4] = x[1];
            sol.rr[5] = x[2];
            sol.qv[0] = Q[(0, 0)]; /* xx */
            sol.qv[1] = Q[(1, 1)]; /* yy */
            sol.qv[2] = Q[(2, 2)]; /* zz */
            sol.qv[3] = Q[(0, 1)]; /* xy */
            sol.qv[4] = Q[(1, 2)]; /* yz */
            sol.qv[5] = Q[(0, 2)]; /* zx */
            break;
        }
    }
}

/// compute receiver position, velocity, clock bias by single-point positioning
///
/// # Arguments
/// - `obs`: observation data vector
/// - `n`: number of observation data
/// - `nav`: navigation data
/// - `opt`: processing options
/// - `azel`: azimuth/elevation angle
/// - `ssat`: satellite status
///
/// # Returns
/// status
pub fn pntpos(
    obs: &Vec<Obs>,
    n: usize,
    nav: &mut Nav,
    opt: &PrcOpt,
    sol: &mut SolOpt,
    azel: Option<&mut DMatrix<f64>>,
    ssat: Option<&mut Vec<Ssat>>,
) -> Result<bool, String> {
    let opt_ = opt.clone();
    let mut rs = DMatrix::zeros(6, n);
    let mut dts = DMatrix::zeros(2, n);
    let mut var = DVector::zeros(n);
    let mut azel_ = DMatrix::zeros(2, n);
    let mut resp = DVector::zeros(n);
    let mut dop = [0.0; 4];
    let mut vsat = [0; MAXOBS];
    let mut svh = [0; MAXOBS];

    sol.stat = 0;

    if n <= 0 {
        return Err(String::from("no observation data"));
    }

    sol.time = obs[0].time;

    // if opt_.mode != PMODE_SINGLE {
    //     opt_.ionoopt = IONOOPT_BRDC;
    //     opt_.tropopt = TROPOPT_SAAS;
    // }

    satposs(sol.time, obs, n, nav, &mut rs, &mut dts, &mut var, &mut svh);
    /*for i in 0..n{
        println!("{:20.8} {:20.8} {:20.8}",rs[(0,i)],rs[(1,i)],rs[(2,i)]);
    }*/
    let stat = estpos(
        obs, n, &rs, &dts, &var, &mut svh, nav, &opt_, sol, &mut azel_, &mut vsat, &mut resp,
    );

    if stat == Ok(true) {
        estvel(
            obs, n, &mut rs, &mut dts, nav, &opt_, sol, &mut azel_, &vsat,
        );
    }

    if let Some(azel) = azel {
        azel.copy_from(&azel_);
    }

    dops(n, &azel_, opt.elmin, &mut dop);

    if let Some(ssat) = ssat {
        for s in ssat.iter_mut() {
            s.vs = 0;
            s.azel = [0.0, 0.0];
            s.resp[0] = 0.0;
            s.resc[0] = 0.0;
            s.snr[0] = 0;
            s.dop = [0.0, 0.0, 0.0, 0.0];
        }

        for i in 0..n {
            let sat_idx = obs[i].sat - 1;
            ssat[sat_idx].azel[0] = azel_[(0, i)];
            ssat[sat_idx].azel[1] = azel_[(1, i)];
            ssat[sat_idx].snr[0] = obs[i].snr[0];

            if vsat[i] == 0 {
                continue;
            }

            ssat[sat_idx].vs = 1;
            ssat[sat_idx].resp[0] = resp[i];
            ssat[sat_idx].dop.copy_from_slice(&dop);
        }
    }

    stat
}

// search next observation data index
fn nextobsf(obss: &Obss, i: &mut usize, rcv: u8) -> usize {
    let mut n = 0;

    // Find the first observation with the given receiver number
    while *i < obss.n {
        if obss.obs[*i].rcv == rcv {
            break;
        }
        *i += 1;
    }

    // Count the number of consecutive observations with the same receiver number
    while *i + n < obss.n {
        let tt = timediff(obss.obs[*i + n].time.clone(), obss.obs[*i].time.clone());
        if obss.obs[*i + n].rcv != rcv || tt > DTTOL {
            break;
        }
        n += 1;
    }

    n
}

/// Extract the obs of a single epoch from all obs and return number
pub fn inputobs(obss: &Obss, obs: &mut Vec<Obs>) -> i32 {
    unsafe {
        let mut temp = IOB;
        let nu = nextobsf(obss, &mut temp, 1);
        if nu <= 0 {
            return -1;
        }
        IOB = temp;
        for i in 0..nu {
            obs.push(obss.obs[IOB + i].clone());
        }
        IOB += nu;
        return nu as i32;
    }
}
