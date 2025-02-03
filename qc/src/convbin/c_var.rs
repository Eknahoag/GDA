#![allow(non_snake_case)]
use std::ffi::CStr;
use std::os::raw::{c_char, c_double, c_float, c_int, c_uchar, c_uint, c_ushort};
use std::slice;
use libc::c_short;
use crate::basic::var::*;

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GTime {
    pub time: i64,
    pub sec: c_double,
}

#[allow(non_snake_case)]
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Obsd {
    pub time: GTime,
    pub sat: c_uchar,
    pub rcv: c_uchar,
    pub SNR: [c_ushort; NFREQ + NEXOBS],
    pub LLI: [c_uchar; NFREQ + NEXOBS],
    pub code: [c_uchar; NFREQ + NEXOBS],
    pub L: [c_double; NFREQ + NEXOBS],
    pub P: [c_double; NFREQ + NEXOBS],
    pub D: [c_float; NFREQ + NEXOBS],
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Obs {
    pub n: c_int,
    pub nmax: c_int,
    pub data: *mut Obsd,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Eph {
    pub sat: c_int,
    pub iode: c_int,
    pub iodc: c_int,
    pub sva: c_int,
    pub svh: c_int,
    pub week: c_int,
    pub code: c_int,
    pub flag: c_int,
    pub toe: GTime,
    pub toc: GTime,
    pub ttr: GTime,
    pub A: c_double,
    pub e: c_double,
    pub i0: c_double,
    pub OMG0: c_double,
    pub omg: c_double,
    pub M0: c_double,
    pub deln: c_double,
    pub OMGd: c_double,
    pub idot: c_double,
    pub crc: c_double,
    pub crs: c_double,
    pub cuc: c_double,
    pub cus: c_double,
    pub cic: c_double,
    pub cis: c_double,
    pub toes: c_double,
    pub fit: c_double,
    pub f0: c_double,
    pub f1: c_double,
    pub f2: c_double,
    pub tgd: [c_double; 6],
    pub Adot: c_double,
    pub ndot: c_double,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Geph {
    pub sat: c_int,
    pub iode: c_int,
    pub frq: c_int,
    pub svh: c_int,
    pub sva: c_int,
    pub age: c_int,
    pub toe: GTime,
    pub tof: GTime,
    pub pos: [c_double; 3],
    pub vel: [c_double; 3],
    pub acc: [c_double; 3],
    pub taun: c_double,
    pub gamn: c_double,
    pub dtaun: c_double,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Seph {
    pub sat: c_int,              // 卫星编号
    pub t0: GTime,               // 参考历元时间（GPST）
    pub tof: GTime,              // 消息帧时间（GPST）
    pub sva: c_int,              // SV 精度（URA 索引）
    pub svh: c_int,              // SV 健康状态（0: 正常）
    pub pos: [c_double; 3],      // 卫星位置 (m) (ECEF)
    pub vel: [c_double; 3],      // 卫星速度 (m/s) (ECEF)
    pub acc: [c_double; 3],      // 卫星加速度 (m/s^2) (ECEF)
    pub af0: c_double,           // 卫星时钟偏移 (s)
    pub af1: c_double,           // 卫星时钟漂移 (s/s)
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Peph {
    pub time: GTime,                     // 时间（GPST）
    pub index: c_int,                    // 多文件的星历索引
    pub pos: [[c_double; 4]; MAXSAT],    // 卫星位置/时钟（ECEF）(m|s)
    pub std: [[c_float; 4]; MAXSAT],     // 卫星位置/时钟标准差 (m|s)
    pub vel: [[c_double; 4]; MAXSAT],    // 卫星速度/时钟频率 (m/s|s/s)
    pub vst: [[c_float; 4]; MAXSAT],     // 卫星速度/时钟频率标准差 (m/s|s/s)
    pub cov: [[c_float; 3]; MAXSAT],     // 卫星位置协方差 (m^2)
    pub vco: [[c_float; 3]; MAXSAT],     // 卫星速度协方差 (m^2)
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Pclk {
    pub time: GTime,                   // 时间（GPST）
    pub index: c_int,                  // 多文件的时钟索引
    pub clk: [[c_double; 1]; MAXSAT],  // 卫星时钟 (s)
    pub std: [[c_float; 1]; MAXSAT],   // 卫星时钟标准差 (s)
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Alm {
    pub sat: c_int,          // 卫星编号
    pub svh: c_int,          // SV 健康状态（0: 正常）
    pub svconf: c_int,       // AS 和 SV 配置
    pub week: c_int,         // GPS/QZS: GPS 周，GAL: Galileo 周
    pub toa: GTime,          // Toa
    // SV 轨道参数
    pub A: c_double,         // 卫星轨道半长轴
    pub e: c_double,         // 卫星轨道离心率
    pub i0: c_double,        // 卫星轨道倾角
    pub OMG0: c_double,      // 卫星轨道升交点赤经
    pub omg: c_double,       // 卫星轨道近地点角距
    pub M0: c_double,        // 卫星轨道平近点角
    pub OMGd: c_double,      // 升交点赤经变化率
    pub toas: c_double,      // Toa（周内秒）
    pub f0: c_double,        // SV 时钟参数 (af0)
    pub f1: c_double,        // SV 时钟参数 (af1)
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Tec {
    pub time: GTime,               // 历元时间（GPST）
    pub ndata: [c_int; 3],         // TEC 网格数据大小 {nlat, nlon, nhgt}
    pub rb: c_double,              // 地球半径（km）
    pub lats: [c_double; 3],       // 纬度起点/间隔（度）
    pub lons: [c_double; 3],       // 经度起点/间隔（度）
    pub hgts: [c_double; 3],       // 高度起点/间隔（km）
    pub data: *mut c_double,       // TEC 网格数据（TECU）
    pub rms: *mut c_float,         // RMS 值（TECU）
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Erpd {
    pub mjd: c_double,        // MJD（天）
    pub xp: c_double,         // 极移偏移量（弧度）
    pub yp: c_double,         // 极移偏移量（弧度）
    pub xpr: c_double,        // 极移偏移率（弧度/天）
    pub ypr: c_double,        // 极移偏移率（弧度/天）
    pub ut1_utc: c_double,    // UT1-UTC（秒）
    pub lod: c_double,        // 日长变化（秒/天）
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Erp {
    pub n: c_int,                // 数据数量
    pub nmax: c_int,             // 最大数据数量
    pub data: *mut Erpd,         // 地球自转参数数据指针
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Pcv {
    pub sat: c_int,                        // 卫星编号（0:接收机）
    pub r#type: [c_char; MAXANT],          // 天线类型
    pub code: [c_char; MAXANT],            // 序列号或卫星代码
    pub ts: GTime,                         // 有效时间开始
    pub te: GTime,                         // 有效时间结束
    pub off: [[c_double; 3]; NFREQ],       // 相位中心偏移 e/n/u 或 x/y/z (m)
    pub var: [[c_double; 19]; NFREQ],      // 相位中心变化 (m)
    // el=90,85,...,0 或 天顶角=0,1,2,3,... (度)
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbsfcorr {
    pub t0: GTime,            // 适用时间（TOF）
    pub prc: c_double,        // 伪距修正（PRC）（米）
    pub rrc: c_double,        // 距离率修正（RRC）（米/秒）
    pub dt: c_double,         // 距离率修正时间增量（秒）
    pub iodf: c_int,          // IODF（快速修正日期问题）
    pub udre: c_short,      // UDRE+1
    pub ai: c_short,        // 衰减因子指示符
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbslcorr {
    pub t0: GTime,            // 修正时间
    pub iode: c_int,          // IODE（星历日期问题）
    pub dpos: [c_double; 3],  // 位置增量（米）（ECEF）
    pub dvel: [c_double; 3],  // 速度增量（米/秒）（ECEF）
    pub daf0: c_double,       // 时钟偏移增量（秒）
    pub daf1: c_double,       // 时钟漂移增量（秒/秒）
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbssatp {
    pub sat: c_int,           // 卫星编号
    pub fcorr: Sbsfcorr,      // 快速修正
    pub lcorr: Sbslcorr,      // 长期修正
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbssat {
    pub iodp: c_int,                      // IODP（掩码日期问题）
    pub nsat: c_int,                      // 卫星数量
    pub tlat: c_int,                      // 系统延迟（秒）
    pub sat: [Sbssatp; MAXSAT],           // 卫星修正
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbsigp {
    pub t0: GTime,        // 修正时间
    pub lat: c_short,     // 纬度（度）
    pub lon: c_short,     // 经度（度）
    pub give: c_short,    // GIVI+1
    pub delay: c_float,   // 垂直延迟估计（米）
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sbsion {
    pub iodi: c_int,                        // IODI（电离层修正日期问题）
    pub nigp: c_int,                        // IGP 数量
    pub igp: [Sbsigp; MAXNIGP],             // 电离层修正
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Sta {
    pub name: [c_char; MAXANT],
    pub marker: [c_char; MAXANT],
    pub antdes: [c_char; MAXANT],
    pub antsno: [c_char; MAXANT],
    pub rectype: [c_char; MAXANT],
    pub recver: [c_char; MAXANT],
    pub recsno: [c_char; MAXANT],
    pub antsetup: c_int,
    pub itrf: c_int,
    pub deltype: c_int,
    pub pos: [c_double; 3],
    pub del: [c_double; 3],
    pub hgt: c_double,
    pub glo_cp_align: c_int,
    pub glo_cp_bias: [c_double; 4],
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Dgps {
    pub t0: GTime,
    pub prc: c_double,
    pub rrc: c_double,
    pub iod: c_int,
    pub udre: c_double,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Ssr {
    pub t0: [GTime; 6],
    pub udi: [c_double; 6],
    pub iod: [c_int; 6],
    pub iode: c_int,
    pub iodcrc: c_int,
    pub ura: c_int,
    pub refd: c_int,
    pub deph: [c_double; 3],
    pub ddeph: [c_double; 3],
    pub dclk: [c_double; 3],
    pub hrclk: c_double,
    pub cbias: [c_float; MAXCODE as usize],
    pub pbias: [c_double; MAXCODE as usize],
    pub stdpb: [c_float; MAXCODE as usize],
    pub yaw_ang: c_double,
    pub yaw_rate: c_double,
    pub update: c_uchar,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct RtcmResC {
    pub r#type: c_int,
    pub staid: c_int,
    pub stah: c_int,
    pub seqno: c_int,
    pub time: GTime,
    pub time_s: GTime,
    pub obs: Obs,
    pub eph: Eph,
    pub geph: Geph,
    pub sta: Sta,
    pub dgps: *mut Dgps,
    pub ssr: [Ssr; MAXSAT],
    pub obsflag: c_int,
    pub ephsat: c_int,
    pub ephset: c_int,
    pub len: c_int,
}

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Nav {
    pub n: c_int,                             // 广播星历数量
    pub nmax: c_int,                          // 广播星历最大数量
    pub ng: c_int,                            // GLONASS 星历数量
    pub ngmax: c_int,                         // GLONASS 星历最大数量
    pub ns: c_int,                            // SBAS 星历数量
    pub nsmax: c_int,                         // SBAS 星历最大数量
    pub ne: c_int,                            // 精密星历数量
    pub nemax: c_int,                         // 精密星历最大数量
    pub nc: c_int,                            // 精密时钟数量
    pub ncmax: c_int,                         // 精密时钟最大数量
    pub na: c_int,                            // 星历数据数量
    pub namax: c_int,                         // 星历数据最大数量
    pub nt: c_int,                            // TEC 网格数据数量
    pub ntmax: c_int,                         // TEC 网格数据最大数量
    pub eph: *mut Eph,                        // GPS/QZS/GAL/BDS/IRN 星历
    pub geph: *mut Geph,                      // GLONASS 星历
    pub seph: *mut Seph,                      // SBAS 星历
    pub peph: *mut Peph,                      // 精密星历
    pub pclk: *mut Pclk,                      // 精密时钟
    pub alm: *mut Alm,                        // 星历数据
    pub tec: *mut Tec,                        // TEC 网格数据
    pub erp: Erp,                             // 地球自转参数
    pub utc_gps: [c_double; 8],               // GPS delta-UTC 参数
    pub utc_glo: [c_double; 8],               // GLONASS UTC 时间参数
    pub utc_gal: [c_double; 8],               // Galileo UTC 参数
    pub utc_qzs: [c_double; 8],               // QZS UTC 参数
    pub utc_cmp: [c_double; 8],               // 北斗 UTC 参数
    pub utc_irn: [c_double; 9],               // IRNSS UTC 参数
    pub utc_sbs: [c_double; 4],               // SBAS UTC 参数
    pub ion_gps: [c_double; 8],               // GPS 电离层模型参数
    pub ion_gal: [c_double; 4],               // Galileo 电离层模型参数
    pub ion_qzs: [c_double; 8],               // QZSS 电离层模型参数
    pub ion_cmp: [c_double; 8],               // 北斗电离层模型参数
    pub ion_irn: [c_double; 8],               // IRNSS 电离层模型参数
    pub glo_fcn: [c_int; 32],                 // GLONASS FCN + 8
    pub cbias: [[c_double; 3]; MAXSAT],       // 卫星 DCB (0:P1-P2,1:P1-C1,2:P2-C2)
    pub rbias: [[[c_double; 3]; 2]; MAXRCV],  // 接收机 DCB (0:P1-P2,1:P1-C1,2:P2-C2)
    pub pcvs: [Pcv; MAXSAT],                  // 卫星天线 pcv
    pub sbssat: Sbssat,                       // SBAS 卫星修正
    pub sbsion: [Sbsion; MAXBAND + 1],          // SBAS 电离层修正
    pub dgps: [Dgps; MAXSAT],                 // DGPS 修正
    pub ssr: [Ssr; MAXSAT],                   // SSR 修正
}


#[repr(C)]
#[derive(Debug, Clone)]
pub struct RtcmC {
    pub r#type: c_int,
    pub staid: c_int,
    pub stah: c_int,
    pub seqno: c_int,
    pub outtype: c_int,
    pub time: GTime,
    pub time_s: GTime,
    pub obs: Obs,
    pub nav: Nav,
    pub sta: Sta,
    pub dgps: *mut Dgps,
    pub ssr: [Ssr; MAXSAT],
    pub msg: [c_char; 128],
    pub msgtype: [c_char; 256],
    pub msmtype: [[c_char; 128]; 7],
    pub obsflag: c_int,
    pub ephsat: c_int,
    pub ephset: c_int,
    pub cp: [[c_double; NFREQ + NEXOBS]; MAXSAT],
    pub lock: [[c_ushort; NFREQ + NEXOBS]; MAXSAT],
    pub loss: [[c_ushort; NFREQ + NEXOBS]; MAXSAT],
    pub lltime: [[GTime; NFREQ + NEXOBS]; MAXSAT],
    pub nbyte: c_int,
    pub nbit: c_int,
    pub len: c_int,
    pub buff: [c_uchar; 1200],
    pub word: c_uint,
    pub nmsg2: [c_uint; 100],
    pub nmsg3: [c_uint; 400],
    pub opt: [c_char; 256],
}

impl From<Eph> for crate::basic::var::Eph {
    fn from(value: Eph) -> Self {
        crate::basic::var::Eph {
            sat: value.sat as usize,
            iode: value.iode,
            iodc: value.iodc,
            sva: value.sva,
            svh: value.svh,
            week: value.week,
            code: value.code,
            flag: value.flag,
            toe: value.toe.into(),
            toc: value.toc.into(),
            ttr: value.ttr.into(),
            a: value.A,
            e: value.e,
            i0: value.i0,
            omg0: value.OMG0,
            omg: value.omg,
            m0: value.M0,
            deln: value.deln,
            omgd: value.OMGd,
            idot: value.idot,
            crc: value.crc,
            crs: value.crs,
            cuc: value.cuc,
            cus: value.cus,
            cic: value.cic,
            cis: value.cis,
            toes: value.toes,
            fit: value.fit,
            f0: value.f0,
            f1: value.f1,
            f2: value.f2,
            tgd: {
                let mut new_tgd = [0.0; 6];
                for (i, &val) in value.tgd.iter().enumerate() {
                    new_tgd[i] = val as f64;
                }
                new_tgd
            },
            adot: value.Adot,
            ndot: value.ndot,
        }
    }
}

impl From<Geph> for crate::basic::var::Geph {
    fn from(value: Geph) -> Self {
        crate::basic::var::Geph {
            sat: value.sat as usize,
            iode: value.iode,
            frq: value.frq,
            svh: value.svh,
            sva: value.sva,
            age: value.age,
            toe: value.toe.into(),
            tof: value.tof.into(),
            pos: value.pos.map(|x| x as f64),
            vel: value.vel.map(|x| x as f64),
            acc: value.acc.map(|x| x as f64),
            taun: value.taun,
            gamn: value.gamn,
            dtaun: value.dtaun,
        }
    }
}

impl From<Obsd> for crate::basic::var::Obs {
    fn from(value: Obsd) -> Self {
        crate::basic::var::Obs {
            sat: value.sat as usize,
            time: value.time.into(),
            rcv: value.rcv,
            snr: value.SNR.map(|x| x as i32),
            lli: value.LLI.map(|x| x as i32),
            code: value.code.map(|x| x as u8),
            l: value.L.map(|x| x as f64),
            p: value.P.map(|x| x as f64),
            d: value.D.map(|x| x as f64),
        }
    }
}

impl From<Obs> for Obss {
    fn from(value: Obs) -> Self {
        let obs_vec = unsafe {
            slice::from_raw_parts(value.data, value.n as usize)
                .iter()
                .map(|&obsd| crate::basic::var::Obs::from(obsd))
                .collect()
        };

        Obss {
            n: value.n as usize,
            nmax: value.nmax as usize,
            obs: obs_vec,
        }
    }
}

fn cstr_to_string(cstr: &[c_char]) -> String {
    let cstr = unsafe { CStr::from_ptr(cstr.as_ptr()) };
    cstr.to_string_lossy().into_owned()
}

impl From<Sta> for crate::basic::var::Sta {
    fn from(sta_c: Sta) -> Self {
        crate::basic::var::Sta {
            name: cstr_to_string(&sta_c.name),
            marker: cstr_to_string(&sta_c.marker),
            antdes: cstr_to_string(&sta_c.antdes),
            antsno: cstr_to_string(&sta_c.antsno),
            rectype: cstr_to_string(&sta_c.rectype),
            recver: cstr_to_string(&sta_c.recver),
            recsno: cstr_to_string(&sta_c.recsno),
            antsetup: sta_c.antsetup,
            itrf: sta_c.itrf,
            deltype: sta_c.deltype,
            pos: sta_c.pos,
            del: sta_c.del,
            hgt: sta_c.hgt,
            glo_cp_align: sta_c.glo_cp_align,
            glo_cp_bias: sta_c.glo_cp_bias,
        }
    }
}






