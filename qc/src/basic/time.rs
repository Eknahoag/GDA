use crate::basic::var::*;
use chrono::{FixedOffset, Datelike, NaiveDate, NaiveDateTime, NaiveTime, Timelike, Utc};
use std::error::Error;
use std::time::{SystemTime, UNIX_EPOCH};

const GPST0: [f64; 6] = [1980.0, 1.0, 6.0, 0.0, 0.0, 0.0];
const BDT0: [f64; 6] = [2006.0, 1.0, 1.0, 0.0, 0.0, 0.0];
pub fn timestr_rnx() -> String {
    let now = Utc::now();
    let formatted_time = format!(
        "{:04}{:02}{:02} {:02}{:02}{:02} UTC",
        now.year(),
        now.month(),
        now.day(),
        now.hour(),
        now.minute(),
        now.second()
    );
    formatted_time
}

pub fn timeget() -> i64 {
    let system_time = SystemTime::now();
    let t = system_time.duration_since(UNIX_EPOCH).expect("SystemTime before UNIX EPOCH!");
    t.as_secs() as i64
}

pub fn timeadd(mut t: GTime, sec: f64) -> GTime {
    t.sec += sec;
    let tt = t.sec.floor();
    t.time += tt as i64;
    t.sec -= tt;
    t
}

pub fn timediff(t1: GTime, t2: GTime) -> f64 {
    (t1.time - t2.time) as f64 + (t1.sec - t2.sec)
}

pub fn adjweek(t: GTime, t0: GTime) -> GTime {
    let tt = timediff(t, t0);
    if tt < -302400.0 {
        return timeadd(t, 604800.0);
    }
    if tt > 302400.0 {
        return timeadd(t, -604800.0);
    }
    t
}

pub fn adjday(t: GTime, t0: GTime) -> GTime {
    let tt = timediff(t, t0);
    if tt < -43200.0 {
        return timeadd(t, 86400.0);
    }
    if tt > 43200.0 {
        return timeadd(t, -86400.0);
    }
    return t;
}

pub fn bdt2gpst(t: GTime) -> GTime {
    return timeadd(t, 14.0);
}

pub fn gpst2bdt(t: GTime) -> GTime {
    return timeadd(t, -14.0);
}

pub fn bdt2time(week: i32, sec: f64) -> GTime {
    let mut t = epoch2time(&BDT0);
    let mut s = sec;
    if sec < -1E9 || 1E9 < sec {
        s = 0.0;
    }
    t.time += 86400 * 7 * week as i64 + s as i64;
    t.sec = s - s.floor();
    return t;
}

pub fn time2bdt(t: GTime, week: &mut i32) -> f64 {
    let t0 = epoch2time(&BDT0);
    let sec = t.time - t0.time;
    let w = sec / 86400 / 7;
    if *week != 0 { *week = w as i32 }
    (sec - w * 86400 * 7) as f64 + t.sec
}

pub fn gpst2utc(t: GTime) -> GTime {
    let mut tu: GTime;
    for i in 0..LEAPS.len() {
        if LEAPS[i][0] <= 0.0 {
            break;
        }
        tu = timeadd(t, LEAPS[i][6]);
        if timediff(tu, epoch2time(&LEAPS[i][..6].try_into().unwrap())) >= 0.0 {
            return tu;
        }
    }
    t
}

pub fn gpst2time(week: i32, sec: f64) -> GTime {
    let mut t = epoch2time(&GPST0);
    let mut s = sec;
    if sec < -1E9 || 1E9 < sec {
        s = 0.0;
    }
    t.time += 86400 * 7 * week as i64 + s as i64;
    t.sec = s - s.floor();
    return t;
}

pub fn utc2gpst(t: GTime) -> GTime {
    for leap in LEAPS.iter() {
        if leap[0] <= 0.0 {
            break;
        }
        if timediff(t, epoch2time(&leap[..6].try_into().unwrap())) >= 0.0 {
            return timeadd(t, -leap[6]);
        }
    }
    t
}

pub fn time2gpst(t: GTime, week: Option<&mut i32>) -> f64 {
    let t0 = epoch2time(&GPST0);
    let sec = t.time - t0.time;
    let w = (sec / (86400 * 7)) as i32;

    if let Some(week_ref) = week {
        *week_ref = w;
    }
    (sec - (w as i64 * 86400 * 7)) as f64 + t.sec
}

pub fn str2time(s: &str, t: &mut GTime) -> Result<(), Box<dyn Error>> {
    let mut parts: Vec<String> = s.split_whitespace().map(|s| s.to_string()).collect();
    // 检查年份是否只有两位数
    if parts[0].len() == 2 {
        // 将两位数年份转换为四位数，假设21世纪
        parts[0] = format!("20{}", parts[0]);
    }
    // 重新组合字符串
    let ss = parts.join(" ");

    let date_time = NaiveDateTime::parse_from_str(&ss, "%Y %m %d %H %M %S%.f")?;

    // Convert NaiveDateTime to seconds since epoch
    let epoch_date = NaiveDate::from_ymd_opt(1970, 1, 1).ok_or("Invalid epoch date")?;
    let epoch_time = NaiveTime::from_hms_opt(0, 0, 0).ok_or("Invalid epoch time")?;
    let epoch = NaiveDateTime::new(epoch_date, epoch_time);

    let duration = date_time.signed_duration_since(epoch);

    t.time = duration.num_seconds();
    t.sec = 0.0; // Assuming the input format does not have fractional seconds

    Ok(())
}

pub fn time2epoch(t: GTime, ep: &mut [f64; 6]) {
    const MDAY: [i32; 48] = [
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30,
        31, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31,
        30, 31,
    ];

    let days = (t.time / 86400) as i32;
    let sec = (t.time - (days as i64 * 86400)) as i32;
    let mut day = days % 1461;
    let mut mon = 0;

    while mon < 48 {
        if day >= MDAY[mon] {
            day -= MDAY[mon];
        } else {
            break;
        }
        mon += 1;
    }

    ep[0] = 1970.0 + (days / 1461 * 4) as f64 + (mon / 12) as f64;
    ep[1] = (mon % 12 + 1) as f64;
    ep[2] = (day + 1) as f64;
    ep[3] = (sec / 3600) as f64;
    ep[4] = (sec % 3600 / 60) as f64;
    ep[5] = (sec % 60) as f64 + t.sec;
}

pub fn epoch2time(ep: &[f64; 6]) -> GTime {
    const DOY: [i32; 12] = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];

    let mut time = GTime { time: 0, sec: 0.0 };
    let year = ep[0] as i32;
    let mon = ep[1] as i32;
    let day = ep[2] as i32;

    if year < 1970 || year > 2099 || mon < 1 || mon > 12 {
        return time;
    }

    // leap year if year%4==0 in 1901-2099
    let days = (year - 1970) * 365 + (year - 1969) / 4 + DOY[(mon - 1) as usize] + day - 2
        + if year % 4 == 0 && mon >= 3 { 1 } else { 0 };

    let sec = ep[5].floor() as i32;
    time.time = (days as i64) * 86400 + (ep[3] as i64) * 3600 + (ep[4] as i64) * 60 + sec as i64;
    time.sec = ep[5] - sec as f64;

    time
}

pub fn time2str(t: GTime) -> String {
    let mut ep: [f64; 6] = [0.0; 6];
    let mut t = t;

    if 1.0 - t.sec < 0.5 {
        t.time += 1;
        t.sec = 0.0;
    }

    time2epoch(t, &mut ep);

    let formatted_str = format!(
        "{:04}/{:02}/{:02} {:02}:{:02}:{:02}",
        ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
    );

    formatted_str
}

pub fn screent(time: GTime, ts: GTime, te: GTime, tint: f64) -> bool {
    let cond_tint = tint <= 0.0 || (time2gpst(time, None) + DTTOL) % tint <= DTTOL * 2.0;
    let cond_ts = ts.time == 0 || timediff(time, ts) >= -DTTOL;
    let cond_te = te.time == 0 || timediff(time, te) < DTTOL;

    cond_tint && cond_ts && cond_te
}

// 判断是否为北京时间零点
pub fn bjt_0() -> bool {
    // 获取当前UTC时间
    let now_utc = Utc::now();

    // 将UTC时间转换为北京时间（东八区）
    let beijing_tz = FixedOffset::east_opt(8 * 3600).expect("invalid time zone set");
    let now_beijing = now_utc.with_timezone(&beijing_tz);

    // 检查当前北京时间是否为零点
    if now_beijing.hour() == 0 && now_beijing.minute() == 0 && now_beijing.second() == 0 {
        return true;
    }
    false
}

pub fn date2gtime(date: &str) -> GTime {
    let mut epe = [0.0; 6];

    if let Some((date_str, time_str)) = date.split_once(' ') {
        date_str
            .split('/')
            .zip(&mut epe[..3])
            .for_each(|(s, e)| *e = s.parse().unwrap_or(0.0));
        time_str
            .split(':')
            .zip(&mut epe[3..])
            .for_each(|(s, e)| *e = s.parse().unwrap_or(0.0));
    }

    epoch2time(&epe)
}
