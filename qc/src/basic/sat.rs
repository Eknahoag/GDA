use crate::basic::var::*;

/// convert satellite system+prn to satellite number
pub fn satno(sys: usize, prn: usize) -> usize {
    if prn <= 0 {
        return 0;
    }

    match sys {
        SYS_GPS if (MINPRNGPS..=MAXPRNGPS).contains(&prn) => prn,
        SYS_GLO if (MINPRNGLO..=MAXPRNGLO).contains(&prn) => NSATGPS + prn,
        SYS_GAL if (MINPRNGAL..=MAXPRNGAL).contains(&prn) => NSATGPS + NSATGLO + prn,
        SYS_QZS if (MINPRNQZS..=MAXPRNQZS).contains(&prn) => NSATGPS + NSATGLO + NSATGAL + prn - MINPRNQZS + 1,
        SYS_CMP if (MINPRNCMP..=MAXPRNCMP).contains(&prn) => {
            NSATGPS + NSATGLO + NSATGAL + NSATQZS + prn
        }
        SYS_IRN if (MINPRNIRN..=MAXPRNIRN).contains(&prn) => {
            NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + prn
        }
        _ => 0,
    }
}

/// convert satellite number to PRN
pub fn getprn(sat: usize) -> usize {
    let mut prn = sat;
    if prn <= 0 {
        0
    } else if prn <= NSATGPS {
        return prn;
    } else if { prn -= NSATGPS; prn } <= NSATGLO {
        return prn;
    }else if {prn-=NSATGLO;prn}<=NSATGAL {
        return prn;
    }else if {prn-=NSATGAL;prn}<=NSATQZS {
        return prn;
    }else if {prn-=NSATQZS;prn}<=NSATCMP {
        return prn;
    }else if {prn-=NSATCMP;prn}<=NSATIRN {
        return prn;
    }
    else {
        return 0;
    }
}

/// convert satellite number to satellite id
pub fn satno2id(sat: usize) -> String {
    let prn=getprn(sat);
    match satsys(sat) {
        SYS_GPS => format!("G{:02}", prn),
        SYS_GLO => format!("R{:02}", prn),
        SYS_GAL => format!("E{:02}", prn),
        SYS_QZS => format!("J{:02}", prn),
        SYS_CMP => format!("C{:02}", prn),
        SYS_IRN => format!("I{:02}", prn),
        _ => String::new(),
    }
}

/// convert satellite id to satellite number
pub fn satid2no(id: &str) -> usize {
    let sys: usize;

    if let Ok(num) = id.parse::<usize>() {
        // Assume MIN/MAX constants are defined somewhere
        if (MINPRNGPS..=MAXPRNGPS).contains(&num) {
            sys = SYS_GPS;
        } else if (MINPRNQZS..=MAXPRNQZS).contains(&num) {
            sys = SYS_QZS;
        } else {
            return 0;
        }
        return satno(sys, num);
    }

    let code = id.chars().next().unwrap();
    let mut prn = id[1..].parse::<usize>().unwrap();

    sys = match code {
        'G' => SYS_GPS,
        'R' => SYS_GLO,
        'E' => SYS_GAL,
        'J' => {
            prn += MINPRNQZS - 1;
            SYS_QZS
        }
        'C' => SYS_CMP,
        'I' => SYS_IRN,
        _ => return 0,
    };
    satno(sys, prn)
}

/// get satellite system
pub fn satsys(sat: usize) -> usize {
    if sat <= 0 {
        SYS_NONE
    } else if sat <= NSATGPS {
        SYS_GPS
    } else if sat <= NSATGPS + NSATGLO {
        SYS_GLO
    } else if sat <= NSATGPS + NSATGLO + NSATGAL {
        SYS_GAL
    } else if sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS {
        SYS_QZS
    } else if sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP {
        SYS_CMP
    } else if sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN {
        SYS_IRN
    } else {
        SYS_NONE
    }
}