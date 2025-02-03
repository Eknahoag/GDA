#![allow(unused_imports)]
#![allow(dead_code)]

use std::{env, io};
use std::fs::File;
use std::io::BufRead;
use std::path::Path;
use clap::{Parser, Subcommand};
use gda::basic::ntrip::ntripqc;
use gda::basic::time::{date2gtime, str2time, timeget};
use gda::basic::var::{GTime, MAXINFILE, Nav, Obss, QcOpt, RnxOpt, Sta, SYS_CMP, SYS_GAL, SYS_GLO, SYS_GPS, SYS_IRN, SYS_QZS};
use gda::convbin::rnxout::{ntrip2rtcm, rtcm2rnx};
use gda::qc::qc_fun::{obsqc, rtcmqc};
use gda::convbin::decode::Decoder;
use gda::basic::read::readobsnav;
use gda::qc::cors::con_zmq;
use gda::qc::qc_var::{Outtype, QC, RedisParams};

#[derive(Parser, Debug)]
#[command(
    name = "GDA",
    about = "A GNSS toolbox for data processing like decoding and quality analysis.",
    after_help = "If you have any questions, please reach out to: ghk40041@whu.edu.cn"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Decode RTCM3 to RINEX
    Decode {
        /// Input RTCM3 File
        #[clap(short = 'i', long = "inp")]
        ifile: String,

        /// Output OBS File
        #[clap(short = 'o', long = "obs")]
        obsfile: Option<String>,

        /// Output NAV File
        #[clap(short = 'n', long = "nav")]
        navfile: Option<String>,

        /// Time of Start
        #[clap(short = 's', long = "ts")]
        ts: Option<String>,

        /// Time of End
        #[clap(short = 'e', long = "te")]
        te: Option<String>,

        /// Sample
        #[clap(short = 't', long = "ti")]
        ti: Option<f64>,

        /// Approximate Time
        #[clap(short = 'r', long = "tr")]
        tr: Option<String>,

        /// Output File Version
        #[clap(short = 'v', long = "version",default_value="3.04")]
        version: f64,

        /// System Selection
        #[clap(short = 'y', long = "sys", use_value_delimiter = true)]
        navsys: Option<Vec<char>>,
    },

    /// Store Bin Data to File
    Binstore {
        /// NTRIP Server Hosting
        #[clap(short = 's', long = "host")]
        host: String,

        /// NTRIP Server Port
        #[clap(short = 'p', long = "port")]
        port: String,

        /// NTRIP Mount Point
        #[clap(short = 'm', long = "mountpoint")]
        mountpoint: String,

        /// NTRIP Username
        #[clap(short = 'u', long = "username")]
        username: String,

        /// NTRIP Password
        #[clap(short = 'w', long = "password")]
        password: String,

        /// Output File, defaults to the name of the mount point in the project directory.
        #[clap(short = 'o', long = "oup")]
        outfile: Option<String>,
    },

    /// QC by NTRIP Message to *File or Redis
    Ntrip {
        /// NTRIP Server Hosting
        #[clap(short = 's', long = "host")]
        host: String,

        /// NTRIP Server Port
        #[clap(short = 'p', long = "port")]
        port: String,

        /// NTRIP Mount Point
        #[clap(short = 'm', long = "mountpoint")]
        mountpoint: String,

        /// NTRIP Username
        #[clap(short = 'u', long = "username")]
        username: String,

        /// NTRIP Password
        #[clap(short = 'w', long = "password")]
        password: String,

        /// QC File Output Path, saved in the current directory by default with the name of the mount point
        #[clap(short = 'q', long = "qcpath")]
        qcpath: Option<String>,

        /// Redis Server Address
        #[clap(short = 'i', long = "redisip")]
        redisip: Option<String>,

        /// Redis Server Port
        #[clap(short = 'o', long = "redisport")]
        redisport: Option<String>,

        /// Redis Server Authentication Commands
        #[clap(short = 'a', long = "auth")]
        redisauth: Option<String>,

        /// Additional Documents for Ephemeris Supplement
        #[clap(short = 'd', long = "addition")]
        addition: Option<String>,
    },
    /// QC by ZMQ to File or *Redis
    Zmq {
        /// Subscription Address
        #[clap(short = 's', long = "host")]
        host: String,
        /// Server Port
        #[clap(short = 'p', long = "port")]
        port: String,
        /// Site Configuration
        #[clap(short = 't', long = "stations", use_value_delimiter = true)]
        station: Option<Vec<String>>,
        /// Redis Server Address
        #[clap(short = 'i', long = "redisip")]
        redisip: Option<String>,
        /// Redis Server Port
        #[clap(short = 'o', long = "redisport")]
        redisport: Option<String>,
        /// Redis Server Authentication Commands
        #[clap(short = 'a', long = "auth")]
        redisauth: Option<String>,
        /// Default Output Path for Analysis Files
        #[clap(short = 'q', long = "qcpath")]
        qcpath: Option<String>,
        /// Supplementary ephemeris
        #[clap(short = 'd', long = "addition")]
        addition: Option<String>,
    },
    /// QC by RINEX File
    Obsqc {
        /// Input file
        #[clap(short = 'i', long = "input", use_value_delimiter = true)]
        ifile: Vec<String>,

        /// Output file
        #[clap(short = 'o', long = "output")]
        ofile: String,

        /// Time of Start
        #[clap(short = 's', long = "ts")]
        ts: Option<String>,

        /// Time of End
        #[clap(short = 'e', long = "te")]
        te: Option<String>,

        /// Approximate Time
        #[clap(short = 't', long = "ti")]
        ti: Option<f64>,

        /// Output Details (0, 1, 2)
        #[clap(short = 'q', long = "detail", default_value = "0")]
        q: i32,
    },
    /// QC by RTCM3 File
    Rtcmqc {
        /// Input File
        #[clap(short = 'i', long = "input")]
        ifile: String,

        /// Output File
        #[clap(short = 'o', long = "output")]
        ofile: String,

        /// Time of Start
        #[clap(short = 's', long = "ts")]
        ts: Option<String>,

        /// Time of End
        #[clap(short = 'e', long = "te")]
        te: Option<String>,

        /// Sample
        #[clap(short = 't', long = "ti")]
        ti: Option<f64>,

        /// Approximate Time
        #[clap(short = 'r', long = "tr")]
        tr: Option<String>,

        /// Output Details (0, 1, 2)
        #[clap(short = 'q', long = "detail", default_value = "0")]
        q: i32,
    },
}

#[tokio::main]
async fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Decode {
            ifile,
            obsfile,
            navfile,
            ts,
            te,
            ti,
            tr,
            version,
            navsys
        } => {
            if ifile.is_empty() {
                eprintln!("Please set input file!");
                return;
            }
            if obsfile.is_none() && navfile.is_none() {
                eprintln!("Please set output file!");
                return;
            }
            println!("input file:{}", ifile);
            if obsfile.is_some() { println!("output file:{}", obsfile.clone().unwrap()) }
            if navfile.is_some() { println!("output file:{}", navfile.clone().unwrap()) }

            let mut opt = QcOpt::default();

            opt.tr = if tr.is_none() {
                GTime {
                    time: timeget(),
                    sec: 0.0,
                }
            } else {
                date2gtime(&tr.unwrap())
            };
            if ts.is_some() { opt.ts = date2gtime(&ts.unwrap()) }
            if te.is_some() { opt.te = date2gtime(&te.unwrap()) }
            if ti.is_some() { opt.ti = ti.unwrap() }

            if let Some(sys)=navsys{
                for value in sys  {
                    match value {
                        'G'=>opt.navsys|=SYS_GPS,
                        'R'=>opt.navsys|=SYS_GLO,
                        'E'=>opt.navsys|=SYS_GAL,
                        'C'=>opt.navsys|=SYS_CMP,
                        'Q'=>opt.navsys|=SYS_QZS,
                        'I'=>opt.navsys|=SYS_IRN,
                        _=>continue
                    }
                }
            }else {
                opt.navsys=SYS_GPS | SYS_GLO | SYS_GAL | SYS_CMP | SYS_QZS
            }
            opt.rnxver=(version*100.0) as i32;

            let _=rtcm2rnx(&ifile,obsfile,navfile,opt);
        }

        Commands::Binstore {
            host,
            port,
            mountpoint,
            username,
            password,
            outfile,
        }=>{
            let file = outfile.clone().unwrap_or_else(|| {
                let current_dir = env::current_dir().expect("Unable to get the current directory");
                format!("{}\\{}",current_dir.to_str().unwrap(),mountpoint).to_string()
            });
            if host.is_empty() { println!("Please input host");return; }
            if port.is_empty() { println!("Please input port");return; }
            if mountpoint.is_empty() { println!("Please input mountpoint");return; }

            ntrip2rtcm(&host,&port,&mountpoint,&username,&password,&file).await
        }

        Commands::Ntrip {
            host,
            port,
            mountpoint,
            username,
            password,
            qcpath,
            redisip,
            redisport,
            redisauth,
            addition,
        } => {
            let path = qcpath.clone().unwrap_or_else(|| {
                let current_dir = env::current_dir().expect("Unable to get the current directory");
                current_dir.to_str().unwrap().to_string()
            });

            if host.is_empty() { println!("Please input host");return; }
            if port.is_empty() { println!("Please input port");return; }
            if mountpoint.is_empty() { println!("Please input mountpoint");return; }

            // Set the output mode, preferentially file output
            let outtype = if redisip.is_none() {
                Outtype::File(path)
            } else if qcpath.is_none() {
                Outtype::Redis(RedisParams { ip: redisip, port: redisport, auth: redisauth })
            } else {
                Outtype::All(path, RedisParams { ip: redisip, port: redisport, auth: redisauth })
            };

            let mut supple_nav = None;

            // Set the supplemental ephemeris file to allow replacement if ephemeris data is not received
            if addition.is_some() {
                let mut obss: Obss = Obss::new();
                let mut navs = Nav::new();
                let mut sta = Sta::default();
                let mut nepoch = 0;
                let index = [0; MAXINFILE];
                //let ifile=vec![Arg::parse().addition.unwrap()];
                let ifile = vec![String::from(addition.unwrap())];
                let t0 = GTime::default();
                let _ = readobsnav(
                    t0,
                    t0,
                    0.0,
                    ifile,
                    index,
                    &mut obss,
                    &mut navs,
                    &mut sta,
                    &mut nepoch,
                );
                if navs.n + navs.ng > 0 { supple_nav = Some(navs) };
            }

            ntripqc(&host, &port, &mountpoint, &username, &password, supple_nav, outtype).await
        }

        Commands::Zmq {
            host,
            port,
            station,
            redisip,
            redisport,
            redisauth,
            qcpath,
            addition,
        } => {
            let path = qcpath.clone().unwrap_or_else(|| {
                let current_dir = env::current_dir().expect("Unable to get the current directory");
                current_dir.to_str().unwrap().to_string()
            });
            if host.is_empty() { println!("Please input host");return; }
            if port.is_empty() { println!("Please input port");return; }

            // Set the output mode, with redis output preferred
            let outtype = if qcpath.is_none() {
                Outtype::Redis(RedisParams { ip: redisip, port: redisport, auth: redisauth })
            } else if redisip.is_none() {
                Outtype::File(path.clone())
            } else {
                Outtype::All(path.clone(), RedisParams { ip: redisip, port: redisport, auth: redisauth })
            };

            let mut supple_nav = None;

            if addition.is_some() {
                let mut obss: Obss = Obss::new();
                let mut navs = Nav::new();
                let mut sta = Sta::default();
                let mut nepoch = 0;
                let index = [0; MAXINFILE];
                //let ifile=vec![Arg::parse().addition.unwrap()];
                let ifile = vec![String::from(addition.unwrap())];
                let t0 = GTime::default();
                let _ = readobsnav(
                    t0,
                    t0,
                    0.0,
                    ifile,
                    index,
                    &mut obss,
                    &mut navs,
                    &mut sta,
                    &mut nepoch,
                );
                if navs.n + navs.ng > 0 { supple_nav = Some(navs) };
            }

            con_zmq(&host, &port, station, supple_nav, outtype).await
        }

        Commands::Obsqc {
            ifile,
            ofile,
            ts,
            te,
            ti,
            q,
        } => {
            let mut opt = QcOpt::default();
            if ifile.len() == 0 {
                eprintln!("Please set input file!");
                return;
            } else if ifile.len() > 2 {
                eprintln!("Too many input file!");
                return;
            }
            if ofile.is_empty() {
                eprintln!("Please set output file!");
                return;
            }
            for i in 0..ifile.len() {
                println!("input file:{}", ifile[i]);
            }
            println!("output file:{}", ofile);

            if ts.is_some() { opt.ts = date2gtime(&ts.unwrap()) }
            if te.is_some() { opt.te = date2gtime(&te.unwrap()) }
            if ti.is_some() { opt.ti = ti.unwrap() }
            opt.detail = q;

            obsqc(&mut opt, &ifile, &ofile)
        }
        Commands::Rtcmqc {
            ifile,
            ofile,
            ts,
            te,
            ti,
            tr,
            q,
        } => {
            let mut opt = QcOpt::default();
            if ifile.is_empty() {
                eprintln!("Please set input file!");
                return;
            }
            if ofile.is_empty() {
                eprintln!("Please set output file!");
                return;
            }
            println!("input file:{}", ifile);
            println!("output file:{}", ofile);

            if ts.is_some() { opt.ts = date2gtime(&ts.unwrap()) }
            if te.is_some() { opt.te = date2gtime(&te.unwrap()) }
            if ti.is_some() { opt.ti = ti.unwrap() }
            opt.tr = if tr.is_none() {
                GTime {
                    time: timeget(),
                    sec: 0.0,
                }
            } else {
                date2gtime(&tr.unwrap())
            };
            opt.detail = q;

            rtcmqc(&mut opt, &ifile, &ofile);
        }
    }
}