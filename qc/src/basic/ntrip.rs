use std::fs::{self, File};
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;
use tokio::net::TcpStream;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::sync::mpsc;
use tokio::time::{sleep, Duration, timeout};
use std::time::Instant;
use base64::{engine::general_purpose, Engine as _};
use chrono::prelude::*;
use crate::basic::time::{bjt_0, time2str};
use crate::basic::var::{ MAXSAT, Nav, Obss, PrcOpt, Rtk};
use crate::convbin::decode::Decoder;
use crate::qc::qc_fun::{qc_cal_rt, qc_obs_rt, qc_out_rt, rt_rtcm2qc};
use crate::qc::qc_var::{Outtype, QC};

const RTCM_TYPE: &[u16] = &[
    //OBS L1 L1/L2
    1002, 1004, 1010, 1012,
    //STA POSITION
    1005, 1006,
    //ANT/RCV INFO
    1007, 1008, 1033,
    //NAV
    1019, 1020, 1029, 1041, 1044, 1045, 1046, 1042, 63,
    // SSR
    1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068,
    1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251,
    1258, 1259, 1260, 1261, 1262, 1263,
    1252, 1253, 1254, 1255, 1256, 1257,
    11, 12, 13, 14,
    //GPS MSM
    1074, 1075, 1076, 1077,
    //GLO MSM
    1084, 1085, 1086, 1087,
    //GAL MSM
    1094, 1095, 1096, 1097,
    //SBS MSM
    1104, 1105, 1106, 1107,
    //QZS MSM
    1114, 1115, 1116, 1117,
    //BDS MSM
    1124, 1125, 1126, 1127,
    //IRN MSM
    1134, 1135, 1136, 1137,
    //PROPRIETARY
    4076,
    1230,
    4073,
];

#[allow(dead_code)]
// Functions to convert byte arrays to hexadecimal strings
fn bytes_to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x} ", b)).collect()
}
#[allow(dead_code)]
// Creates a real-time QC (Quality Control) file with a timestamp based on the provided mount point and QC path.
fn create_realtime_qcfile(mountpoint: &str, qcpath: &str) -> String {
    if mountpoint.is_empty() || qcpath.is_empty() {
        return String::from("");
    }

    // Get the current UTC time
    let now_utc = Utc::now();

    // Convert UTC time to Beijing time
    let bjt = now_utc.with_timezone(&FixedOffset::east_opt(8 * 3600).expect("invalid time zone set"));
    let current_time = bjt.format("%Y%m%d_%H%M%S").to_string();

    // Generate the filename
    let filename = format!("{}_{}.qc", mountpoint, current_time);

    // Check if the path exists, create if it does not exist
    let path = Path::new(qcpath);
    if !path.exists() {
        fs::create_dir_all(path).expect("Failed to create directory");
    }

    // Generate the full path and filename
    let full_filename = format!("{}/{}", qcpath, filename);

    full_filename
}


fn create_qcfile(mountpoint: &str, qcpath: &str) -> String {
    // Check if either mountpoint or qcpath is empty, return an empty string if true
    if mountpoint.is_empty() || qcpath.is_empty() {
        return String::from("");
    }

    // Generate the filename
    let filename = format!("{}.qc", mountpoint);
    // Check if the path exists, create it if it does not exist
    let path = Path::new(qcpath);
    if !path.exists() {
        fs::create_dir_all(path).expect("Failed to create directory");
    }
    // Generate the full path and filename
    let full_filename = format!("{}/{}", qcpath, filename);

    full_filename
}

/// Extracts RTCM messages from the given byte slice.
///
/// # Arguments
/// * `bytes` - A byte slice containing RTCM messages.
///
/// # Returns
/// A vector of vectors, where each inner vector contains the bytes of a single RTCM message.
///
/// # Remarks
/// This function iterates over the input byte slice to locate and extract RTCM messages that match specific types.
/// The message type and length are parsed according to the RTCM message specification.
pub fn get_single_rtcm(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut res = Vec::new();
    let mut i = 0;

    while i < bytes.len() {
        if bytes[i] == 0xd3 {
            if i + 4 < bytes.len() {
                let type_ = (((bytes[i + 3] as u16) << 8) | (bytes[i + 4] as u16)) >> 4;
                if !RTCM_TYPE.contains(&type_) {
                    i += 1;
                    continue;
                }
                let len = ((((bytes[i + 1] as u16) << 8) | (bytes[i + 2] as u16)) & 0x3FF) as usize + 3 + 3;
                if i + len <= bytes.len() {
                    let segment = &bytes[i..=i + len - 1];
                    res.push(segment.to_vec());
                    i += len;
                } else { break; }
            }
        }
        i += 1;
    }

    res
}

/// Output analysis results of ntrip real-time data
pub fn ntrip_out(
    qc: &mut QC,
    rtk: &mut Rtk,
    popt: &PrcOpt,
    obss: &Obss,
    navs: &Nav,
    qcfile: &str,
) {
    let navflag = if navs.n == 0 || navs.ng == 0 { false } else { true };
    qc_obs_rt(qc, rtk, &popt, &obss.obs, navs);
    qc_cal_rt(qc, popt.navsys);
    if !qcfile.is_empty() {
        let file = File::create(qcfile).expect("Fail to create output file!");
        let mut writer = BufWriter::new(file);
        qc_out_rt(qc, popt.navsys, navflag, &rtk, &mut writer);
    }
    print!("\rProcessing at {}", time2str(obss.obs[0].time));
    io::stdout().flush().unwrap();
    for i in 0..MAXSAT {
        qc.satinfo[i].code = [0; 6];
    }
}

/// Connect to NTRIP server
pub async fn conntrip(
    host: &str,
    port: &str,
    mountpoint: &str,
    username: &str,
    password: &str,
) -> Option<TcpStream> {
    let server = format!("{}:{}", host, port);
    let auth_value = format!(
        "Basic {}",
        general_purpose::STANDARD.encode(format!("{}:{}", username, password))
    );

    let mut retries = 0;
    let max_retries = 10;
    let connect_timeout_secs = 60;
    let retry_delay_secs = 1;

    // HTTP/1.1 and Ntrip 2.0 by default, Connection: close, fallback on demand
    let mut use_http11 = true;
    let mut use_ntrip2 = true;
    let mut use_keep_alive = false; // 默认 Connection: close

    while retries < max_retries {
        retries += 1;
        let start_time = Instant::now();

        // Set HTTP, Ntrip-Version, and Connection fields based on current state
        let http_version = if use_http11 { "1.1" } else { "1.0" };
        let ntrip_version = if use_ntrip2 { "Ntrip/2.0" } else { "Ntrip/1.0" };
        let connection_type = if use_keep_alive { "keep-alive" } else { "close" };

        let request = format!(
            "GET /{} HTTP/{}\r\n\
             Host: {}\r\n\
             Authorization: {}\r\n\
             Ntrip-Version: {}\r\n\
             Connection: {}\r\n\
             \r\n",
            mountpoint, http_version, host, auth_value, ntrip_version, connection_type
        );

        println!(
            "Attempting to connect to NTRIP server: {} with HTTP/{} and Ntrip-Version: {}, Connection: {} (Attempt {}/{})",
            server, http_version, ntrip_version, connection_type, retries, max_retries
        );

        match timeout(Duration::from_secs(connect_timeout_secs), TcpStream::connect(&server)).await {
            Ok(Ok(mut stream)) => {
                println!("Connected to server: {}", server);

                // Send Request
                if let Err(e) = stream.write_all(request.as_bytes()).await {
                    eprintln!("Failed to send request: {}", e);
                    return None;
                }

                // Validate server response, read response header
                let mut response = vec![0; 1024];
                if let Ok(n) = stream.read(&mut response).await {
                    let resp_str = String::from_utf8_lossy(&response[..n]);
                    if resp_str.contains("200 OK") {
                        println!(
                            "Server response OK with HTTP/{} and Ntrip-Version: {}, Connection: {}",
                            http_version, ntrip_version, connection_type
                        );
                        return Some(stream);
                    } else {
                        eprintln!("Unexpected server response: {}", resp_str);

                        // Adapt retry logic to different conditions
                        if use_http11 {
                            println!("Falling back to HTTP/1.0...");
                            use_http11 = false;
                        } else if use_ntrip2 {
                            println!("Falling back to Ntrip/1.0...");
                            use_ntrip2 = false;
                        } else if !use_keep_alive {
                            println!("Trying Connection: keep-alive...");
                            use_keep_alive = true;
                        } else {
                            break; // If all fallback options fail, exit the loop
                        }
                        continue;
                    }
                } else {
                    eprintln!("Failed to read server response");
                }

                return None;
            }
            Ok(Err(e)) => {
                eprintln!("Failed to connect (attempt {}): {}. Retrying...", retries, e);
            }
            Err(_) => {
                eprintln!("Connection timed out (attempt {}). Retrying...", retries);
            }
        }

        // Check if the connection attempt times out
        if start_time.elapsed().as_secs() > connect_timeout_secs {
            eprintln!("Connection attempt took too long. Aborting...");
            return None;
        }

        // Retry after delay
        if retries < max_retries {
            println!("Waiting {} seconds before next attempt...", retry_delay_secs);
            sleep(Duration::from_secs(retry_delay_secs)).await;
        } else {
            eprintln!("Reached maximum retry attempts. Aborting...");
            return None;
        }
    }

    None
}

/// QC on ntrip real-time data
pub async fn ntripqc(
    host: &str,
    port: &str,
    mountpoint: &str,
    username: &str,
    password: &str,
    supple_nav: Option<Nav>,
    mut outtype: Outtype
) {
    let connect= conntrip(host, port, mountpoint, username, password).await;
    if connect.is_none() { return; }
    let mut stream=connect.unwrap();

    let mut navs = Nav::new();
    let mut obss = Obss::new();
    let mut epochtime = 0;
    let mut qc = QC::init();
    let popt = PrcOpt::default();
    let mut rtk = Rtk::init(&popt);

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

    if let Outtype::File(path) = outtype {
        let qcfile = create_qcfile(mountpoint, &path);
        outtype = Outtype::File(qcfile);
    }
    let decoder = Decoder::new();
    unsafe {
        (decoder.free_rtcm)(decoder.rtcm_ptr);
        libc::free(decoder.rtcm_ptr as *mut libc::c_void);
    }
    loop {
        tokio::select! {
            _ = async {
                if bjt_0(){
                    qc=QC::init();
                }
                buf.fill(0);

                stream.read(&mut buf).await.expect("error in stream read");
                let rtcm3 = get_single_rtcm(&buf);
                for data in rtcm3{
                    rt_rtcm2qc(
                        &decoder,
                        mountpoint,
                        data,
                        &mut epochtime,
                        &mut obss,
                        &mut navs,
                        &popt,
                        &mut rtk,
                        &mut qc,
                        &supple_nav,
                        outtype.clone()
                    )
                }
            } => {},

            // Processing user input
            Some(input) = rx.recv() => {
                match input.trim() {
                    // "qc" => {
                    //     println!("QC processing...");
                    //     qc_cal_rt(&mut qc,popt.navsys);
                    //     let navflag=if navs.n==0||navs.ng==0 {false} else {true};
                    //         if let Outtype::File(path) = outtype {
                    //             let qcfile_rt=create_realtime_qcfile(mountpoint,path);
                    //             outtype = &Outtype::File(qcfile_rt);
                    //         }
                    //
                    //     if !qcfile_rt.is_empty(){
                    //         let file=File::create(qcfile_rt).expect("Fail to create output file!");
                    //         let mut writer = BufWriter::new(file);
                    //         qc_out_rt(&qc,popt.navsys,navflag,&rtk,&mut writer);
                    //     }
                    //
                    // },
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


