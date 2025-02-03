use zmq::{Context, Message};
use std::io;
use std::io::BufRead;
use std::time::Duration;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use redis::{Client, Commands, Connection, RedisResult};
use serde_json::json;
use tokio::sync::mpsc;
use once_cell::sync::Lazy;
use crate::basic::sat::satno2id;
use crate::basic::time::{bjt_0, time2str};
use crate::basic::var::{Nav, NFREOBS, NUMSYS, Obs, OBS_CODES, Obss, PrcOpt, Rtk};
use crate::convbin::decode::Decoder;
use crate::qc::qc_fun::{codefreq, namebysyslog, rt_rtcm2qc};
use crate::qc::qc_var::{Outtype, QC};

// mod dreamnet {
//     include!(concat!(env!("OUT_DIR"), "/protobuf_dreamnetmsg.rs"));
// }

struct StationData {
    qc: QC,
    epochtime: i64,
    obss: Obss,
    navs: Nav,
    rtk: Rtk,
}

impl StationData {
    fn new(popt: &PrcOpt) -> Self {
        StationData {
            qc: QC::init(),
            epochtime: 0,
            obss: Obss::new(),
            navs: Nav::new(),
            rtk: Rtk::init(popt),
        }
    }
}

static REDIS_CONNECTION: Lazy<Mutex<Option<Arc<Mutex<Connection>>>>> = Lazy::new(|| Mutex::new(None));

/// Creating a Redis Connection
fn create_redis_connection(url: String) -> RedisResult<Connection> {
    let client = Client::open(url)?;

    client.get_connection_with_timeout(Duration::from_secs(5))
}

/// Getting a global Redis connection
fn get_redis_connection(url: String) -> RedisResult<Arc<Mutex<Connection>>> {
    let mut connection = REDIS_CONNECTION.lock().unwrap();

    if connection.is_none() {
        match create_redis_connection(url) {
            Ok(conn) => {
                *connection = Some(Arc::new(Mutex::new(conn)));
                println!("redis connection success")
            }
            Err(e) => {
                return Err(e);
            }
        }
    }

    Ok(Arc::clone(connection.as_ref().unwrap()))
}


pub async fn con_zmq(
    ip: &str,
    port: &str,
    stations: Option<Vec<String>>,
    supple_nav: Option<Nav>,
    outtype: Outtype,
) {
    if ip.is_empty() || port.is_empty() {
        eprintln!("Empty parameter");
        return;
    }
    let url = format!("tcp://{}:{}", ip, port);
    let context = Context::new();
    let subscriber = context.socket(zmq::SUB).expect("Failed to create socket");
    subscriber.connect(&url).expect("Failed in connect");
    println!("zmq for {}", url);
    subscriber.set_rcvtimeo(5000).expect("Failed to set receive timeout");

    if stations.is_none() {
        // Subscribe to all threads
        subscriber.set_subscribe(b"").expect("Failed to subscribe to all topics");
    } else {
        // Subscribe to the specified topics in turn
        for topic in stations.unwrap() {
            subscriber.set_subscribe(topic.as_bytes()).expect(&format!("Failed to subscribe to the topic:{}", topic));
        }
    }

    // Calling the decode library
    let decoder = Decoder::new();
    unsafe {
        (decoder.free_rtcm)(decoder.rtcm_ptr);
        libc::free(decoder.rtcm_ptr as *mut libc::c_void);
    }

    // Define a HashMap for storing site data
    let mut station_data_map: HashMap<String, StationData> = HashMap::new();

    let popt = PrcOpt::default();

    let (tx, mut rx) = mpsc::channel::<String>(100);
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

    loop {
        tokio::select! {
            _ = async {
                // timer::start_timing();
                // Receive theme frames
                let mut topic = Message::new();
                match subscriber.recv(&mut topic, 0) {
                    Ok(_) => {
                        // Theme received successfully
                        let topic_str = topic.as_str().unwrap_or("empty topic");
                        // println!("Received topic: {}", topic_str);

                        // Generate new output types based on themes
                        let new_outtype = if let Outtype::File(ref path) = outtype {
                            let file = format!("{}/{}", path, topic_str);
                            Outtype::File(file)
                        } else if let Outtype::All(ref path, ref r) = outtype {
                            let file = format!("{}/{}", path, topic_str);
                            Outtype::All(file, r.clone())
                        } else {
                            outtype.clone()  // If `outtype` is neither `File` nor `All`, it remains unchanged.
                        };

                        // Receive data frame
                        let mut msg = Message::new();
                        match subscriber.recv(&mut msg, 0) {
                            Ok(_) => {
                                // Data received successfully
                                let rtcm_bin: Vec<u8> = msg.to_vec();

                                // Getting or initializing data for the site
                                let station_data = station_data_map.entry(topic_str.to_string())
                                    .or_insert_with(|| StationData::new(&popt));

                                if bjt_0() {
                                    station_data.qc = QC::init();
                                }

                                // qc progress
                                rt_rtcm2qc(
                                    &decoder,
                                    topic_str,
                                    rtcm_bin,
                                    &mut station_data.epochtime,
                                    &mut station_data.obss,
                                    &mut station_data.navs,
                                    &popt,
                                    &mut station_data.rtk,
                                    &mut station_data.qc,
                                    &supple_nav,
                                    new_outtype,
                                );
                            }
                            Err(zmq::Error::EAGAIN) => {
                                println!("No data received for topic: {}, retrying...", topic_str);
                            }
                            Err(e) => {
                                eprintln!("Error receiving data for topic {}: {}", topic_str, e);
                            }
                        }
                    }
                    Err(zmq::Error::EAGAIN) => {
                        println!("No topic received, retrying...");
                    }
                    Err(e) => {
                        eprintln!("Error receiving topic: {}", e);
                    }
                }
            } => {},

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
            }
        }
    }
}

// output rt_qc result to redis
pub fn qc2redis(
    site: &str,
    qc: &QC,
    obss: &Obss,
    redis_url: Option<String>,
    redis_port: Option<String>,
    redis_auth: Option<String>,
) -> RedisResult<()> {
    // connect redis
    let url = format!("redis://:{}@{}:{}/", redis_auth.unwrap_or_default(), redis_url.unwrap_or_default(), redis_port.unwrap_or_default());
    let con = get_redis_connection(url)?;
    let mut con = con.lock().unwrap();
    // customize the field structure for output to redis
    let key = format!("StationQC:Station_{}", site);
    let mut station_data = HashMap::new();
    // Adding General Fields
    station_data.insert("fstep", json!(time2str(qc.totinfo.firstepoch)));
    station_data.insert("lstep", json!(time2str(qc.totinfo.lastepoch)));
    station_data.insert("sample", json!(qc.totinfo.sample));
    station_data.insert("eprt", json!(qc.totinfo.ep_rt));
    station_data.insert("obsrt", json!(qc.totinfo.obs_rt));
    station_data.insert("obsrt10", json!(qc.totinfo.obs_el_rt));
    station_data.insert("oslps", json!(qc.totinfo.oslps));

    let mut mp_data = HashMap::new();
    for i in 0..NUMSYS {
        if qc.sysinfo[i].codenum == 0 { continue; }
        let sysname = namebysyslog(i);
        let mut sys_mp = HashMap::new();
        for j in 0..NFREOBS {
            if qc.sysinfo[i].totmp[j] == 0.0 { continue; }
            sys_mp.insert(format!("mp{}", codefreq(qc.sysinfo[i].code[j])), qc.sysinfo[i].totmp[j]);
        }
        // Storage by System
        mp_data.insert(sysname, sys_mp);
    }
    // Convert MP nested tables to json deposit
    station_data.insert("MP", json!(mp_data));

    let mut snr_list = Vec::new();
    for i in 0..obss.obs.len() {
        let obs: Obs = obss.obs[i].clone();
        let id = satno2id(obs.sat);
        for j in 0..NFREOBS {
            let snr = obs.snr[j] / 1000;
            if snr == 0 { continue; }
            let signame = OBS_CODES[obs.code[j] as usize];

            // Constructing Nested SNR Information for Individual Satellites
            let mut snr_info = HashMap::new();
            snr_info.insert("code", signame);
            let snr_str = snr.to_string();
            snr_info.insert("snr", &snr_str);
            // Serialize snr_info to JSON string
            let snr_info_json = serde_json::to_string(&snr_info).expect("Failed to serialize snr_info");
            // Store prn and its corresponding snr
            let mut prn_info = HashMap::new();
            prn_info.insert("prn", id.clone());
            prn_info.insert("snr", snr_info_json);

            snr_list.push(prn_info)
        }
    }
    // Convert SNR list to json deposit
    station_data.insert("SNR", json!(snr_list));

    // Serialize to JSON and store in Redis
    let station_data_json = serde_json::to_string(&station_data).expect("Failed to serialize station data");
    let _: () = con.set(&key, station_data_json)?;

    Ok(())
}
