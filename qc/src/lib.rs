pub mod basic;
pub mod qc;
pub mod convbin;
pub mod timer;

#[cfg(test)]
mod tests {
    use std::path::Path;
    use std::fs::File;
    use std::io::{BufWriter, Write};
    use tokio::time::{Duration, timeout};
    use tokio::io::AsyncReadExt;
    use crate::basic::ntrip::{get_single_rtcm, conntrip};
    use crate::basic::time::date2gtime;
    use crate::basic::var::{Nav, Obss, PrcOpt, QcOpt, Rtk, SYS_CMP, SYS_GAL, SYS_GLO, SYS_GPS, SYS_QZS};
    use crate::convbin::decode::Decoder;
    use crate::convbin::rnxout::rtcm2rnx;
    use crate::qc::qc_fun::{obsqc, rt_rtcm2qc, rtcmqc};
    use crate::qc::qc_var::{Outtype, QC};

    #[test]
    fn test_decode() {
        let ifile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.2023090binRTCM3";
        let obsfile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.obs".to_string();
        let navfile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.nav".to_string();
        let mut opt = QcOpt::default();
        opt.tr = date2gtime("2023/03/31 00:00:00");
        opt.navsys = SYS_GPS | SYS_GLO | SYS_GAL | SYS_CMP | SYS_QZS;
        opt.rnxver = 304;
        let _ = rtcm2rnx(ifile, Some(obsfile), Some(navfile), opt);
    }

    #[tokio::test]
    async fn test_conntrip() {
        let host = "ntrip-server.navfirst.com";
        let port = "9095";
        let mountpoint = "JACK03";
        let username = "JACK03";
        let password = "JACK03";
        let result = timeout(Duration::from_secs(60), async {
            let connect = conntrip(host, port, mountpoint, username, password).await;
            if connect.is_none() { return; }
        }).await;
        if result.is_err() { panic!() }
    }

    #[tokio::test]
    async fn test_ntripqc() {
        let host = "ntrip-server.navfirst.com";
        let port = "9095";
        let mountpoint = "JACK02";
        let username = "JACK02";
        let password = "JACK02";
        let file_path = "D:\\code\\rust\\qc\\sample\\JACK02.qc";

        let outtype = Outtype::File(file_path.to_string());
        let popt = PrcOpt::default();
        let result = timeout(Duration::from_secs(20), async {
            let connect = conntrip(host, port, mountpoint, username, password).await;
            if connect.is_none() { return; }
            let mut stream = connect.unwrap();
            let mut buf = vec![0; 2048];
            let decoder = Decoder::new();
            unsafe {
                (decoder.free_rtcm)(decoder.rtcm_ptr);
                libc::free(decoder.rtcm_ptr as *mut libc::c_void);
            }

            loop {
                stream.read(&mut buf).await.expect("error in stream read");
                let rtcm3 = get_single_rtcm(&buf);
                for data in rtcm3 {
                    rt_rtcm2qc(
                        &decoder,
                        mountpoint,
                        data,
                        &mut 0,
                        &mut Obss::new(),
                        &mut Nav::new(),
                        &PrcOpt::default(),
                        &mut Rtk::init(&popt),
                        &mut QC::init(),
                        &None,
                        outtype.clone(),
                    )
                }

                if Path::new(file_path).exists() {
                    let metadata = std::fs::metadata(file_path).unwrap();
                    if metadata.len() > 0 {
                        return;
                    }
                }
            }
        }).await;

        if result.is_err() { panic!() }
    }

    #[tokio::test]
    async fn test_binstore() {
        let host = "ntrip-server.navfirst.com";
        let port = "9095";
        let mountpoint = "JACK02";
        let username = "JACK02";
        let password = "JACK02";
        let outfile = "D:\\code\\rust\\qc\\sample\\JACK02";

        let file = File::create(&outfile).expect("Fail to create output obs file!");
        let mut writer = BufWriter::new(file);
        let result = timeout(Duration::from_secs(20), async {
            let connect = conntrip(host, port, mountpoint, username, password).await;
            if connect.is_none() { return; }
            let mut stream = connect.unwrap();
            let mut buf = vec![0; 2048];

            loop {
                let n = stream.read(&mut buf).await.expect("error in stream read");
                if n > 0
                {
                    writer.write_all(&buf[..n]).expect("Fail to store rtcm");
                    writer.flush().expect("Failed to flush writer");
                    writer.get_ref().sync_all().expect("Failed to sync data");
                    return;
                }

                if Path::new(outfile).exists() {
                    let metadata = std::fs::metadata(outfile).unwrap();
                    if metadata.len() > 0 {
                        return;
                    }
                }
            }
        }).await;

        if result.is_err() { panic!() }
    }


    #[test]
    fn test_obsqc() {
        let mut ifile = Vec::new();
        ifile.push("D:\\code\\rust\\qc\\sample\\PXCEDL36_test.obs".to_string());
        ifile.push("D:\\code\\rust\\qc\\sample\\PXCEDL36_test.nav".to_string());
        let ofile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.qc";
        let mut opt = QcOpt::default();
        opt.detail = 2;
        obsqc(&mut opt, &ifile, ofile)
    }

    #[test]
    fn test_rtcmqc() {
        let ifile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.2023090binRTCM3";
        let ofile = "D:\\code\\rust\\qc\\sample\\PXCEDL36.qc";
        let mut opt = QcOpt::default();
        opt.tr = date2gtime("2023/03/31 00:00:00");
        opt.detail = 2;
        rtcmqc(&mut opt, ifile, ofile)
    }
}
 