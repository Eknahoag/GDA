
use std::time::{Instant, Duration};
use once_cell::sync::Lazy;
use std::sync::Mutex;

pub static START_TIME: Lazy<Mutex<Option<Instant>>> = Lazy::new(|| Mutex::new(None));
pub static TOTAL_DURATION: Lazy<Mutex<Duration>> = Lazy::new(|| Mutex::new(Duration::new(0, 0)));
pub static COUNT: Lazy<Mutex<u64>> = Lazy::new(|| Mutex::new(0));

pub fn start_timing() {
    let mut start_time = START_TIME.lock().unwrap();
    *start_time = Some(Instant::now());
}

pub fn stop_timing() {
    let mut count = COUNT.lock().unwrap();
    let mut total_duration = TOTAL_DURATION.lock().unwrap();
    let start_time = START_TIME.lock().unwrap();

    if let Some(start) = *start_time {
        *total_duration += start.elapsed();
        *count += 1;
    }
}

pub fn get_average_time() -> Option<Duration> {
    let count = COUNT.lock().unwrap();
    let total_duration = TOTAL_DURATION.lock().unwrap();

    if *count > 0 {
        Some(*total_duration / *count as u32)
    } else {
        None
    }
}
