use std::cmp::Ordering;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Lines};
use std::ops::Mul;
use std::path::Path;
use num_traits::Float;

pub fn sqr<T>(x: T) -> T
where
    T: Mul<Output = T> + Copy,
{
    x * x
}

pub fn norm<T>(vector: &[T]) -> T
where
    T: Float,
{
    vector.iter()
        .map(|&x| x * x)
        .fold(T::zero(), |acc, x| acc + x)
        .sqrt()
}

pub trait Sqrt {
    fn sqrt(self) -> Self;
}

impl Sqrt for f64 {
    fn sqrt(self) -> f64 {
        f64::sqrt(self)
    }
}

impl Sqrt for f32 {
    fn sqrt(self) -> f32 {
        f32::sqrt(self)
    }
}

impl Sqrt for usize {
    fn sqrt(self) -> usize {
        f32::sqrt(self as f32) as usize
    }
}
pub fn sqvar<T>(covar: T) -> T
where
    T: Mul<Output = T> + PartialOrd + Copy + From<f64> + Into<f64> + Sqrt,
{
    if covar < T::from(0.0) {
        T::from((-covar.into()).sqrt())
    } else {
        covar.sqrt()
    }
}

pub fn openfile<P: AsRef<Path>>(path: P) -> io::Result<Lines<BufReader<File>>> {
    let file = File::open(path)?;
    Ok(BufReader::new(file).lines())
}

pub fn rms(data: &[f32]) -> f32 {
    let sum_of_squares: f32 = data.iter().map(|&x| x * x).sum();
    let mean_of_squares = sum_of_squares / data.len() as f32;
    mean_of_squares.sqrt()
}

// Calculate the five-number summary of a dataset (minimum, first quartile, median, third quartile, maximum)
// Returns None if the input data is empty
pub fn cal_box(data: &mut Vec<f64>) -> Option<[f64; 5]> {
    // Check if the input data is empty
    if data.is_empty() {
        return None;
    }

    // Sorting data
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    // Define the function to calculate the percentile
    let percentile = |data: &Vec<f64>, p: f64| -> f64 {
        let n = data.len() as f64;
        let pos = p * (n - 1.0);
        let low = pos.floor() as usize;
        let high = pos.ceil() as usize;
        if low == high {
            data[low]
        } else {
            let low_val = data[low];
            let high_val = data[high];
            low_val + (high_val - low_val) * (pos - low as f64)
        }
    };

    // Special case handling: if the length of data is less than 4
    let n = data.len();
    if n == 1 {
        let single_value = data[0];
        return Some([single_value, single_value, single_value, single_value, single_value]);
    } else if n == 2 {
        let min_val = data[0];
        let max_val = data[1];
        let median = (min_val + max_val) / 2.0;
        return Some([min_val, median, median, median, max_val]);
    } else if n == 3 {
        let min_val = data[0];
        let median = data[1];
        let max_val = data[2];
        return Some([min_val, min_val, median, max_val, max_val]);
    }

    // Calculate first quartile (Q1), median (Q2), third quartile (Q3)
    let q1 = percentile(&data, 0.25);
    let q3 = percentile(&data, 0.75);
    let iqr = q3 - q1; // Interquartile range
    let lower_bound = q1 - 1.5 * iqr; // Lower bound for outliers
    let upper_bound = q3 + 1.5 * iqr; // Upper bound for outliers

    // Filter out the outliers
    let filtered_data: Vec<f64> = data.iter()
        .filter(|&&x| x >= lower_bound && x <= upper_bound)
        .cloned()
        .collect();

    // Check if there are any non-outlier values
    if filtered_data.is_empty() {
        return None; // If all values are outliers
    }

    // Calculate new minimum and maximum values from filtered data
    let min_val = *filtered_data.first().unwrap();
    let max_val = *filtered_data.last().unwrap();
    let median = percentile(&filtered_data, 0.5); // Median from filtered data

    Some([q1, median, q3, min_val, max_val])
}