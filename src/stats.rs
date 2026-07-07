/// z-score const for a two-sided 95% confint
pub const WILSON_Z_95: f64 = 1.959963984540054;

/// Wilson score confidence interval for a binomial proportion
pub fn wilson_interval(successes: usize, n: usize, z: f64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 0.0);
    }
    let n = n as f64;
    let p = successes as f64 / n;
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denom;
    let margin = (z / denom) * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt();
    ((center - margin).max(0.0), (center + margin).min(1.0))
}

/// Standard normal survival function, P(Z >= x)
pub fn normal_survival(x: f64) -> f64 {
    (1.0 - normal_cdf(x)).clamp(0.0, 1.0)
}

fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    // Abramowitz and Stegun 7.1.26 approximation
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + P * x);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-x * x).exp();
    sign * y
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-4, "expected {b}, got {a}");
    }

    #[test]
    fn test_wilson_reference_values() {
        // Known ref: 8 out of 10 at 95% -> ~(0.4901, 0.9433)
        let (lo, hi) = wilson_interval(8, 10, WILSON_Z_95);
        approx(lo, 0.490_1);
        approx(hi, 0.943_3);
    }

    #[test]
    fn test_wilson_contains_point_estimate() {
        for &(s, n) in &[(1usize, 5usize), (50, 100), (999, 1000), (3, 7)] {
            let p = s as f64 / n as f64;
            let (lo, hi) = wilson_interval(s, n, WILSON_Z_95);
            assert!(lo <= p && p <= hi, "p={p} not in [{lo}, {hi}]");
            assert!((0.0..=1.0).contains(&lo) && (0.0..=1.0).contains(&hi));
        }
    }

    #[test]
    fn test_wilson_edges() {
        // p = 0: lower bd pinned at 0, high bd positive
        let (lo, hi) = wilson_interval(0, 20, WILSON_Z_95);
        approx(lo, 0.0);
        assert!(hi > 0.0 && hi < 1.0);

        // p = 1: high bound pinned at 1, low bound below 1
        let (lo, hi) = wilson_interval(20, 20, WILSON_Z_95);
        approx(hi, 1.0);
        assert!(lo < 1.0 && lo > 0.0);

        // n = 0: degenerate
        assert_eq!(wilson_interval(0, 0, WILSON_Z_95), (0.0, 0.0));
    }

    #[test]
    fn test_wilson_widens_as_n_shrinks() {
        // Same prop (0.9), fewer trials, wider interval
        let (lo_small, hi_small) = wilson_interval(9, 10, WILSON_Z_95);
        let (lo_large, hi_large) = wilson_interval(900, 1000, WILSON_Z_95);
        assert!(hi_small - lo_small > hi_large - lo_large);
    }

    #[test]
    fn test_normal_survival_reference_values() {
        approx(normal_survival(0.0), 0.5);
        approx(normal_survival(1.959_963_984_540_054), 0.025);
    }
}
