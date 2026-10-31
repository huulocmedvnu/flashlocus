//! Port of R's Brent_fmin (src/library/stats/src/optimize.c), used by optim(method = "Brent").

pub fn fmin<F: FnMut(f64) -> f64>(ax: f64, bx: f64, f: &mut F, tol: f64) -> f64 {
    let c = (3.0 - 5f64.sqrt()) * 0.5;
    let eps = f64::EPSILON.sqrt();
    let (mut a, mut b) = (ax, bx);
    let mut v = a + c * (b - a);
    let (mut w, mut x) = (v, v);
    let (mut d, mut e) = (0.0f64, 0.0f64);
    let mut fx = f(x);
    let (mut fv, mut fw) = (fx, fx);
    let tol3 = tol / 3.0;
    loop {
        let xm = (a + b) * 0.5;
        let tol1 = eps * x.abs() + tol3;
        let t2 = tol1 * 2.0;
        if (x - xm).abs() <= t2 - (b - a) * 0.5 {
            break;
        }
        let (mut p, mut q, mut r) = (0.0f64, 0.0f64, 0.0f64);
        if e.abs() > tol1 {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.0;
            if q > 0.0 { p = -p } else { q = -q }
            r = e;
            e = d;
        }
        if p.abs() >= (q * 0.5 * r).abs() || p <= q * (a - x) || p >= q * (b - x) {
            e = if x < xm { b - x } else { a - x };
            d = c * e;
        } else {
            d = p / q;
            let u = x + d;
            if u - a < t2 || b - u < t2 {
                d = tol1;
                if x >= xm { d = -d }
            }
        }
        let u = if d.abs() >= tol1 { x + d } else if d > 0.0 { x + tol1 } else { x - tol1 };
        let fu = f(u);
        if fu <= fx {
            if u < x { b = x } else { a = x }
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if u < x { a = u } else { b = u }
            if fu <= fw || w == x {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u; fv = fu;
            }
        }
    }
    x
}
