//! xoshiro256++ with splitmix64 seeding; Box-Muller normals. Deterministic across platforms.
pub struct Rng { s: [u64; 4], spare: Option<f64> }
impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut x = seed; let mut s = [0u64; 4];
        for v in s.iter_mut() { x = x.wrapping_add(0x9E3779B97F4A7C15); let mut z = x; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); *v = z ^ (z >> 31) }
        Rng { s, spare: None }
    }
    pub fn next_u64(&mut self) -> u64 {
        let r = (self.s[0].wrapping_add(self.s[3])).rotate_left(23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1]; self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].rotate_left(45);
        r
    }
    pub fn uniform(&mut self) -> f64 { (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) }
    pub fn below(&mut self, n: usize) -> usize { (self.uniform() * n as f64) as usize % n.max(1) }
    pub fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() { return z }
        let (mut u, v) = (self.uniform(), self.uniform());
        if u < 1e-300 { u = 1e-300 }
        let r = (-2.0 * u.ln()).sqrt(); let th = 2.0 * std::f64::consts::PI * v;
        self.spare = Some(r * th.sin());
        r * th.cos()
    }
}
