//! CPU topology (Linux /sys) for thread pinning: one logical CPU per physical core, sockets interleaved,
//! so memory-bound loci spread over both memory controllers and each locus' LD is first-touched locally.
pub fn physical_cores() -> Vec<usize> {
    let mut seen = std::collections::BTreeMap::<(usize, usize), usize>::new();
    let mut i = 0;
    loop {
        let base = format!("/sys/devices/system/cpu/cpu{i}/topology");
        let (Ok(pk), Ok(co)) = (std::fs::read_to_string(format!("{base}/physical_package_id")), std::fs::read_to_string(format!("{base}/core_id"))) else { break };
        let key = (pk.trim().parse().unwrap_or(0), co.trim().parse().unwrap_or(0));
        seen.entry(key).or_insert(i);
        i += 1;
    }
    if seen.is_empty() { return (0..std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1)).collect() }
    let mut by_pkg: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for ((pk, _), cpu) in seen { by_pkg.entry(pk).or_default().push(cpu) }
    let lists: Vec<Vec<usize>> = by_pkg.into_values().collect();
    let mut out = vec![]; let mx = lists.iter().map(|l| l.len()).max().unwrap_or(0);
    for k in 0..mx { for l in &lists { if k < l.len() { out.push(l[k]) } } }
    out
}

pub fn pool(threads: usize, pin: bool) -> rayon::ThreadPool {
    let cores = physical_cores();
    rayon::ThreadPoolBuilder::new().num_threads(threads).start_handler(move |i| {
        if pin && !cores.is_empty() { core_affinity::set_for_current(core_affinity::CoreId { id: cores[i % cores.len()] }); }
    }).build().unwrap()
}
