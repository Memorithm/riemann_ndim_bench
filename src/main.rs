use riemann_ndim_bench::{PiRadialGeometry, SpectralPoint};

fn main() {
    println!("Riemann N-dimensional bench - foundation v0");
    println!("This stage checks symmetry and radial-coordinate identities only.\n");

    println!("sigma      mirror      raw_radius       normalized       log_normalized");
    println!("--------------------------------------------------------------------------");

    for sigma in [0.1_f64, 0.3, 0.49, 0.5, 0.51, 0.7, 0.9] {
        let point = SpectralPoint::new(sigma, 0.0);
        let mirror = point.critical_line_reflection();
        let raw = PiRadialGeometry::raw_radius(sigma);
        let normalized = PiRadialGeometry::normalized_radius(sigma);
        let log_normalized = PiRadialGeometry::log_normalized_radius(sigma);

        println!(
            "{sigma:<10.3} {mirror_sigma:<10.3} {raw:<16.10} {normalized:<16.10} {log_normalized:<+16.10}",
            mirror_sigma = mirror.sigma,
        );
    }

    println!(
        "\ncritical raw radius = {:.15}",
        PiRadialGeometry::critical_radius()
    );
}
