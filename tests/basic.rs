// basic.rs
// :PROPERTIES:
// :header-args: :tangle tests/basic.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*basic.rs][basic.rs:1]]
use gchemol_lattice::Lattice;

use approx::*;
use vecfx::*;

#[test]
fn test_lattice_construct() {
    let mut lat = Lattice::default();
    let loc = [1.0, 2.0, 3.0];
    lat.set_origin(loc);
    assert_eq!(loc, lat.origin());

    let mut lat = Lattice::new([[15.3643, 0., 0.], [4.5807, 15.5026, 0.], [0., 0., 17.4858]]);

    let [a, b, c] = lat.lengths();
    assert_eq!(false, lat.is_orthorhombic());

    assert_relative_eq!(a, 15.3643, epsilon = 1e-4);
    assert_relative_eq!(b, 16.1652, epsilon = 1e-4);
    assert_relative_eq!(c, 17.4858, epsilon = 1e-4);

    let [alpha, beta, gamma] = lat.angles();
    assert_relative_eq!(alpha, 90.0, epsilon = 1e-4);
    assert_relative_eq!(beta, 90.0, epsilon = 1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon = 1e-4);

    let mut lat = Lattice::from_params(a, b, c, alpha, beta, gamma);
    assert_eq!([a, b, c], lat.lengths());
    assert_eq!([alpha, beta, gamma], lat.angles());
}

#[test]
fn test_lattice_volume() {
    let vts = [[5., 0., 0.], [5., 5., 0.], [1., 0., 5.]];

    let mut lat = Lattice::new(vts);
    assert_eq!(vts, lat.vectors());
    assert_eq!(vts[0], lat.vector_a());
    assert_eq!(vts[1], lat.vector_b());
    assert_eq!(vts[2], lat.vector_c());

    assert_relative_eq!(125.0, lat.volume(), epsilon = 1e-4);
    lat.scale_by(4.);
    assert_relative_eq!(8000.0, lat.volume(), epsilon = 1e-4);
}

#[test]
fn test_lattice_frac_cart() {
    // ovito/tests/files/LAMMPS/multi_sequence_1.dump
    let mut lat = Lattice::new([[5.09, 0.00, 0.00], [0.00, 6.74, 0.00], [0.00, 0.00, 4.53]]);

    let fs = lat.to_frac([2.1832, 1.6850, 3.8505]);
    assert_relative_eq!(fs[0], 0.4289, epsilon = 1e-3);
    assert_relative_eq!(fs[1], 0.2500, epsilon = 1e-3);
    assert_relative_eq!(fs[2], 0.8500, epsilon = 1e-3);
    let fs = lat.to_frac([6.9068, 5.0550, 0.6795]);
    assert_relative_eq!(fs[0], 1.3569, epsilon = 1e-3);
    let fs = lat.to_frac([4.3618, 5.0550, 1.5855]);
    assert_relative_eq!(fs[2], 0.3500, epsilon = 1e-3);

    let coords = lat.to_cart([0.4289, 0.2500, 0.8500]);
    assert_relative_eq!(coords[0], 2.1832, epsilon = 1e-3);
    assert_relative_eq!(coords[1], 1.6850, epsilon = 1e-3);
    assert_relative_eq!(coords[2], 3.8505, epsilon = 1e-3);
}

#[test]
fn test_mic_distance() {
    // Setup lattice
    // a = b = c = 4, alpha = beta = gamma = 60
    let cell = [
        [4.00000000, 0.00000000, 0.00000000],
        [2.00000000, 3.46410162, 0.00000000],
        [2.00000000, 1.15470054, 3.26598632],
    ];
    let mut lattice = Lattice::new(cell);

    // Safe distance range where Tuckermann algorithm will work
    let safe_r_max = 0.5 * lattice.widths().min();
    assert_relative_eq!(safe_r_max, 1.4142, epsilon = 1e-4);
    let pi = [-0.0000000, 0.0000000, -0.0000000];
    let pj = [-0.0743502, 2.5356374, -2.0623249];
    let dij = lattice.distance(pi, pj);
    assert_relative_eq!(dij, 1.2270, epsilon = 1e-4);

    // When true mic distance out of the safe range for Tuckermann algorithm.
    let pj = [-0.0834941, 1.8252187, -1.5169388];
    let dij = lattice.distance(pi, pj);
    assert_relative_eq!(dij, 1.8167, epsilon = 1e-4);
}

#[test]
fn test_mic_vector() {
    let mut lat = Lattice::new([
        [7.055, 0., 0.],
        [0., 6.795, 0.],
        [-1.14679575, 0., 5.65182701],
    ]);

    // mic vector
    let expected = Vector3f::from([-0.48651737, 0.184824, -1.31913642]);
    let pmic = lat.apply_mic([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic, epsilon = 1e-4);

    let pmic = lat.apply_mic([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic, epsilon = 1e-4);
}
// basic.rs:1 ends here
