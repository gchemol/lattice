// [[file:../lattice.note::4d507fb4][4d507fb4]]
use gchemol_lattice::Lattice;

use approx::*;
use vecfx::*;

#[test]
fn test_lattice_construct() {
    let mut lat = Lattice::default();
    let loc = [1.0, 2.0, 3.0];
    lat.set_origin(loc);

    let lat = Lattice::new([[18.256, 0., 0.], [0., 20.534, 0.], [0., 0., 15.084]]);
    assert_eq!(true, lat.is_orthorhombic());

    let lat = Lattice::new([[15.3643, 0., 0.], [4.5807, 15.5026, 0.], [0., 0., 17.4858]]);
    let [a, b, c] = lat.lengths();
    assert_eq!(false, lat.is_orthorhombic());
    assert_relative_eq!(a, 15.3643, epsilon = 1e-4);
    assert_relative_eq!(b, 16.1652, epsilon = 1e-4);
    assert_relative_eq!(c, 17.4858, epsilon = 1e-4);

    let [alpha, beta, gamma] = lat.angles();
    assert_relative_eq!(alpha, 90.0, epsilon = 1e-4);
    assert_relative_eq!(beta, 90.0, epsilon = 1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon = 1e-4);

    let lat = Lattice::from_params(a, b, c, alpha, beta, gamma);
    assert_eq!([a, b, c], lat.lengths());
    assert_eq!([alpha, beta, gamma], lat.angles());

    // scale_by_a, scale_by_b, scale_by_c
    {
        let mut l1 = lat.clone();
        let [a, b, c] = l1.lengths();
        l1.scale_by_a(2.0);
        let [a2, b, c] = l1.lengths();

        assert_eq!(a2, a * 2.0);
        let mut l1 = lat.clone();
        l1.scale_by_b(2.0);
        let [_, b2, _] = l1.lengths();
        assert_eq!(b2, b * 2.0);

        let mut l1 = lat.clone();
        l1.scale_by_c(2.0);
        let [_, _, c2] = l1.lengths();
        assert_eq!(c2, c * 2.0);
    }

    let vectors = lat.vectors();
    let _ = Lattice::new(vectors);
    let mat = lat.matrix();
    let _ = Lattice::from_matrix(mat);
}

#[test]
fn test_lattice_volume() {
    let vts = [[5., 0., 0.], [5., 5., 0.], [1., 0., 5.]];
    let mut lat = Lattice::new(vts);
    let mat: Matrix3f = vts.into();
    assert_eq!(mat, lat.matrix());
    let va: [f64; 3] = lat.vector_a().into();
    let vb: [f64; 3] = lat.vector_b().into();
    let vc: [f64; 3] = lat.vector_c().into();
    assert_eq!(vts[0], va);
    assert_eq!(vts[1], vb);
    assert_eq!(vts[2], vc);

    assert_relative_eq!(125.0, lat.volume(), epsilon = 1e-4);
    lat.scale_by(4.);
    assert_relative_eq!(8000.0, lat.volume(), epsilon = 1e-4);
}

#[test]
fn test_lattice_frac_cart() {
    // ovito/tests/files/LAMMPS/multi_sequence_1.dump
    let lat = Lattice::new([[5.09, 0.00, 0.00], [0.00, 6.74, 0.00], [0.00, 0.00, 4.53]]);

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
// adopted from lumol
fn test_wrap() {
    // Cubic unit cell
    let cell = Lattice::from_params(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    let wrapped: Vector3f = cell.wrap([9.0, 18.0, -6.0]).into();
    assert_relative_eq!(wrapped, Vector3f::from([9.0, 8.0, 4.0]), epsilon = 1e-4);

    // Orthorhombic unit cell
    let cell = Lattice::from_params(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
    let wrapped = cell.wrap([1.0, 1.5, 6.0]);
    assert_relative_eq!(wrapped, Vector3f::from([1.0, 1.5, 1.0]), epsilon = 1e-4);

    // Triclinic unit cell
    let cell = Lattice::from_params(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
    let wrapped = cell.wrap([1.0, 1.5, 6.0]);
    assert_relative_eq!(wrapped, Vector3f::from([1.0, 1.5, 1.0]), epsilon = 1e-4);
}
// 4d507fb4 ends here
