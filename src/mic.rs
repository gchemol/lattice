// imports

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*imports][imports:1]]
use vecfx::*;

#[cfg(test)]
use approx::*;

use crate::Lattice;
// imports:1 ends here

// Tuckerman algorithm

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*Tuckerman algorithm][Tuckerman algorithm:1]]
impl Lattice {
    /// Return the approximated mic vector using Tuckerman's algorithm.
    ///
    /// Reference
    /// ---------
    /// - Tuckerman, M. E. Statistical Mechanics: Theory and Molecular
    /// Simulation, 1 edition.; Oxford University Press: Oxford ; New York,
    /// 2010.
    pub(crate) fn apply_mic_tuckerman(&mut self, p: [f64; 3]) -> Vector3f {
        // apply minimum image convention on the scaled coordinates
        let mut fcoords = self.to_frac(p);

        let mut image = [1.0; 3];
        for i in 0..3 {
            image[i] = -1.0 * fcoords[i].round();
            fcoords[i] += image[i];
        }

        // transform back to cartesian coordinates
        let pij = self.to_cart(fcoords);

        Vector3f::from(pij)
    }

    /// Return the distance between two points computed using the minimum image
    /// convention.
    ///
    /// Reference
    /// ---------
    /// - Tuckerman, M. E. Statistical Mechanics: Theory and Molecular
    /// Simulation, 1 edition.; Oxford University Press: Oxford ; New York,
    /// 2010.
    pub(crate) fn distance_tuckerman(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let pij = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];

        let pmic = self.apply_mic_tuckerman(pij);
        pmic.norm()
    }
}
// Tuckerman algorithm:1 ends here

// brute force

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*brute force][brute force:1]]
impl Lattice {
    /// Return the minimal number of images for neighborhood search on each cell
    /// direction within cutoff radius
    fn n_min_images(&mut self, radius: f64) -> [isize; 3] {
        let mut ns = [0; 3];

        for (i, &w) in self.widths().iter().enumerate() {
            let n = (radius / w).ceil();
            ns[i] = n as isize;
        }

        ns
    }

    /// Return the shortest distance between `pi` (point i) and the periodic
    /// images of `pj` (point j). This algorithm will loop over all relevant
    /// images
    pub(crate) fn distance_brute_force(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let v = Vector3f::from(pj) - Vector3f::from(pi);
        let pmic = self.apply_mic_brute_force(v.into());

        pmic.norm()
    }

    /// Return the relevant periodic images required for neighborhood search
    /// within cutoff radius
    pub(crate) fn relevant_images(&mut self, radius: f64) -> Vec<Vector3f> {
        let ns = self.n_min_images(radius);
        let na = ns[0] as isize;
        let nb = ns[1] as isize;
        let nc = ns[2] as isize;

        let mut images = vec![];
        for i in (-na)..(na + 1) {
            for j in (-nb)..(nb + 1) {
                for k in (-nc)..(nc + 1) {
                    let v = Vector3f::from([i as f64, j as f64, k as f64]);
                    images.push(v);
                }
            }
        }

        images
    }

    /// Return the mic vector and its length. This algorithm will loop over all
    /// relevant images.
    pub(crate) fn apply_mic_brute_force(&mut self, p: [f64; 3]) -> Vector3f {
        // Calculate the cutoff radius for relevant images.
        // Use the value from Tuckerman algorithm as cutoff radius, since it is
        // always larger than the real distance using minimum image convention
        let p = self.apply_mic_tuckerman(p);
        let cutoff = p.norm();
        let [na, nb, nc] = self.n_min_images(cutoff);

        // search the MIC point with minimum length among relevant image points
        let points = vec![[p[0], p[1], p[2]]];
        let relevant_points = self.replicate(&points, -na..=na, -nb..=nb, -nc..=nc);
        let distances2: Vec<_> = relevant_points
            .iter()
            .map(|&p| Vector3f::from(p).norm_squared())
            .collect();
        let i = distances2.imin();
        let pmin = relevant_points[i];

        pmin.into()
    }
}
// brute force:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*test][test:1]]
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
    let dij_naive = lattice.distance_tuckerman(pi, pj);
    // it is ok: 1.2270 < 1.4142
    assert_relative_eq!(dij_naive, 1.2270, epsilon = 1e-4);
    let dij_brute = lattice.distance_brute_force(pi, pj);
    assert_relative_eq!(dij_naive, dij_brute, epsilon = 1e-4);

    // When true mic distance out of the safe range for Tuckermann algorithm.
    let pj = [-0.0834941, 1.8252187, -1.5169388];
    let dij_brute = lattice.distance_brute_force(pi, pj);
    assert_relative_eq!(dij_brute, 1.8167, epsilon = 1e-4);
    // tuckerman algo will fail since: 1.8167 > 1.4142
    let dij_naive = lattice.distance_tuckerman(pi, pj);
    assert!(dij_naive > dij_brute);
    let dij = lattice.distance(pi, pj);
    assert_relative_eq!(dij_brute, dij, epsilon = 1e-4);
}

#[test]
fn test_mic_vector() {
    let mut lat = Lattice::new([
        [7.055000000, 0.000000, 0.00000000],
        [0.000000000, 6.795000, 0.00000000],
        [-1.14679575, 0.000000, 5.65182701],
    ]);

    // mic vector
    let expected = Vector3f::from([-0.48651737, 0.184824, -1.31913642]);
    let pmic = lat.apply_mic_tuckerman([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic, epsilon = 1e-4);

    let pmic = lat.apply_mic([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic, epsilon = 1e-4);
}

#[test]
fn test_mic_distance_2() {
    let mut lat = Lattice::new([[5.0, 0.0, 0.0], [1.0, 5.0, 0.0], [1.0, 1.0, 5.0]]);

    // the shortest distance: 2.61383
    let d = lat.distance_tuckerman([0.; 3], [-0.94112, -4.34823, 2.53058]);
    assert_relative_eq!(2.66552, d, epsilon = 1e-4);
    let d = lat.distance_brute_force([0.; 3], [-0.94112, -4.34823, 2.53058]);
    assert_relative_eq!(2.61383, d, epsilon = 1e-4);

    // the shortest distance: 2.53575
    let d = lat.distance_tuckerman([0.; 3], [-2.46763, 0.57717, 0.08775]);
    assert_relative_eq!(2.59879, d, epsilon = 1e-4);
    let d = lat.distance_brute_force([0.; 3], [-2.46763, 0.57717, 0.08775]);
    assert_relative_eq!(2.53575, d, epsilon = 1e-4);
}

#[test]
fn test_neighborhood() {
    let mut lat = Lattice::new([[18.256, 0., 0.], [0., 20.534, 0.], [0., 0., 15.084]]);
    assert_eq!(true, lat.is_orthorhombic());

    assert_eq!([1, 1, 1], lat.n_min_images(9.));
    assert_eq!([2, 1, 2], lat.n_min_images(19.));
    assert_eq!([2, 1, 2], lat.n_min_images(20.));
    assert_eq!([2, 2, 2], lat.n_min_images(20.6));

    let expected = [
        Vector3f::new(-1.0, -1.0, -1.0),
        Vector3f::new(-1.0, -1.0, 0.0),
        Vector3f::new(-1.0, -1.0, 1.0),
        Vector3f::new(-1.0, 0.0, -1.0),
        Vector3f::new(-1.0, 0.0, 0.0),
        Vector3f::new(-1.0, 0.0, 1.0),
        Vector3f::new(-1.0, 1.0, -1.0),
        Vector3f::new(-1.0, 1.0, 0.0),
        Vector3f::new(-1.0, 1.0, 1.0),
        Vector3f::new(0.0, -1.0, -1.0),
        Vector3f::new(0.0, -1.0, 0.0),
        Vector3f::new(0.0, -1.0, 1.0),
        Vector3f::new(0.0, 0.0, -1.0),
        Vector3f::new(0.0, 0.0, 0.0),
        Vector3f::new(0.0, 0.0, 1.0),
        Vector3f::new(0.0, 1.0, -1.0),
        Vector3f::new(0.0, 1.0, 0.0),
        Vector3f::new(0.0, 1.0, 1.0),
        Vector3f::new(1.0, -1.0, -1.0),
        Vector3f::new(1.0, -1.0, 0.0),
        Vector3f::new(1.0, -1.0, 1.0),
        Vector3f::new(1.0, 0.0, -1.0),
        Vector3f::new(1.0, 0.0, 0.0),
        Vector3f::new(1.0, 0.0, 1.0),
        Vector3f::new(1.0, 1.0, -1.0),
        Vector3f::new(1.0, 1.0, 0.0),
        Vector3f::new(1.0, 1.0, 1.0),
    ];

    let images = lat.relevant_images(3.0);
    assert_eq!(expected.len(), images.len());
    assert_eq!(expected[1][2], images[1][2]);
}
// test:1 ends here
