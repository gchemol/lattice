// imports

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*imports][imports:1]]
use vecfx::*;

#[cfg(test)]
use approx::*;

use crate::Lattice;
// imports:1 ends here

// base
// The periodic image when periodic boundary conditions are employed.

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*base][base:1]]
use std::f64;

#[derive (Debug, Clone)]
pub struct PeriodicImage {
    /// cartesian positions of the particle image
    pub position: Vector3f,
    /// scaled displacment vector relative to origin cell
    pub image   : Vector3f,
}
// base:1 ends here

// distance

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*distance][distance:1]]
impl Lattice {
    /// Return the shortest vector by applying the minimum image convention.
    pub fn apply_mic(&mut self, p: [f64; 3]) -> PeriodicImage {
        if self.is_orthorhombic() {
            self.apply_mic_tuckerman(p)
        } else {
            self.apply_mic_brute_force(p)
        }
    }

    /// Return the mic vector using Tuckerman's algorithm.
    ///
    /// Reference
    /// ---------
    /// - Tuckerman, M. E. Statistical Mechanics: Theory and Molecular
    /// Simulation, 1 edition.; Oxford University Press: Oxford ; New York,
    /// 2010.
    pub fn apply_mic_tuckerman(&mut self, p: [f64; 3]) -> PeriodicImage {
        // apply minimum image convention on the scaled coordinates
        let mut fcoords = self.to_frac(p);

        let mut image = [1.0; 3];
        for i in 0..3 {
            image[i] = -1.0 * fcoords[i].round();
            fcoords[i] += image[i];
        }

        // transform back to cartesian coordinates
        let pij = self.to_cart(fcoords);

        PeriodicImage {
            position: Vector3f::from(pij),
            image   : image.into(),
        }
    }

    // FIXME: remove type conversion
    /// Return the mic vector. This algorithm will loop over all relevant images
    pub fn apply_mic_brute_force(&mut self, p: [f64; 3]) -> PeriodicImage {
        // The cutoff radius for finding relevant images.
        // Use the value from Tuckerman algorithm as cutoff radius, since it is
        // always larger than the real distance using minimum image convention
        let cutoff = self.apply_mic_tuckerman(p).position.norm();
        let relevant_images = self.relevant_images(cutoff);

        // tuple = (distance, position, image)
        let mut target = (f64::MAX,
                          Vector3f::from([0.0; 3]),
                          Vector3f::from([0.0; 3]));
        for image in relevant_images {
            let dd = self.to_cart(image.into());
            let ip = [
                p[0] + dd[0],
                p[1] + dd[1],
                p[2] + dd[2],
            ];

            let v = Vector3f::from(ip);
            let d = v.norm();
            if d < target.0 {
                target = (d, v, image);
            }
        }

        PeriodicImage {
            position: target.1,
            image   : target.2,
        }
    }

    /// Return the relevant periodic images required for neighborhood search
    /// within cutoff radius
    pub fn relevant_images(&mut self, radius: f64) -> Vec<Vector3f> {
        let ns = self.n_min_images(radius);
        let na = ns[0] as isize;
        let nb = ns[1] as isize;
        let nc = ns[2] as isize;

        let mut images = vec![];
        for i in (-na)..(na+1) {
            for j in (-nb)..(nb+1) {
                for k in (-nc)..(nc+1) {
                    let v = Vector3f::from([i as f64, j as f64, k as f64]);
                    images.push(v);
                }
            }
        }

        images
    }

    /// Return the minimal number of images for neighborhood search on each cell
    /// direction within cutoff radius
    fn n_min_images(&mut self, radius: f64) -> [usize; 3]{
        let mut ns = [0; 3];

        for (i, &w) in self.widths().iter().enumerate() {
            let n = (radius / w).ceil();
            ns[i] = n as usize;
        }

        ns
    }

    /// Return the distance between two points computed using the minimum image
    /// convention.
    ///
    /// Reference
    /// ---------
    /// - Tuckerman, M. E. Statistical Mechanics: Theory and Molecular
    /// Simulation, 1 edition.; Oxford University Press: Oxford ; New York,
    /// 2010.
    fn distance_tuckerman(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let pij = [pj[0] - pi[0],
                   pj[1] - pi[1],
                   pj[2] - pi[2]];

        let pmic = self.apply_mic_tuckerman(pij);
        pmic.position.norm()
    }

    /// Return the shortest distance between `pi` (point i) and the periodic
    /// images of `pj` (point j). This algorithm will loop over all relevant
    /// images
    fn distance_brute_force(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let v = Vector3f::from(pj) - Vector3f::from(pi);
        let pmic = self.apply_mic_brute_force(v.into());

        pmic.position.norm()
    }

    // TODO: return the nearest periodic image?
    /// Return the shortest distance between `pi` (point i) and the periodic
    /// images of `pj` (point j) under the minimum image convention
    pub fn distance(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let pmic = self.apply_mic([pj[0] - pi[0],
                                   pj[1] - pi[1],
                                   pj[2] - pi[2]]);
        pmic.position.norm()
    }
}
// distance:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*test][test:1]]
#[test]
fn test_mic_vector() {
    let mut lat = Lattice::new([
        [7.055, 0., 0.],
        [0., 6.795, 0.],
        [-1.14679575, 0., 5.65182701],
    ]);

    // mic vector
    let expected = Vector3f::from([-0.48651737, 0.184824, -1.31913642]);

    let pmic = lat.apply_mic_tuckerman([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic.position, epsilon = 1e-4);

    assert_relative_eq!(
        pmic.image,
        Vector3f::from([-1.0, 0.0, -1.0]),
        epsilon = 1e-4
    );

    let pmic = lat.apply_mic([5.42168688, 0.184824, 4.33269058]);
    assert_relative_eq!(expected, pmic.position, epsilon = 1e-4);
    assert_relative_eq!(
        pmic.image,
        Vector3f::from([-1.0, 0.0, -1.0]),
        epsilon = 1e-4
    );
}

#[test]
fn test_lattice_mic_distance() {
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
fn test_lattice_neighborhood() {
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
