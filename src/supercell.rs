// base

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*base][base:1]]
use crate::Lattice;
use guts::itertools::*;
use vecfx::Vector3f;

impl Lattice {
    /// Create a supercell along three cell directions.
    pub fn replicate(
        &self,
        ra: impl Iterator<Item = isize> + Clone,
        rb: impl Iterator<Item = isize> + Clone,
        rc: impl Iterator<Item = isize> + Clone,
    ) -> impl Iterator<Item = Vector3f> {
        iproduct!(ra, rb, rc).map(|(i, j, k)| Vector3f::from([i as f64, j as f64, k as f64]))
    }
}

// #[cfg(feature = "adhoc")]
// /// Helper struct for periodic image
// #[derive(Clone, Copy, Debug)]
// pub struct Image(pub(crate) isize, pub(crate) isize, pub(crate) isize);

// #[cfg(feature = "adhoc")]
// impl Image {
//     /// Return fractional translation vector for moving a point into this image.
//     pub fn translation_vector(&self) -> Vector3f {
//         [self.0 as f64, self.1 as f64, self.2 as f64].into()
//     }

//     /// Return image location relative to origin cell.
//     pub fn location(&self) -> [isize; 3] {
//         [self.0, self.1, self.2]
//     }
// }
// base:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*test][test:1]]
#[test]
fn test_supercell() {
    // Setup lattice
    // a = b = c = 4, alpha = beta = gamma = 60
    let cell = [
        [4.00000000, 0.00000000, 0.00000000],
        [2.00000000, 3.46410162, 0.00000000],
        [2.00000000, 1.15470054, 3.26598632],
    ];

    let lattice = Lattice::new(cell);
    let cell_images = lattice.replicate(-1..=1, -1..=1, -1..=1);
    assert_eq!(cell_images.count(), 27);
}
// test:1 ends here
