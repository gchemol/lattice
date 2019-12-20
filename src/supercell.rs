// base

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*base][base:1]]
use crate::Lattice;
use guts::itertools::*;
use vecfx::Vector3f;

impl Lattice {
    /// Create a supercell along three cell directions.
    pub fn replicate<T: Into<Vector3f> + Copy>(
        &self,
        points: &[T],
        ra: impl Iterator<Item = isize> + Clone,
        rb: impl Iterator<Item = isize> + Clone,
        rc: impl Iterator<Item = isize> + Clone,
    ) -> Vec<Vector3f> {
        iproduct!(ra, rb, rc)
            .flat_map(|(i, j, k)| {
                let v = [i as f64, j as f64, k as f64];
                let tv = self.to_cart(v);
                points.iter().map(move |&p| tv + p.into())
            })
            .collect()
    }

    #[cfg(feature = "adhoc")]
    /// Create a supercell along three cell directions.
    pub fn replicate_images(
        &self,
        ra: impl Iterator<Item = isize> + Clone,
        rb: impl Iterator<Item = isize> + Clone,
        rc: impl Iterator<Item = isize> + Clone,
    ) -> impl Iterator<Item = Image> {
        iproduct!(ra, rb, rc).map(|(i, j, k)| Image(i, j, k))
    }
}

#[cfg(feature = "adhoc")]
/// Helper struct for periodic image
#[derive(Clone, Copy, Debug)]
pub struct Image(pub isize, pub isize, pub isize);
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
    let points = vec![[0.0; 3]];
    let new_points = lattice.replicate(&points, -1..=1, -1..=1, -1..=1);
    assert_eq!(new_points.len(), 27);
}
// test:1 ends here
