// base

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*base][base:1]]
use crate::Lattice;
use guts::itertools::*;

type Point = [f64; 3];

impl Lattice {
    /// Create a supercell along three cell directions.
    pub fn replicate(
        &self,
        points: &[Point],
        ra: impl Iterator<Item = isize> + Clone,
        rb: impl Iterator<Item = isize> + Clone,
        rc: impl Iterator<Item = isize> + Clone,
    ) -> Vec<Point> {
        iproduct!(ra, rb, rc)
            .flat_map(|(i, j, k)| {
                let v = [i as f64, j as f64, k as f64];
                let [tx, ty, tz] = self.to_cart(v);
                points
                    .iter()
                    .map(move |[x0, y0, z0]| [x0 + tx, y0 + ty, z0 + tz])
            })
            .collect()
    }
}
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

    let mut lattice = Lattice::new(cell);
    let points = vec![[0.0; 3]];
    let new_points = lattice.replicate(&points, -1..=1, -1..=1, -1..=1);
    assert_eq!(new_points.len(), 27);
}
// test:1 ends here
