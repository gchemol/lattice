// mic.rs
// :PROPERTIES:
// :header-args: :tangle tests/mic.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*mic.rs][mic.rs:1]]
use lattice::Lattice;

use approx::*;
use vecfx::*;

#[test]
fn test_mic_distance() {
    // Setup lattice
    // a = b = c = 4, alpha = beta = gamma = 60
    let cell = [
        [4.00000000, 0.00000000, 0.00000000],
        [2.00000000, 3.46410162, 0.00000000],
        [2.00000000, 1.15470054, 3.26598632],
    ];
    let lattice = Lattice::new(cell);

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
// mic.rs:1 ends here
