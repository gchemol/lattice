// core

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*core][core:1]]
use vecfx::*;

// matrix inversion
pub(crate) fn get_inv_matrix(matrix: Matrix3f) -> Matrix3f {
    matrix.try_inverse().expect("bad matrix")
}

// cell volume
pub(crate) fn get_cell_volume(mat: Matrix3f) -> f64 {
    let va = mat.column(0);
    let vb = mat.column(1);
    let vc = mat.column(2);
    va.dot(&vb.cross(&vc))
}

// return cell length parameters
pub(crate) fn get_cell_lengths(mat: Matrix3f) -> [f64; 3] {
    [
        mat.column(0).norm(),
        mat.column(1).norm(),
        mat.column(2).norm(),
    ]
}

// return cell angle parameters in degrees
pub(crate) fn get_cell_angles(mat: Matrix3f) -> [f64; 3] {
    let va = mat.column(0);
    let vb = mat.column(1);
    let vc = mat.column(2);
    [
        vb.angle(&vc).to_degrees(),
        va.angle(&vc).to_degrees(),
        va.angle(&vb).to_degrees(),
    ]
}
// core:1 ends here
