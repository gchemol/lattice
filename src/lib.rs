// header

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*header][header:1]]
//===============================================================================#
//   DESCRIPTION:  Represents 3D periodic lattice
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-29 14:27>
//       UPDATED:  <2019-12-18 Wed 21:14>
//===============================================================================#
// header:1 ends here

// imports

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*imports][imports:1]]
use guts::prelude::*;
use vecfx::*;
// imports:1 ends here

// mods

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*mods][mods:1]]
mod mic;
mod utils;
mod supercell;

use crate::utils::*;
// mods:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*base][base:1]]
/// Periodic 3D lattice
#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct Lattice {
    /// internal translation matrix
    matrix: Matrix3f,
    /// Lattice origin
    origin: Vector3f,
    /// inverse of lattice matrix
    inv_matrix: Matrix3f,
}

impl Default for Lattice {
    fn default() -> Self {
        let matrix = Matrix3f::identity();
        let inv_matrix = get_inv_matrix(&matrix);
        Lattice {
            matrix,
            inv_matrix,
            origin: Vector3f::zeros(),
        }
    }
}
// base:1 ends here

// api

// [[file:~/Workspace/Programming/gchemol-rs/lattice/lattice.note::*api][api:1]]
impl Lattice {
    pub fn new<T: Into<[[f64; 3]; 3]>>(tvs: T) -> Self {
        let matrix = Matrix3f::from(tvs.into());
        let inv_matrix = get_inv_matrix(&matrix);
        Lattice {
            matrix,
            inv_matrix,
            ..Default::default()
        }
    }

    /// Construct lattice from lattice parameters
    /// Unit cell angles in degrees, lengths in Angstrom
    pub fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let alpha = alpha.to_radians();
        let beta = beta.to_radians();
        let gamma = gamma.to_radians();

        let acos = alpha.cos();
        let bcos = beta.cos();
        let gcos = gamma.cos();
        let gsin = gamma.sin();
        let v = (1. - acos.powi(2) - bcos.powi(2) - gcos.powi(2) + 2.0 * acos * bcos * gcos).sqrt();

        let va = [a, 0.0, 0.0];

        let vb = [b * gcos, b * gsin, 0.0];

        let vc = [c * bcos, c * (acos - bcos * gcos) / gsin, c * v / gsin];

        Lattice::new([va, vb, vc])
    }

    /// Return the perpendicular widths of the cell along three directions. i.e.
    /// the distance between opposite faces of the unit cell
    pub fn widths(&self) -> [f64; 3] {
        let volume = self.volume();
        let [van, vbn, vcn] = self.lengths();

        let wa = volume / (vbn * vcn);
        let wb = volume / (vcn * van);
        let wc = volume / (van * vbn);

        [wa, wb, wc]
    }

    /// Return the volume of the unit cell
    /// the cache will be updated if necessary
    pub fn volume(&self) -> f64 {
        get_cell_volume(self.matrix)
    }

    /// Set cell origin in Cartesian coordinates
    pub fn set_origin(&mut self, loc: [f64; 3]) {
        self.origin = Vector3f::from(loc);
    }

    /// Lattice length parameters: a, b, c
    pub fn lengths(&self) -> [f64; 3] {
        let lengths = get_cell_lengths(self.matrix);

        [lengths[0], lengths[1], lengths[2]]
    }

    /// Lattice angle parameters in degrees
    pub fn angles(&self) -> [f64; 3] {
        let angles = get_cell_angles(self.matrix);

        [angles[0], angles[1], angles[2]]
    }

    // FIXME: cell widths
    /// Scale Lattice by a positive constant
    pub fn scale_by(&mut self, v: f64) {
        debug_assert!(v > 0.);
        self.matrix *= v;
        self.inv_matrix = get_inv_matrix(&self.matrix);
    }

    /// Get cell origin in Cartesian coordinates
    pub fn origin(&self) -> [f64; 3] {
        self.origin.into()
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac(&self, p: [f64; 3]) -> [f64; 3] {
        let v = Vector3f::from(p);
        let fs = self.inv_matrix * (v - self.origin);
        fs.into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    pub fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        let v = Vector3f::from(p);
        let fs = self.matrix * v + self.origin;

        fs.into()
    }

    /// Lattice vector a
    pub fn vector_a(&self) -> [f64; 3] {
        self.matrix.column(0).transpose().into()
    }

    /// Lattice vector b
    pub fn vector_b(&self) -> [f64; 3] {
        self.matrix.column(1).transpose().into()
    }

    /// Lattice vector c
    pub fn vector_c(&self) -> [f64; 3] {
        self.matrix.column(2).transpose().into()
    }

    /// Lattice vectors
    pub fn vectors(&self) -> [[f64; 3]; 3] {
        self.matrix.into()
    }

    /// Check if lattice is orthorhombic
    pub fn is_orthorhombic(&self) -> bool {
        let diag = self.matrix.diagonal();
        let m = Matrix3f::from_diagonal(&diag);
        m == self.matrix
    }

    /// Wrap a point to unit cell, obeying the periodic boundary conditions.
    pub fn wrap(&self, vec: [f64; 3]) -> [f64; 3] {
        let [fx, fy, fz] = self.to_frac(vec);
        let fcoords_wrapped = [fx - fx.floor(), fy - fy.floor(), fz - fz.floor()];
        self.to_cart(fcoords_wrapped)
    }

    /// Return the shortest distance between `pi` (point i) and the periodic
    /// images of `pj` (point j) under the minimum image convention
    ///
    /// Parameters
    /// ----------
    /// * pi, pj: Cartesian coordinates of point i and point j
    pub fn distance(&self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let pmic = self.apply_mic([pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]]);
        pmic.norm()
    }

    /// Return the shortest vector by applying the minimum image convention.
    pub fn apply_mic(&self, p: [f64; 3]) -> Vector3f {
        // Tuckerman algorithm works well for Orthorombic cell
        let v_naive = self.apply_mic_tuckerman(p);
        if self.is_orthorhombic() {
            v_naive
        } else {
            let r_max = 0.5 * self.widths().min();
            if v_naive.norm() < r_max {
                v_naive
            } else {
                self.apply_mic_brute_force(p)
            }
        }
    }
}
// api:1 ends here
