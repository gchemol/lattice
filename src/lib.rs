// [[file:../lattice.note::*header][header:1]]
//===============================================================================#
//   DESCRIPTION:  Represents 3D periodic lattice
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-29 14:27>
//       UPDATED:  <>
//===============================================================================#
// header:1 ends here

// [[file:../lattice.note::*imports][imports:1]]
use gchemol_gut::prelude::*;
use vecfx::*;
// imports:1 ends here

// [[file:../lattice.note::*mods][mods:1]]
mod mic;
mod supercell;
mod utils;

use crate::utils::*;
// mods:1 ends here

// [[file:../lattice.note::*base][base:1]]
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

// [[file:../lattice.note::f072864d][f072864d]]
impl Lattice {
    /// Construct `Lattice` from three lattice vectors.
    pub fn new<T: Into<Vector3f> + Copy>(tvs: [T; 3]) -> Self {
        let vectors = [tvs[0].into(), tvs[1].into(), tvs[2].into()];
        let matrix = Matrix3f::from_columns(&vectors);
        Self::from_matrix(matrix)
    }

    /// Construct `Lattice` from lattice matrix (3x3).
    pub fn from_matrix<T: Into<Matrix3f>>(tvs: T) -> Self {
        let matrix = tvs.into();
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
    pub fn set_origin<T: Into<Vector3f>>(&mut self, loc: T) {
        self.origin = loc.into()
    }

    /// Lattice length parameters: a, b, c
    pub fn lengths(&self) -> [f64; 3] {
        get_cell_lengths(self.matrix).into()
    }

    /// Lattice angle parameters in degrees
    pub fn angles(&self) -> [f64; 3] {
        get_cell_angles(self.matrix).into()
    }

    /// Scale Lattice by a positive constant `v`
    pub fn scale_by(&mut self, v: f64) {
        assert!(v.is_sign_positive(), "invalid scale factor: {v}");
        self.matrix *= v;
        self.inv_matrix = get_inv_matrix(&self.matrix);
    }

    /// Scale Lattice in `a` direction by a positive constant `v`
    pub fn scale_by_a(&mut self, v: f64) {
        self.scale_by_abc(v, 0)
    }

    /// Scale Lattice in `b` direction by a positive constant `v`
    pub fn scale_by_b(&mut self, v: f64) {
        self.scale_by_abc(v, 1)
    }

    /// Scale Lattice in `c` direction by a positive constant `v`
    pub fn scale_by_c(&mut self, v: f64) {
        self.scale_by_abc(v, 2)
    }

    /// Scale Lattice in length direction by a positive constant
    fn scale_by_abc(&mut self, v: f64, i: usize) {
        assert!(v.is_sign_positive(), "invalid scale factor: {v}");
        let mut x = self.matrix.column_mut(i);
        x *= v;
        self.inv_matrix = get_inv_matrix(&self.matrix);
    }

    /// Get cell origin in Cartesian coordinates
    pub fn origin(&self) -> Vector3f {
        self.origin
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac<T: Into<Vector3f>>(&self, p: T) -> Vector3f {
        self.inv_matrix * (p.into() - self.origin)
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    pub fn to_cart<T: Into<Vector3f>>(&self, p: T) -> Vector3f {
        self.matrix * p.into() + self.origin
    }

    /// Lattice vector a
    pub fn vector_a(&self) -> Vector3f {
        self.matrix.column(0).into()
    }

    /// Lattice vector b
    pub fn vector_b(&self) -> Vector3f {
        self.matrix.column(1).into()
    }

    /// Lattice vector c
    pub fn vector_c(&self) -> Vector3f {
        self.matrix.column(2).into()
    }

    /// Three lattice vectors.
    pub fn vectors(&self) -> [Vector3f; 3] {
        [self.vector_a(), self.vector_b(), self.vector_c()]
    }

    /// Lattice vectors
    pub fn matrix(&self) -> Matrix3f {
        self.matrix
    }

    /// inverse of lattice matrix
    pub fn inv_matrix(&self) -> Matrix3f {
        self.inv_matrix
    }

    /// Check if lattice is orthorhombic
    pub fn is_orthorhombic(&self) -> bool {
        let diag = self.matrix.diagonal();
        let m = Matrix3f::from_diagonal(&diag);
        m == self.matrix
    }

    /// Wrap a point in cartesian coordinates into unit cell, obeying the
    /// periodic boundary conditions. Returns cartesian coordinates.
    pub fn wrap<T: Into<Vector3f>>(&self, vec: T) -> Vector3f {
        let f = self.to_frac(vec);
        let fcoords_wrapped = self.wrap_frac(f);
        self.to_cart(fcoords_wrapped)
    }

    /// Wrap a point in fractional coordinates into unit cell, obeying the
    /// periodic boundary conditions. Returns fractional coordinates.
    pub fn wrap_frac<T: Into<Vector3f>>(&self, f: T) -> Vector3f {
        let f = f.into();
        let fcoords_wrapped = [f.x - f.x.floor(), f.y - f.y.floor(), f.z - f.z.floor()];
        fcoords_wrapped.into()
    }

    /// Return the shortest distance between `pi` (point i) and the periodic
    /// images of `pj` (point j) under the minimum image convention
    ///
    /// Parameters
    /// ----------
    /// * pi, pj: Cartesian coordinates of point i and point j
    pub fn distance<T: Into<Vector3f>>(&self, pi: T, pj: T) -> f64 {
        let p = pj.into() - pi.into();
        let pmic = self.apply_mic(p);
        pmic.norm()
    }

    /// Return the shortest vector obeying the minimum image convention.
    pub fn apply_mic<T: Into<[f64; 3]>>(&self, p: T) -> Vector3f {
        let p = p.into();
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
// f072864d ends here
