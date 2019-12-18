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
//       UPDATED:  <2019-12-18 Wed 09:30>
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

    /// Cached inverse of lattice matrix
    inv_matrix: Option<Matrix3f>,

    /// Cached volume of the unit cell.
    volume: Option<f64>,

    /// The perpendicular widths of the unit cell on each direction,
    /// i.e. the distance between opposite faces of the unit cell
    widths: Option<[f64; 3]>,

    /// Cached cell lengths parameters
    lengths: Option<[f64; 3]>,

    /// Cached cell angles parameters
    angles: Option<[f64; 3]>,
}

impl Default for Lattice {
    fn default() -> Self {
        Lattice {
            matrix: Matrix3f::identity(),
            origin: Vector3f::zeros(),

            inv_matrix: None,
            volume: None,
            widths: None,
            lengths: None,
            angles: None,
        }
    }
}

impl Lattice {
    pub fn new<T: Into<[[f64; 3]; 3]>>(tvs: T) -> Self {
        Lattice {
            matrix: Matrix3f::from(tvs.into()),
            ..Default::default()
        }
    }

    /// using a cache to reduce the expensive matrix inversion calculations
    fn inv_matrix(&mut self) -> Matrix3f {
        // make a readonly reference
        let matrix = self.matrix;
        let im = self
            .inv_matrix
            .get_or_insert_with(|| matrix.try_inverse().expect("bad matrix"));

        *im
    }

    fn get_cell_widths(&mut self) -> [f64; 3] {
        let volume = self.volume();
        let [van, vbn, vcn] = self.lengths();

        let wa = volume / (vbn * vcn);
        let wb = volume / (vcn * van);
        let wc = volume / (van * vbn);

        [wa, wb, wc]
    }

    /// Return the perpendicular widths of the cell along three directions.
    pub fn widths(&mut self) -> [f64; 3] {
        if let Some(ws) = self.widths {
            return ws;
        } else {
            let ws = self.get_cell_widths();
            self.widths = Some(ws);

            ws
        }
    }

    /// Return the volume of the unit cell
    /// the cache will be updated if necessary
    pub fn volume(&mut self) -> f64 {
        // make a read-only reference
        let mat = self.matrix;
        let volume = self.volume.get_or_insert_with(|| get_cell_volume(mat));

        *volume
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

    /// Set cell origin in Cartesian coordinates
    pub fn set_origin(&mut self, loc: [f64; 3]) {
        self.origin = Vector3f::from(loc);
    }

    /// Lattice length parameters: a, b, c
    pub fn lengths(&mut self) -> [f64; 3] {
        let mat = self.matrix;
        let lengths = self.lengths.get_or_insert_with(|| get_cell_lengths(mat));

        [lengths[0], lengths[1], lengths[2]]
    }

    /// Lattice angle parameters in degrees
    pub fn angles(&mut self) -> [f64; 3] {
        let mat = self.matrix;
        let angles = self.angles.get_or_insert_with(|| get_cell_angles(mat));

        [angles[0], angles[1], angles[2]]
    }

    // FIXME: cell widths
    /// Scale Lattice by a positive constant
    pub fn scale_by(&mut self, v: f64) {
        debug_assert!(v > 0.);
        self.matrix *= v;

        // reset caches
        self.inv_matrix = None;
        self.volume = None;
        self.widths = None;
        self.lengths = None;
        self.angles = None;
    }

    /// Get cell origin in Cartesian coordinates
    pub fn origin(&self) -> [f64; 3] {
        self.origin.into()
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac(&mut self, p: [f64; 3]) -> [f64; 3] {
        let im = self.inv_matrix();
        let v = Vector3f::from(p);
        let fs = im * (v - self.origin);
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
}
// base:1 ends here
