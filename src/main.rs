#![allow(unused_macros, dead_code)]
use fehler::throws;
use std::fmt::{self, Debug, Write as FmtWrite};
use std::io::{self, prelude::*};
use std::ops::{Deref, Index, IndexMut};

const TAB: &'static str = "    ";

type Dec = fraction::GenericFraction<usize>;

type Error = io::Error;

trait ToTex {
    fn to_tex(&self, f: &mut dyn FmtWrite) -> fmt::Result;

    fn to_tex_str(&self) -> String {
        let mut s = String::new();
        self.to_tex(&mut s).unwrap();
        s
    }
}

impl ToTex for Dec {
    fn to_tex(&self, f: &mut dyn FmtWrite) -> fmt::Result {
        let fac = match (
            self.is_sign_positive(),
            self.numer().unwrap(),
            self.denom().unwrap(),
        ) {
            (true, n, 1) => n.to_string(),
            (true, n, d) => format!(r"\frac{{{}}}{{{}}}", n, d),
            (false, n, 1) => format!("-{}", n),
            (false, n, d) => format!(r"-\frac{{{}}}{{{}}}", n, d),
        };
        f.write_str(&fac)
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Rowop {
    Mult(usize, Dec),
    Add(Dec, usize, usize),
}

impl ToTex for Rowop {
    fn to_tex(&self, w: &mut dyn FmtWrite) -> fmt::Result {
        use Rowop::*;
        match self {
            Mult(row, by) => {
                let f = if by.is_sign_positive() {
                    by.to_tex_str()
                } else {
                    format!("({})", by.to_tex_str())
                };
                write!(w, r"\mult{{{}}}{{\scriptstyle \cdot {}}}", row, f)
            }
            Add(fac, row1, row2) => write!(
                w,
                r"\add[{sign}{}]{{{}}}{{{}}}",
                fac.to_tex_str(),
                row1,
                row2,
                sign = if fac.is_sign_positive() { "+" } else { "" }
            ),
        }
    }
}

#[derive(PartialEq)]
struct Row(Vec<Dec>);

impl Row {
    fn leading_zeroes(&self) -> usize {
        self.iter().take_while(|&&e| e == Dec::from(0)).count()
    }

    fn leading_non_zero(&self) -> Option<(usize, Dec)> {
        let zs = self.leading_zeroes();
        if zs >= self.0.len() {
            None
        } else {
            Some((zs, self.0[zs]))
        }
    }

    fn iter(&self) -> std::slice::Iter<'_, Dec> {
        self.0.iter()
    }

    fn iter_mut(&mut self) -> std::slice::IterMut<'_, Dec> {
        self.0.iter_mut()
    }
}

impl Deref for Row {
    type Target = [Dec];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Debug for Row {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            self.to_tex(f)
        } else {
            f.debug_list().entries(self.iter()).finish()
        }
    }
}

impl<Idx: std::slice::SliceIndex<[Dec]>> Index<Idx> for Row {
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {
        self.0.index(index)
    }
}

impl ToTex for [Dec] {
    fn to_tex(&self, w: &mut dyn FmtWrite) -> fmt::Result {
        write!(
            w,
            "{}",
            self.iter()
                .map(ToTex::to_tex_str)
                .collect::<Vec<_>>()
                .join(" & ")
        )
    }
}

impl From<Vec<Dec>> for Row {
    fn from(v: Vec<Dec>) -> Self {
        Row(v)
    }
}

fn main() {
    println!("Hello, world!");
}

struct Matrix<W: Write> {
    width: usize,
    height: usize,
    rows: Vec<Row>,
    rhs: usize,
    output: W,
}

impl<W: Write> Matrix<W> {
    fn new(rows: Vec<Row>, rhs: usize, output: W) -> Self {
        let height = rows.len();
        let width = rows[0].len();

        assert!(
            rows.iter().all(|r| r.len() == width),
            "All rows need to have the same length."
        );

        Self {
            width,
            height,
            rows,
            rhs,
            output,
        }
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.rows.swap(i, j);
    }

    fn mult(&mut self, row: usize, by: Dec) {
        self[row].iter_mut().for_each(|v| *v *= by);
    }

    fn add(&mut self, factor: Dec, row: usize, to: usize) {
        assert_ne!(row, to, "Failed to meet precondition: `row != to`");
        let row = &self[row] as *const Row;

        // This is safe because we assert that the row that is being mutated
        // isn't the same as the row we're iterating immutably.
        unsafe {
            (*row)
                .iter()
                .zip(self[to].iter_mut())
                .for_each(|(&a, b)| *b += factor * a);
        }
    }

    fn sort_by_leading_zeroes(&mut self) {
        self.rows.sort_by_key(Row::leading_zeroes)
    }

    #[throws]
    fn gauss_pass(&mut self, row: usize) {
        if let Some((lzs, lnz)) = self[row].leading_non_zero() {
            if lnz != Dec::from(1) {
                let r = lnz.recip();
                self.write_matrix(&[Rowop::Mult(row, r)])?;
                self.mult(row, r);
            }
            self.write_begin()?;
            let mut ops = Vec::new();
            for j in row + 1..self.height {
                let next_row = &self[j];
                if let Some((nlzs, nlnz)) = next_row.leading_non_zero() {
                    if lzs == nlzs {
                        let fac = -nlnz;
                        ops.push(Rowop::Add(fac, row, j));
                        self.add(fac, row, j);
                    }
                } else {
                    break;
                }
            }
            self.write_end(&ops)?;
        }
    }

    #[throws]
    fn gauss_jordan_pass(&mut self, row: usize) {
        if let Some((lzs, _)) = self[row].leading_non_zero() {
            for j in (0..row).rev() {
                let next_row = &self[j];
                let col_val = next_row[lzs];
                if col_val != Dec::from(0) {
                    self.write_matrix(&[Rowop::Add(-col_val, row, j)])?;
                    self.add(-col_val, row, j);
                }
            }
        }
    }

    #[throws]
    fn gauss_solve(&mut self) {
        for i in 0..self.height {
            self.sort_by_leading_zeroes();
            self.gauss_pass(i)?;
        }
        self.write_matrix(&[])?;
    }

    #[throws]
    fn gauss_jordan_solve(&mut self) {
        self.gauss_solve()?;
        for i in (1..self.height).rev() {
            self.gauss_jordan_pass(i)?;
        }
        self.write_matrix(&[])?;
    }

    fn transform_output<V: Write>(self, output: V) -> Matrix<V> {
        Matrix {
            width: self.width,
            height: self.height,
            rows: self.rows,
            rhs: self.rhs,
            output,
        }
    }

    fn into_output(self) -> W {
        self.output
    }
}

// Implementation of formatting methods.
impl<W: Write> Matrix<W> {
    #[throws]
    fn write_begin(&mut self) {
        if self.rhs == 0 {
            writeln!(self, r"\begin{{gmatrix}}[p]")?;
            self.write_rows()?;
        } else {
            self.write_lhs()?;
            self.write_begin_rhs()?;
        }
    }

    #[throws]
    fn write_end(&mut self, ops: &[Rowop]) {
        if !ops.is_empty() {
            self.write_rowops(ops)?;
        }
        writeln!(self, r"\end{{gmatrix}} \\")?;
    }

    #[throws]
    fn write_rows(&mut self) {
        let upto = self.height - 1;
        for row in &self.rows[..upto] {
            writeln!(self.output, r"{tab}{} \\", row.to_tex_str(), tab = TAB)?;
        }
        writeln!(self, "{tab}{}", self.rows[upto].to_tex_str(), tab = TAB)?;
    }

    #[throws]
    fn write_rowops(&mut self, ops: &[Rowop]) {
        writeln!(self, r"{tab}\rowops", tab = TAB)?;
        for op in ops {
            writeln!(self, "{tab}{}", op.to_tex_str(), tab = TAB)?;
        }
    }

    #[throws]
    fn write_matrix(&mut self, ops: &[Rowop]) {
        self.write_begin()?;
        self.write_end(ops)?;
    }

    // Formatting methods for bisected matrix.
    #[throws]
    fn write_lhs(&mut self) {
        writeln!(self, r"\left.\begin{{gmatrix}}[L]")?;
        self.write_rows_lhs()?;
        writeln!(self, r"\end{{gmatrix}}\right|\;")?;
    }

    #[throws]
    fn write_rows_lhs(&mut self) {
        let rhs = self.rhs;
        let mut iter = self.rows.iter().map(|row| &row[..rhs]);
        for row in iter.by_ref().take(self.height - 1) {
            writeln!(self.output, r"{tab}{} \\", row.to_tex_str(), tab = TAB)?;
        }
        writeln!(
            self.output,
            r"{tab}{}",
            iter.next().unwrap().to_tex_str(),
            tab = TAB
        )?;
    }

    #[throws]
    fn write_begin_rhs(&mut self) {
        writeln!(self, r"\begin{{gmatrix}}[R]")?;
        self.write_rows_rhs()?;
    }

    #[throws]
    fn write_rows_rhs(&mut self) {
        let rhs = self.rhs;
        let mut iter = self.rows.iter().map(|row| &row[rhs..]);
        for row in iter.by_ref().take(self.height - 1) {
            writeln!(self.output, r"{tab}{} \\", row.to_tex_str(), tab = TAB)?;
        }
        writeln!(
            self.output,
            r"{tab}{}",
            iter.next().unwrap().to_tex_str(),
            tab = TAB
        )?;
    }
}

impl<W: Write> Index<usize> for Matrix<W> {
    type Output = Row;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.rows[idx]
    }
}

impl<W: Write> IndexMut<usize> for Matrix<W> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.rows[idx]
    }
}

impl<W: Write> PartialEq for Matrix<W> {
    fn eq(&self, other: &Matrix<W>) -> bool {
        self.width == other.width && self.height == other.height && self.rows == other.rows
    }
}

impl<W: Write> Write for Matrix<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.output.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.output.flush()
    }
}

impl<W: Write> Debug for Matrix<W> {
    #[throws(fmt::Error)]
    fn fmt(&self, f: &mut fmt::Formatter) {
        for row in &self.rows {
            f.write_str(
                &row.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join("  "),
            )?;
        }
    }
}

macro_rules! matrix {
    ($($($x:expr),+);+) => {
        {
            let mut width = None;
            let mut rows = Vec::new();
            $(
                let mut row = Vec::new();
                $(
                    row.push(Dec::from($x));
                )+
                match width {
                    None => {
                        width.replace(row.len());
                    }
                    Some(w) => if !(w == row.len()) { panic!("Not all rows have the same length") }
                }
                rows.push(row.into());
            )+
            Matrix {
                width: width.unwrap(),
                height: rows.len(),
                rows,
                rhs: 0,
                output: std::io::sink(),
            }
        }
    };
    // ($($($x:expr),+ ;; ($y:expr),+);+) => {

    // }
}

macro_rules! rows {
    ($($($x:expr),+);+) => {
        {
            vec![$(Row::from(vec![$(Dec::from($x)),+])),+]
        }
    };
    ($($($x:expr),+);+;) => {
        {
            rows!($($($x),+);+)
        }
    };
    // ($($($x:expr),+ ;; ($y:expr),+);+) => {

    // }
}

#[cfg(test)]
mod tests {
    #![allow(unused_must_use)]
    use super::*;
    #[test]
    fn test_leading_non_zero() {
        let mut m = matrix! {
            1, 2, 3, 4;
            0, 0, 0, 0;
            0, 6, 0, 8
        };

        assert_eq!(m[0].leading_non_zero(), Some((0, Dec::from(1))));
        assert_eq!(m[1].leading_non_zero(), None);
        assert_eq!(m[2].leading_non_zero(), Some((1, Dec::from(6))));

        m.sort_by_leading_zeroes();

        assert_eq!(
            m,
            matrix! {
                1, 2, 3, 4;
                0, 6, 0, 8;
                0, 0, 0, 0
            }
        );
    }

    #[test]
    fn test_gauss_pass() {
        let mut m = matrix! {
            0, 2, 4, 6, 8;
            0, 3, 6, 9, 13;
            0, 5.4, 10.8, 9, 13;
            0, 0, 1, 0, 0
        };

        m.gauss_pass(0).unwrap();

        assert_eq!(
            m,
            matrix! {
                0, 1, 2, 3, 4;
                0, 0, 0, 0, 1;
                0, 0, 0, -7.2, -8.6;
                0, 0, 1, 0, 0
            }
        );
    }

    #[test]
    fn test_gauss_jordan_pass() {
        let mut m = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        m.gauss_solve();

        assert_eq!(
            m,
            matrix! {
                1, 0, 1.5;
                0, 1, -0.75;
                0, 0, 1;
                0, 0, 0
            }
        );

        m.gauss_jordan_pass(2);

        assert_eq!(
            m,
            matrix! {
                1, 0, 0;
                0, 1, 0;
                0, 0, 1;
                0, 0, 0
            }
        );
    }

    #[test]
    fn test_gauss_jordan_solve() {
        let mut m = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        m.gauss_jordan_solve();

        assert_eq!(
            m,
            matrix! {
                1, 0, 0;
                0, 1, 0;
                0, 0, 1;
                0, 0, 0
            }
        );
    }

    #[test]
    fn test_gauss_solve() {
        let mut m = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        let mut m2 = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        m.sort_by_leading_zeroes();
        m.gauss_pass(0);

        assert_eq!(
            m,
            matrix! {
                1, 0, 1.5;
                0, -2, 1.5;
                0, 1, -2;
                0, 1, -2
            }
        );

        m.sort_by_leading_zeroes();
        m.gauss_pass(1);
        assert_eq!(
            m,
            matrix! {
                1, 0, 1.5;
                0, 1, -0.75;
                0, 0, -1.25;
                0, 0, -1.25
            }
        );

        m.sort_by_leading_zeroes();
        m.gauss_pass(2);
        assert_eq!(
            m,
            matrix! {
                1, 0, 1.5;
                0, 1, -0.75;
                0, 0, 1;
                0, 0, 0
            }
        );

        m.sort_by_leading_zeroes();
        m.gauss_pass(3);
        assert_eq!(
            m,
            matrix! {
                1, 0, 1.5;
                0, 1, -0.75;
                0, 0, 1;
                0, 0, 0
            }
        );

        m2.gauss_solve();

        assert_eq!(m2, m);
        assert_eq!(
            m2,
            matrix! {
                1, 0, 1.5;
                0, 1, -0.75;
                0, 0, 1;
                0, 0, 0
            }
        );
    }

    #[test]
    fn test_output() -> std::io::Result<()> {
        let f = std::fs::File::create("./output.txt")?;
        let m = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        let mut m = m.transform_output(f);

        m.gauss_solve();

        Ok(())
    }

    #[test]
    fn test_bisected_matrix() {
        let m = matrix! {
            1, 2, 3;
            2, 3, 4;
            1, 2, 3
        };

        let mut m = m.transform_output(io::stdout());
        m.rhs = 2;
        m.gauss_pass(0);
    }

    #[test]
    fn test_rows_macro() {
        let rows: Vec<Row> = rows! {
            1,2,3;
            2,3,4;
            3,4,5
        };

        let rows2: Vec<Row> = rows! {
            1,2,3;
            2,3,4;
            3,4,5;
        };

        assert_eq!(
            rows,
            vec![
                Row::from(vec![1.into(), 2.into(), 3.into()]),
                Row::from(vec![2.into(), 3.into(), 4.into()]),
                Row::from(vec![3.into(), 4.into(), 5.into()]),
            ]
        );
        assert_eq!(rows, rows2);
    }
}
