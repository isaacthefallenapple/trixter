#![allow(unused_macros, dead_code)]
use std::fmt::{self, Debug, Write as FmtWrite};
use std::io::{self, prelude::*};
use std::ops::{Index, IndexMut};

type Dec = fraction::GenericFraction<usize>;

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

impl Debug for Row {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            self.to_tex(f)
        } else {
            f.debug_list().entries(self.iter()).finish()
        }
    }
}

impl Index<usize> for Row {
    type Output = Dec;

    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}

impl ToTex for Row {
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

struct Matrix<'a> {
    width: usize,
    height: usize,
    rows: Vec<Row>,
    rowops: Vec<Rowop>,
    output: Box<dyn Write + 'a>,
}

impl<'a> Matrix<'a> {
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

    fn gauss_pass(&mut self, row: usize) -> io::Result<()> {
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
        Ok(())
    }

    fn gauss_jordan_pass(&mut self, row: usize) -> io::Result<()> {
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
        Ok(())
    }

    fn gauss_solve(&mut self) -> io::Result<()> {
        for i in 0..self.height {
            self.sort_by_leading_zeroes();
            self.gauss_pass(i)?;
        }
        self.write_matrix(&[])
    }

    fn gauss_jordan_solve(&mut self) -> io::Result<()> {
        self.gauss_solve()?;
        for i in (1..self.height).rev() {
            self.gauss_jordan_pass(i)?;
        }
        self.write_matrix(&[])?;
        Ok(())
    }

    fn set_output(&mut self, w: impl Write + 'a) {
        self.output = Box::new(w);
    }

    fn write_begin(&mut self) -> io::Result<()> {
        writeln!(self.output, r"\begin{{gmatrix}}[p]")?;
        self.write_rows()
    }

    fn write_end(&mut self, ops: &[Rowop]) -> io::Result<()> {
        if !ops.is_empty() {
            self.write_rowops(ops)?;
        }
        writeln!(self, r"\end{{gmatrix}} \\")
    }

    fn write_rows(&mut self) -> io::Result<()> {
        let upto = self.height - 1;
        for row in &self.rows[..upto] {
            writeln!(self.output, r"  {} \\", row.to_tex_str())?;
        }
        writeln!(self, "  {}", self.rows[upto].to_tex_str())
    }

    fn write_rowops(&mut self, ops: &[Rowop]) -> io::Result<()> {
        writeln!(self, r"  \rowops")?;
        for op in ops {
            writeln!(self, "  {}", op.to_tex_str())?;
        }
        Ok(())
    }

    fn write_matrix(&mut self, ops: &[Rowop]) -> io::Result<()> {
        self.write_begin()?;
        self.write_end(ops)
    }

    fn to_tex_lines(&self) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.height + 3);
        lines.push(r"\begin{gmatrix}[p]".into());
        let upto = self.height - 1;
        for row in &self.rows[..upto] {
            lines.push(format!(r"  {:#?} \\", row));
        }
        lines.push(format!("  {:#?}", self.rows[upto]));
        lines.push(r"  \rowops".into());
        lines
    }
}

impl<'a> Index<usize> for Matrix<'a> {
    type Output = Row;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.rows[idx]
    }
}

impl<'a> IndexMut<usize> for Matrix<'a> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.rows[idx]
    }
}

impl<'a> PartialEq for Matrix<'a> {
    fn eq(&self, other: &Matrix) -> bool {
        self.width == other.width && self.height == other.height && self.rows == other.rows
    }
}

impl<'a> ToTex for Matrix<'a> {
    fn to_tex(&self, w: &mut dyn FmtWrite) -> fmt::Result {
        for line in &self.to_tex_lines() {
            writeln!(w, "{}", line)?;
        }
        Ok(())
    }
}

impl<'a> Write for Matrix<'a> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.output.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.output.flush()
    }
}

impl<'a> Debug for Matrix<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if f.alternate() {
            self.to_tex(f)
        } else {
            for row in &self.rows {
                writeln!(
                    f,
                    "{}",
                    row.iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(" | ")
                )?;
            }
            Ok(())
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
                rowops: Vec::new(),
                output: Box::new(std::io::sink()),
            }
        }
    };
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
        let mut m = matrix! {
            2, 0, 3;
            1, -2, 3;
            0, 1, -2;
            2, 1, 1
        };

        m.set_output(f);

        m.gauss_solve();

        Ok(())
    }
}
