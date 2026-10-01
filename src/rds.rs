//! Minimal writer for R's serialization format (version 3, XDR), as produced
//! by `saveRDS()` (gzip-compressed). Only the object types needed for the
//! alldata tibble are implemented.

use crate::signal::Mapping;
use std::collections::HashMap;
use std::io::{self, Write};

const NILVALUE_SXP: i32 = 254;
const REFSXP: i32 = 255;
const SYMSXP: i32 = 1;
const LISTSXP: i32 = 2;
const CHARSXP: i32 = 9;
const INTSXP: i32 = 13;
const REALSXP: i32 = 14;
const STRSXP: i32 = 16;
const VECSXP: i32 = 19;
const IS_OBJECT: i32 = 1 << 8;
const HAS_ATTR: i32 = 1 << 9;
const HAS_TAG: i32 = 1 << 10;
const ASCII_MASK: i32 = 1 << 6;
const UTF8_MASK: i32 = 1 << 3;
const NA_INTEGER: i32 = i32::MIN;

pub struct RdsWriter<W: Write> {
    w: W,
    symbols: HashMap<&'static str, i32>,
}

impl<W: Write> RdsWriter<W> {
    pub fn new(mut w: W) -> io::Result<Self> {
        w.write_all(b"X\n")?;
        let mut s = RdsWriter { w, symbols: HashMap::new() };
        s.int(3)?; // serialization version
        s.int(0x0004_0403)?; // written by R 4.4.3
        s.int(0x0003_0500)?; // readable by R >= 3.5.0
        s.int(5)?;
        s.w.write_all(b"UTF-8")?;
        Ok(s)
    }

    pub fn finish(self) -> W {
        self.w
    }

    #[inline]
    fn int(&mut self, v: i32) -> io::Result<()> {
        self.w.write_all(&v.to_be_bytes())
    }

    fn len(&mut self, n: usize) -> io::Result<()> {
        if n > i32::MAX as usize {
            self.int(-1)?;
            self.int((n >> 32) as i32)?;
            self.int(n as u32 as i32)
        } else {
            self.int(n as i32)
        }
    }

    fn charsxp(&mut self, s: &str) -> io::Result<()> {
        let gp = if s.is_ascii() { ASCII_MASK } else { UTF8_MASK };
        self.int(CHARSXP | (gp << 12))?;
        self.int(s.len() as i32)?;
        self.w.write_all(s.as_bytes())
    }

    fn symbol(&mut self, name: &'static str) -> io::Result<()> {
        if let Some(&idx) = self.symbols.get(name) {
            return self.int((idx << 8) | REFSXP);
        }
        let idx = self.symbols.len() as i32 + 1;
        self.symbols.insert(name, idx);
        self.int(SYMSXP)?;
        self.charsxp(name)
    }

    fn strsxp<S: AsRef<str>>(&mut self, v: &[S]) -> io::Result<()> {
        self.int(STRSXP)?;
        self.len(v.len())?;
        for s in v {
            self.charsxp(s.as_ref())?;
        }
        Ok(())
    }

    fn intsxp(&mut self, flags: i32, v: impl ExactSizeIterator<Item = i32>) -> io::Result<()> {
        self.int(INTSXP | flags)?;
        self.len(v.len())?;
        for x in v {
            self.int(x)?;
        }
        Ok(())
    }

    fn realsxp(&mut self, v: impl ExactSizeIterator<Item = f64>) -> io::Result<()> {
        self.int(REALSXP)?;
        self.len(v.len())?;
        for x in v {
            self.w.write_all(&x.to_bits().to_be_bytes())?;
        }
        Ok(())
    }

    fn attr_tag(&mut self, name: &'static str) -> io::Result<()> {
        self.int(LISTSXP | HAS_TAG)?;
        self.symbol(name)
    }

    fn end_attrs(&mut self) -> io::Result<()> {
        self.int(NILVALUE_SXP)
    }

    /// Attributes of a tibble: row.names, names, class.
    fn tibble_attrs(&mut self, names: &[&str], nrow: usize, compact_rownames: bool) -> io::Result<()> {
        self.attr_tag("row.names")?;
        if compact_rownames {
            self.intsxp(0, [NA_INTEGER, -(nrow as i32)].into_iter())?;
        } else {
            self.intsxp(0, 1..nrow as i32 + 1)?;
        }
        self.attr_tag("names")?;
        self.strsxp(names)?;
        self.attr_tag("class")?;
        self.strsxp(&["tbl_df", "tbl", "data.frame"])?;
        self.end_attrs()
    }

    fn factor(&mut self, codes: impl ExactSizeIterator<Item = i32>, levels: &[String]) -> io::Result<()> {
        self.intsxp(IS_OBJECT | HAS_ATTR, codes)?;
        self.attr_tag("levels")?;
        self.strsxp(levels)?;
        self.attr_tag("class")?;
        self.strsxp(&["factor"])?;
        self.end_attrs()
    }

    /// The step-04 `alldata` tibble. With one modification X (v18_Br / v18_E):
    /// read_id, flag, chrom, strand, start, end, signalbin(positions, signalX), med_signal, med_signalbin.
    /// With two modifications X, Y (v18_BrE / DS421_v2):
    /// ..., signalbin(positions, signalX, signalY), med_signalXbin, med_signalYbin.
    pub fn write_alldata(&mut self, maps: &[Mapping], chrom_levels: &[String], mods: &[String]) -> io::Result<()> {
        let n = maps.len();
        let nm = mods.len();
        let mut cols: Vec<String> =
            ["read_id", "flag", "chrom", "strand", "start", "end", "signalbin"].iter().map(|s| s.to_string()).collect();
        if nm == 1 {
            cols.push("med_signal".into());
            cols.push("med_signalbin".into());
        } else {
            cols.extend(mods.iter().map(|m| format!("med_signal{m}bin")));
        }
        let mut bin_cols = vec!["positions".to_string()];
        bin_cols.extend(mods.iter().map(|m| format!("signal{m}")));
        let bin_cols: Vec<&str> = bin_cols.iter().map(|s| s.as_str()).collect();

        self.int(VECSXP | IS_OBJECT | HAS_ATTR)?;
        self.len(cols.len())?;
        // read_id
        self.int(STRSXP)?;
        self.len(n)?;
        for m in maps {
            self.charsxp(&m.read_id)?;
        }
        self.intsxp(0, maps.iter().map(|m| m.flag as i32))?;
        self.factor(maps.iter().map(|m| m.chrom as i32 + 1), chrom_levels)?;
        let strand_levels = ["+".to_string(), "-".to_string(), "*".to_string()];
        self.factor(maps.iter().map(|m| if m.minus { 2 } else { 1 }), &strand_levels)?;
        self.realsxp(maps.iter().map(|m| m.start as f64))?;
        self.realsxp(maps.iter().map(|m| m.end as f64))?;
        // signalbin: list of tibbles (positions, signalX[, signalY])
        self.int(VECSXP)?;
        self.len(n)?;
        for m in maps {
            self.int(VECSXP | IS_OBJECT | HAS_ATTR)?;
            self.len(1 + nm)?;
            self.realsxp(m.bins.iter().map(|b| b.0))?;
            for k in 0..nm {
                self.realsxp(m.bins.iter().map(|b| b.1[k]))?;
            }
            self.tibble_attrs(&bin_cols, m.bins.len(), false)?;
        }
        if nm == 1 {
            self.realsxp(maps.iter().map(|m| m.med_signal))?;
        }
        for k in 0..nm {
            self.realsxp(maps.iter().map(|m| m.med_signalbin[k]))?;
        }
        let cols: Vec<&str> = cols.iter().map(|s| s.as_str()).collect();
        self.tibble_attrs(&cols, n, true)
    }
}
