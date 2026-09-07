use pyo3::prelude::*;

#[pymodule(name = "compostr")]
mod compostr_py {
    #[pymodule_export]
    use super::decomposition;
}

#[pymodule(submodule, name = "decomposition")]
mod decomposition {
    use pyo3::exceptions::PyException;
    use pyo3::prelude::*;
    use pyo3::types::{PyBytes, PyList};
    use serde_json;
    use trlib;

    #[pymodule_export]
    use trlib::decomposition::CigarItem;

    #[pyclass(name = "MotifSequenceDecomposition")]
    pub struct PyMotifSequenceDecomposition {
        inner: trlib::decomposition::MotifSequenceDecomposition,
    }

    #[pymethods]
    impl PyMotifSequenceDecomposition {
        /// Returns a vector of CIGAR items representing the alignment of a sequence just composed of motifs from the
        /// motif set against the original sequence.
        pub fn cigar(&self) -> Vec<trlib::decomposition::CigarItem> {
            self.inner.cigar()
        }

        /// Returns a string representation of .cigar()
        fn cigar_string(&self) -> String {
            self.inner.cigar_string()
        }

        /// Given the original sequence, returns the decomposition of it into sub-slices.
        pub fn sequence_items<'py>(
            &self, py: Python<'py>, seq: &Bound<'py, PyBytes>, with_unmapped: bool
        ) -> PyResult<Bound<'py, PyList>> {
            let res = self.inner.sequence_items(seq.as_bytes(), with_unmapped);
            PyList::new(py, res.into_iter().map(|r| PyBytes::new(py, r)))
        }

        // -------------------------------------------------------------------------------------------------------------

        fn to_json(&self, pretty: bool) -> PyResult<String> {
            (if pretty { serde_json::to_string_pretty(&self.inner) } else { serde_json::to_string(&self.inner) })
                .map_err(|err| PyErr::new::<PyException, _>(err.to_string()))
        }
    }

    #[pyclass(name = "MotifSequenceDecomposer")]
    pub struct PyMotifSequenceDecomposer {
        inner: trlib::decomposition::MotifSequenceDecomposer,
    }

    #[pymethods]
    impl PyMotifSequenceDecomposer {
        #[new]
        fn py_new(
            motifs: Vec<Vec<u8>>,
            match_score: i32,
            mismatch_score: i32,
            gap_penalty: i32,
            motif_alignment_score_cutoff: Option<i32>,
        ) -> PyResult<Self> {
            Ok(PyMotifSequenceDecomposer {
                inner:
                    trlib::decomposition::MotifSequenceDecomposer::new(
                        trlib::motif::MotifSet::new(motifs),
                        match_score,
                        mismatch_score,
                        gap_penalty,
                        motif_alignment_score_cutoff,
                    ).map_err(|e| PyException::new_err(e.to_string()))?
            })
        }

        pub fn decompose<'py>(
            &self, py: Python<'py>, seq: &Bound<'py, PyBytes>
        ) -> PyResult<Bound<'py, PyMotifSequenceDecomposition>> {
            Bound::new(py, PyMotifSequenceDecomposition {
                inner: self.inner.decompose(seq.as_bytes()).map_err(|e| PyException::new_err(e.to_string()))?
            })
        }
    }
}
