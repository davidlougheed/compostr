use pyo3::prelude::*;

#[pymodule(name = "compostr")]
mod compostr_py {
    #[pymodule_export]
    use super::decomposition;
    #[pymodule_export]
    use super::encoding;
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
        pub fn cigar(&self) -> PyResult<Vec<trlib::decomposition::CigarItem>> {
            self.inner
                .cigar()
                .map_err(|e| PyErr::new::<PyException, _>(PyException::new_err(e.to_string())))
        }

        /// Returns a string representation of .cigar()
        fn cigar_string(&self) -> PyResult<String> {
            self.inner
                .cigar_string()
                .map_err(|e| PyErr::new::<PyException, _>(PyException::new_err(e.to_string())))
        }

        /// Given the original sequence, returns the decomposition of it into sub-slices.
        pub fn sequence_items<'py>(
            &self,
            py: Python<'py>,
            seq: &Bound<'py, PyBytes>,
            with_unmapped: bool,
        ) -> PyResult<Bound<'py, PyList>> {
            let res = self.inner.sequence_items(seq.as_bytes(), with_unmapped);
            PyList::new(py, res.into_iter().map(|r| PyBytes::new(py, r)))
        }

        // -------------------------------------------------------------------------------------------------------------

        fn to_json(&self, pretty: bool) -> PyResult<String> {
            (if pretty {
                serde_json::to_string_pretty(&self.inner)
            } else {
                serde_json::to_string(&self.inner)
            })
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
            compute_cigars: bool,
            motif_alignment_score_cutoff: Option<i32>,
        ) -> PyResult<Self> {
            Ok(PyMotifSequenceDecomposer {
                inner: {
                    trlib::decomposition::MotifSequenceDecomposer::new(
                        trlib::motif::MotifSet::new(motifs),
                        match_score,
                        mismatch_score,
                        gap_penalty,
                        compute_cigars,
                        motif_alignment_score_cutoff,
                    )
                    .map_err(|e| PyException::new_err(e.to_string()))?
                },
            })
        }

        pub fn decompose<'py>(
            &self,
            py: Python<'py>,
            seq: &Bound<'py, PyBytes>,
        ) -> PyResult<Bound<'py, PyMotifSequenceDecomposition>> {
            Bound::new(
                py,
                PyMotifSequenceDecomposition {
                    inner: self
                        .inner
                        .decompose(seq.as_bytes())
                        .map_err(|e| PyException::new_err(e.to_string()))?,
                },
            )
        }
    }
}

#[pymodule(submodule, name = "encoding")]
mod encoding {
    use pyo3::exceptions::PyException;
    use pyo3::prelude::*;

    #[pyclass(name = "EncodeAlphabet")]
    pub struct PyEncodeAlphabet {
        inner: trlib::encoding::EncodeAlphabet,
    }

    #[pymethods]
    impl PyEncodeAlphabet {}

    #[pyclass(name = "Encoder")]
    pub struct PyEncoder {
        inner: trlib::encoding::Encoder,
    }

    #[pymethods]
    impl PyEncoder {
        #[new]
        fn py_new(alphabet: &str) -> PyResult<Self> {
            Ok(PyEncoder {
                inner: trlib::encoding::Encoder::new_from_str(alphabet),
            })
        }

        fn encode(&mut self, motifs: Vec<Vec<u8>>) -> PyResult<String> {
            let inner_res = self
                .inner
                .encode_vec_u8(motifs)
                .map_err(|e| PyException::new_err(e.to_string()))?;
            Ok(inner_res.0.iter().cloned().collect())
        }

        fn get_alphabet_chars(&self) -> PyResult<String> {
            let inner_res = self
                .inner
                .get_alphabet_chars()
                .map_err(|e| PyException::new_err(e.to_string()))?;
            Ok(inner_res.iter().cloned().collect())
        }
    }
}
