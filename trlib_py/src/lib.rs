use pyo3::prelude::*;

pub mod decomposition;

#[pymodule(name = "decomposition")]
mod decomposition_py {
    #[pymodule_export]
    use crate::decomposition::PyMotifSequenceDecomposition;

    #[pymodule_export]
    use crate::decomposition::PyMotifSequenceDecomposer;
}

#[pymodule]
mod trlib_py {
    #[pymodule_export]
    use super::decomposition_py;
}
