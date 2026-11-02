from ugme import UGME
import numpy as np


def test_UGME_U():
    """
    Dummy example for static TPM.

    Resulting U(t) and TPMs should be dim(100,4,4). Each 4x4 TPM/U(t) should be np.eye(4).
    """
    eye = np.eye(4)
    test_tpm = np.array([eye for _ in range(100)])
    model = UGME(dt=0.1).fit_fetch(test_tpm)
    assert model.check_initial_U()
    assert np.allclose(model.tpm, model.reconstruct_tpm())
    TPMs = model.tpm
    Ut = model.U
    assert np.allclose(test_tpm, TPMs), "UGME procedure has failed. Error in computed TPM is too large."
    assert np.allclose(test_tpm, Ut), "UGME procedure has failed. Error in computed U(t) is too large."
