#!/usr/bin/env python3
"""
Graser et al. arXiv:1003.0133 5-orbital tight-binding model for BaFe2As2
Real DFT-fit parameters (Appendix A, Table I/II), 2D cut at kz=0.

Orbital order (paper's own convention, NOT the common dxy=1 convention):
  0 = dxz, 1 = dyz, 2 = dx2-y2, 3 = dxy, 4 = d3z2-r2

Source: theory/notes/_scratch_ir/graser1003.0133.txt lines 545-626
Note: mn=35 and mn=45 rows in Table II are blank in the paper -> xi35=xi45=0 (not an omission).
kz=0 cut kills all sin(kz) terms and sets cos(kz)=1; this is an explicit simplification
(BaFe2As2 is a 3D compound) documented here, not silently assumed.
"""

import numpy as np
from numpy.linalg import eigh, inv

# ---------------------------------------------------------------------------
# On-site energies (eV), measured from E_F, paper text L.578-580
# ---------------------------------------------------------------------------
eps = np.array([0.0987, 0.0987, -0.3595, 0.2078, -0.7516])  # 1,2,3,4,5 -> idx 0..4

# ---------------------------------------------------------------------------
# Table I: intraorbital hopping parameters (m=1,3,4,5; m=2 by C4 symmetry x<->y)
# ---------------------------------------------------------------------------
t1 = dict(x=-0.0604, y=-0.3005, xx=0.0253, xy=0.2388, xxy=-0.0414,
          xyy=-0.0237, xxyy=0.0158, z=0.0, xz=-0.0101, xxz=0.0126, xyz=0.0421)
t3 = dict(x=0.3378, xy=0.0011, xxy=-0.0947, z=0.1001)  # only nonzero entries per table
# reinterpret cleanly against the printed xi33 formula (x, xy(=cross), xx term):
# xi33 = 2 t3_x (coskx+cosky) + 4 t3_xy_cross coskx cosky + 2 t3_xx (cos2kx+cos2ky)
t3 = dict(x=0.3378, xx=0.0011, xy=-0.0947, z=0.1001)
t4 = dict(x=0.1965, xx=-0.0528, xy=0.1259, xxy=-0.032, xxyy=0.0045, z=0.0563, xz=0.0662, xyz=0.0)
t5 = dict(x=-0.0656, xx=0.0001, xxy=0.01, xxyy=0.0047, z=0.0, xz=-0.0036)

# ---------------------------------------------------------------------------
# Table II: interorbital hopping parameters
# ---------------------------------------------------------------------------
t12 = dict(x=-0.4224, xy=0.1934, xxy=-0.0325, xxyy=0.0158, z=-0.019, xz=0.0524, xyz=-0.0168, xxyz=0.0)
t13 = dict(x=0.1549, xy=0.0589, xxy=0.0005, xxyy=0.0, z=0.0, xz=0.0566, xyz=0.0, xxyz=0.0018)
t14 = dict(x=-0.0526, xy=-0.007, xxy=-0.0055, xxyy=0.0004, z=0.0, xz=-0.0023, xyz=0.0349, xxyz=0.0283)
t24 = dict(x=-0.2845, xy=0.0, xxy=-0.0108)  # only feeds real (sin kz) terms -> vanish at kz=0
t15 = dict(x=0.0, xy=-0.0862, xxy=0.0, xxyy=0.0, z=0.0, xz=0.0, xyz=-0.0203, xxyz=0.0)
t34 = dict(xy=-0.0475, xxy=0.0046)
# t35, t45 rows blank in Table II -> identically zero at all k (documented, not a bug)


def _xi11(kx, ky, kz=0.0):
    """xi_11 (dxz); xi_22 obtained by swapping kx<->ky (C4 symmetry partner)."""
    return (2*t1['x']*np.cos(kx) + 2*t1['y']*np.cos(ky)
            + 4*t1['xx']*np.cos(kx)*np.cos(ky)
            + 2*t1['xy']*(np.cos(2*kx) - np.cos(2*ky))
            + 4*t1['xxy']*np.cos(2*kx)*np.cos(ky)
            + 4*t1['xyy']*np.cos(2*ky)*np.cos(kx)
            + 4*t1['xxyy']*np.cos(2*kx)*np.cos(2*ky)
            + 4*t1['xz']*(np.cos(kx)+np.cos(ky))*np.cos(kz)
            + 4*t1['xxz']*(np.cos(2*kx)-np.cos(2*ky))*np.cos(kz))


def _xi33(kx, ky, kz=0.0):
    return (2*t3['x']*(np.cos(kx)+np.cos(ky))
            + 4*t3['xy']*np.cos(kx)*np.cos(ky)
            + 2*t3['xx']*(np.cos(2*kx)+np.cos(2*ky)))


def _xi44(kx, ky, kz=0.0):
    return (2*t4['x']*(np.cos(kx)+np.cos(ky))
            + 4*t4['xy']*np.cos(kx)*np.cos(ky)
            + 2*t4['xx']*(np.cos(2*kx)+np.cos(2*ky))
            + 4*t4['xxy']*(np.cos(2*kx)*np.cos(ky) + np.cos(2*ky)*np.cos(kx))
            + 4*t4['xxyy']*np.cos(2*kx)*np.cos(2*ky)
            + 2*t4['z']*np.cos(kz)
            + 4*t4['xz']*(np.cos(kx)+np.cos(ky))*np.cos(kz)
            + 8*t4['xyz']*np.cos(kx)*np.cos(ky)*np.cos(kz))


def _xi55(kx, ky, kz=0.0):
    return (2*t5['x']*(np.cos(kx)+np.cos(ky))
            + 2*t5['xx']*(np.cos(2*kx)+np.cos(2*ky))
            + 4*t5['xxy']*(np.cos(2*kx)*np.cos(ky) + np.cos(2*ky)*np.cos(kx))
            + 4*t5['xxyy']*np.cos(2*kx)*np.cos(2*ky)
            + 2*t5['z']*np.cos(kz)
            + 4*t5['xz']*(np.cos(kx)+np.cos(ky))*np.cos(kz))


def _xi12(kx, ky, kz=0.0):
    return (4*t12['x']*np.sin(kx)*np.sin(ky)
            + 4*t12['xy']*(np.sin(2*kx)*np.sin(ky) + np.sin(2*ky)*np.sin(kx))
            + 4*t12['xxyy']*np.sin(2*kx)*np.sin(2*ky)
            + 8*t12['xyz']*np.sin(kx)*np.sin(ky)*np.cos(kz))


def _xi13(kx, ky, kz=0.0):
    return 1j*(2*t13['x']*np.sin(ky) + 4*t13['xy']*np.sin(ky)*np.cos(kx)
               - 4*t13['xxy']*(np.sin(2*ky)*np.cos(kx) - np.cos(2*kx)*np.sin(ky)))


def _xi23(kx, ky, kz=0.0):
    return 1j*(2*t13['x']*np.sin(kx) + 4*t13['xy']*np.sin(kx)*np.cos(ky)
               - 4*t13['xxy']*(np.sin(2*kx)*np.cos(ky) - np.cos(2*ky)*np.sin(kx)))


def _xi14(kx, ky, kz=0.0):
    # imaginary (t14) part only; real (t24, sin kz) part vanishes at kz=0
    return 1j*(2*t14['x']*np.sin(kx) + 4*t14['xy']*np.cos(ky)*np.sin(kx)
               + 4*t14['xxy']*np.sin(2*kx)*np.cos(ky)
               + 4*t14['xz']*np.sin(kx)*np.cos(kz)
               + 8*t14['xyz']*np.cos(ky)*np.sin(kx)*np.cos(kz)
               + 8*t14['xxyz']*np.sin(2*kx)*np.cos(ky)*np.cos(kz))


def _xi24(kx, ky, kz=0.0):
    return -1j*(2*t14['x']*np.sin(ky) + 4*t14['xy']*np.cos(kx)*np.sin(ky)
                + 4*t14['xxy']*np.sin(2*ky)*np.cos(kx)
                + 4*t14['xz']*np.sin(ky)*np.cos(kz)
                + 8*t14['xyz']*np.cos(kx)*np.sin(ky)*np.cos(kz)
                + 8*t14['xxyz']*np.sin(2*ky)*np.cos(kx)*np.cos(kz))


def _xi15(kx, ky, kz=0.0):
    return 1j*(2*t15['x']*np.sin(ky) - 4*t15['xy']*np.sin(ky)*np.cos(kx)
               - 8*t15['xyz']*np.sin(ky)*np.cos(kx)*np.cos(kz))


def _xi25(kx, ky, kz=0.0):
    return -1j*(2*t15['x']*np.sin(kx) - 4*t15['xy']*np.sin(kx)*np.cos(ky)
                - 8*t15['xyz']*np.sin(kx)*np.cos(ky)*np.cos(kz))


def _xi34(kx, ky, kz=0.0):
    return 4*t34['xy']*(np.sin(2*ky)*np.sin(kx) - np.sin(2*kx)*np.sin(ky))


def build_H_graser(kx, ky, kz=0.0):
    """5x5 Hermitian Hamiltonian, Graser et al. BaFe2As2 fit, 2D cut kz=0."""
    H = np.zeros((5, 5), dtype=complex)
    H[0, 0] = eps[0] + _xi11(kx, ky, kz)
    H[1, 1] = eps[1] + _xi11(ky, kx, kz)     # C4 partner: swap kx<->ky
    H[2, 2] = eps[2] + _xi33(kx, ky, kz)
    H[3, 3] = eps[3] + _xi44(kx, ky, kz)
    H[4, 4] = eps[4] + _xi55(kx, ky, kz)

    H[0, 1] = _xi12(kx, ky, kz);  H[1, 0] = np.conj(H[0, 1])
    H[0, 2] = _xi13(kx, ky, kz);  H[2, 0] = np.conj(H[0, 2])
    H[1, 2] = _xi23(kx, ky, kz);  H[2, 1] = np.conj(H[1, 2])
    H[0, 3] = _xi14(kx, ky, kz);  H[3, 0] = np.conj(H[0, 3])
    H[1, 3] = _xi24(kx, ky, kz);  H[3, 1] = np.conj(H[1, 3])
    H[0, 4] = _xi15(kx, ky, kz);  H[4, 0] = np.conj(H[0, 4])
    H[1, 4] = _xi25(kx, ky, kz);  H[4, 1] = np.conj(H[1, 4])
    H[2, 3] = _xi34(kx, ky, kz);  H[3, 2] = np.conj(H[2, 3])
    H[2, 4] = 0.0;                H[4, 2] = 0.0   # xi35 = 0 (blank row in Table II)
    H[3, 4] = 0.0;                H[4, 3] = 0.0   # xi45 = 0 (blank row in Table II)
    return H


if __name__ == "__main__":
    # sanity check: hermiticity + band structure at a few points
    for (kx, ky) in [(0, 0), (np.pi, 0), (np.pi, np.pi), (0.3, 0.7)]:
        H = build_H_graser(kx, ky)
        assert np.allclose(H, H.conj().T, atol=1e-12), f"Not Hermitian at ({kx},{ky})"
        eig = np.linalg.eigvalsh(H)
        print(f"k=({kx:.3f},{ky:.3f}): bands = {np.round(eig,4)}")
