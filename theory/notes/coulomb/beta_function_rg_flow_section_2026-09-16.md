# §X. The Tolmachev logarithm as a rank-one RG flow

*Draft section for P4-D, prepared 2026-09-16. Insertion point: after §1 of
`P4D_coulomb_operator.md` (operator form $\mathfrak K_c=2|u_c\rangle\langle u_c|$) and
before/alongside §3 of `P4D_coulomb_uniqueness.md` (rank-one resolvent reduction). Notation
follows those two files: $\mu^*$ is the bare Morel–Anderson pseudopotential, $\tilde\mu:=\mu^*/\lambda$,
$\varpi$ the Einstein reduced temperature, $N_c=\lfloor\omega_c/2\pi T\rfloor$.*

**Strictness map for this section**: §X.1–X.2 are **analytic** (closed-form identities and
algebra, no numerics). §X.3 is a **leading-log matching argument**, standard in the
Morel–Anderson/Tolmachev/Bogoliubov–Tolmachev–Shirkov literature; it is presented as a
heuristic derivation of the functional form, not a proof that the log-coefficient is exactly
$1$ from the operator alone. §X.4 (Gate C) is **numerical-evidence** (float64, finite-$N$,
no interval arithmetic) — it is the thing that pins the coefficient to exactly $-1$ and rules
out the alternative $-2$. §X.5 (rank-$k$ generalization) is a **scaling remark**, not a
theorem: no multi-channel Coulomb operator has been constructed or certified anywhere in this
repo (cf. `P4D_coulomb_uniqueness.md` (H6), which flags that even the single-atom rank-one
case has not been extended to multi-atom $P$).

---

## X.1 The rank-one norm is an exact log (analytic)

By `P4D_coulomb_uniqueness.md`, event 1.2 (two-line proof, reproduced here for
self-containedness): $\mathfrak K_c(N_c) = P_{N_c}(\mathcal O[1]+I)P_{N_c} = 2\,u_c u_c^\top$,
$u_c[n] = (2n+1)^{-1/2}\mathbb 1_{n<N_c}$. Being rank one, its only nonzero eigenvalue is its
trace,

$$
\|\mathfrak K_c(N_c)\| \;=\; 2\,\|u_c\|^2 \;=\; 2\,S(N_c),\qquad
S(N_c):=\sum_{n=0}^{N_c-1}\frac{1}{2n+1}. \tag{X.1}
$$

$S(N_c)$ has the exact closed form $S(N_c)=\tfrac12\big[\psi(N_c+\tfrac12)-\psi(\tfrac12)\big]$
($\psi$ = digamma), verified against the direct sum to machine precision at $N_c\in\{10,50,200\}$.
Using $\psi(\tfrac12)=-\gamma_E-2\ln2$ and the standard half-integer asymptotic
$\psi(N+\tfrac12)=\ln N+O(N^{-2})$ (no $O(N^{-1})$ term — a special feature of half-integer
digamma), (X.1) gives

$$
\boxed{\;\|\mathfrak K_c(N_c)\| \;=\; \ln N_c \;+\; (\gamma_E+\ln4) \;+\; O(N_c^{-2})\;}
\qquad \gamma_E+\ln4 = 1.963510026\ldots \tag{X.2}
$$

(numerically checked: the residual $2S(N)-\ln N-1.963510\ldots$ is $4.16\times10^{-4}$,
$1.67\times10^{-5}$, $1.04\times10^{-6}$ at $N=10,50,200$ — ratios $\approx 25=(50/10)^2$,
confirming the $O(N^{-2})$ rate rather than $O(N^{-1})$). This reproduces the constant already
recorded in `P4D_coulomb_operator.md` §1 ("$2S=\ln N_c+1.9635$"); (X.2) supplies the closed
form and the rate.

Combined with $N_c=\lfloor\omega_c/2\pi T\rfloor$, at **fixed** $T$,

$$
\|\mathfrak K_c(\omega_c)\| \;=\; \ln\!\Big(\frac{\omega_c}{2\pi T}\Big) + 1.963510\ldots + O(N_c^{-2}).
\tag{X.3}
$$

(X.3) *is* the Tolmachev logarithm: a single rank-one operator whose only spectral weight is a
logarithm of the window size $\omega_c/T$. No sum over channels, no assumption about the
phonon kernel — (X.1)–(X.3) hold for any $\varpi$, any $\lambda$, because $\mathfrak K_c$
does not depend on the phonon spectral form.

---

## X.2 Exact secular equation for the pencil (analytic)

The Tc condition is that the top eigenvalue of the pencil $\lambda H(\varpi) - \mu^*\mathfrak K_c(N_c)$
equal $1$ (`P4D_coulomb_operator.md` §1). Because $\mathfrak K_c=2|u_c\rangle\langle u_c|$ has
rank one, this is a genuine secular equation, not an approximation — this is the content of
Lemma 1 of `P4D_coulomb_uniqueness.md` §3, and is also what `p4d_coulomb_certify.py`
(`mustar_fixed_T`, lines 914–919) solves directly:

$$
1 \;=\; 2\mu^*\, u_c^\top(\lambda H-I)^{-1}u_c \;=:\; 2\mu^*\,R(\varpi,\lambda,N_c). \tag{X.4}
$$

(X.4) is exact for every $N_c$ at which $1\notin\operatorname{spec}(\lambda H)$ and the root
lands above $\mu_2(\lambda H)$ (the regime certified in `P4D_coulomb_uniqueness.md` §3's
theorem). It is the operator-level statement of "one linear solve per window" that Gate C's own
header comment (`p4d_coulomb_certify.py:878–881`) already uses to avoid bisection.

First-order (weak-$\mu^*$) expansion of (X.4) in the top eigenvector $\psi_1$ of $H(\varpi)$
($\mu_1$ its eigenvalue, $c_1:=\langle u_c,\psi_1\rangle$, from `P4D_coulomb_uniqueness.md`
Lemma 1(iii): $\mu_1-\lambda\le 2\mu^* c_1^2$, saturated at leading order) gives

$$
\mu^*(N_c) \;\approx\; \frac{\lambda\mu_1(\varpi)-1}{2\,c_1(N_c)^2}
\qquad\text{(leading order in }\mu^*\text{).} \tag{X.5}
$$

(X.5) makes explicit that all of the $N_c$-dependence of $\mu^*$ at fixed $Tc$ (fixed
$\lambda,\varpi$, hence fixed $\mu_1$) flows through $c_1(N_c)^2=|\langle u_c(N_c),\psi_1\rangle|^2$.
This is the honest analytic content: the log growth of $1/\mu^*$ is a statement about how fast
the rank-one window $u_c(N_c)$ picks up overlap with the phonon operator's own top eigenvector
as the window widens — it is *not* simply a restatement of (X.2), though the two turn out to
track each other (see X.3).

---

## X.3 Leading-log matching (heuristic derivation of the slope, standard in the MA/Tolmachev literature)

Equation (X.5) alone does not fix the *rate* at which $c_1(N_c)^2$ grows with $N_c$ — that
depends on the profile of $\psi_1[n]$ at large $n$, which is not separately certified in this
repo. The standard physics argument (Morel–Anderson 1962; Tolmachev; textbook shell-elimination)
is the following leading-log ansatz, which we adopt here as a *derivation of the functional
form*, not a proof of the coefficient:

- In the window-dominated (weak-coupling-in-$\mu^*$, large-$N_c$) regime, $\psi_1[n]$ tracks the
  same $1/\sqrt{2n+1}$ envelope as $u_n$ over the range $n<N_c$ up to a slowly-varying
  ($N_c$-independent) prefactor $\kappa(\varpi,\lambda)$. Then
  $c_1(N_c) \approx \kappa\sum_{n<N_c}u_n^2 = \kappa\,S(N_c)$, and by (X.2),
  $c_1(N_c)^2 \sim \kappa^2\big[\tfrac12\ln N_c+O(1)\big]^2$ grows *logarithmically in its
  square root* — i.e. $\ln(1/\mu^*)$ via (X.5) grows **linearly** in $\ln N_c$, matching the
  generic ansatz used below. (This step is exactly the "smuggled" high-energy-shell
  elimination that KAY I p.9 footnote 4 objects to on first-principles grounds — see
  `P4D_coulomb_operator.md` §0; we do not re-litigate that objection here, only its algebra.)
- Writing the resulting self-consistent relation in its standard closed form,

$$
\frac{1}{\mu^*(N_c)} \;=\; -a\,\ln N_c + b(\varpi,\lambda), \tag{X.6}
$$

  the coefficient $a$ is fixed by the **normalization built into $\mathfrak K_c$** — the prefactor
  $2$ in (X.1) and the exponent $1/2$ in $u_n=(2n+1)^{-1/2}$ — not by $\varpi,\lambda$. Gate
  C's own header comment states this precisely: *"prefactor 2 in $2|u\rangle\langle u|$ ... vs.
  prefactor 1 → slope $-2$"* (`p4d_coulomb_certify.py:876–877`). The rank-one, prefactor-2
  construction gives $a=1$; this is the content of §X.4.

Equation (X.6) is the Morel–Anderson formula $\mu^*(\omega_c)=\mu/[1+\mu\ln(E_F/\omega_c)]$ in
disguise (invert (X.6): $\mu^*=1/(b-a\ln N_c)$, same one-parameter family), specialized to the
rank-one operator's own normalization.

---

## X.4 Certificate/numerical-evidence: Gate C fixes $a=1$

`p4d_coulomb_certify.py`, Stage G, item G2 (`run_stage_G`, lines 949–1013) computes, at fixed
$\varpi\in\{2,6\}$ and $\lambda$ calibrated so that $\mu^*(N_{c,\mathrm{ref}})=0.100$, the
sequence $\mu^*(N_c)$ for $N_c\in\{16,\ldots,512\}$ (resp. $\{32,\ldots,512\}$) via the exact
secular equation (X.4) (`mustar_fixed_T`), and fits the slope of $1/\mu^*$ against $\ln N_c$:

$$
\left.\frac{d(1/\mu^*)}{d\ln N_c}\right|_{\text{Gate C, float64}} = -0.99955\ \text{to}\ -0.99735
\quad\text{(largest windows, both $\varpi$ values), monotonically approaching } -1.
$$

The alternative normalization ($2|u\rangle\langle u|$ replaced by $|u\rangle\langle u|$,
prefactor 1) is checked and excluded: it gives slope $\approx-2$, not $-1$
(`p4d_coulomb_certify.py:1010-1011`, "the prefactor-1 alternative (slope −2) is EXCLUDED").
This confirms $a=1$ in (X.6) to $<10^{-3}$ at **numerical-evidence** level — a float64 point
computation on finite truncations $N\le 1.5\times4\times512$, not an interval-arithmetic bound.
No certified (mpmath/`iv`) version of Gate C exists in this repo as of this writing; upgrading
it to certificate level (interval-enclosing $R(\varpi,\lambda,N_c)$ via (X.4) and Richardson-in-$N$
on the slope fit) is future work, not yet done.

With $a=1$ pinned down, (X.6) reads, at fixed $T$ and using $N_c=\lfloor\omega_c/2\pi T\rfloor$,

$$
\frac{1}{\mu^*(\omega_c)} \;=\; -\ln\omega_c + C(T), \tag{X.7}
$$

i.e. $\mu^*$ is an **increasing** function of $\omega_c$ at fixed $T$ — widening the Coulomb
window strengthens the effective pseudopotential, as expected physically (less retardation
suppression as $\omega_c\to E_F$).

---

## X.5 The $\beta$-function

Differentiate (X.7) directly:

$$
\frac{d\mu^*}{d\ln\omega_c} \;=\; +(\mu^*)^2 \qquad(\text{fixed }T). \tag{X.8}
$$

To land on the Kondo/BCS-log/KT-transition **normal form** $\beta(g)=-g^2$, use the
Wilsonian convention that the RG "time" runs as the cutoff is **lowered** from a reference
scale $\Lambda$ (the same convention under which the Coulomb pseudopotential is usually said to
be "integrated down" from $E_F$ to the phonon scale): $t:=\ln(\Lambda/\omega_c)$, increasing as
$\omega_c$ decreases. Then $d/dt=-d/d\ln\omega_c$, and (X.8) becomes

$$
\boxed{\;\beta(\mu^*) \;:=\; \frac{d\mu^*}{dt} \;=\; -(\mu^*)^2\;}, \qquad
t=\ln(\Lambda/\omega_c). \tag{X.9}
$$

This is exactly the single-coupling **marginally irrelevant** flow that also governs the 1-loop
Kondo problem ($\beta(g)=-g^2+O(g^3)$, giving the Kondo scale $T_K\sim e^{-1/g}$) and the
vortex-fugacity equation along the KT critical trajectory ($dy/d\ell=y^2$ in the convention
where $\ell$ increases toward the IR, giving the essential-singularity correlation length
$\xi\sim e^{b/\sqrt{T-T_{KT}}}$). Integrating (X.9) from a reference point
$(\mu^*_0,t_0)$ reproduces the familiar resummed form

$$
\mu^*(t) \;=\; \frac{\mu^*_0}{1+\mu^*_0\,(t-t_0)}, \tag{X.10}
$$

the Morel–Anderson renormalized pseudopotential, now derived (at the level of X.3's leading-log
ansatz, pinned by Gate C) from the rank-one operator $\mathfrak K_c=2|u_c\rangle\langle u_c|$
rather than assumed.

**Physical picture.** $\mathfrak K_c$ is a single direction in Hilbert space — one projector,
not a channel sum. Its "size" as seen by the pencil is not a fixed number but the log-divergent
quantity $\|\mathfrak K_c(N_c)\|=2S(N_c)$: each e-fold of the Matsubara window
$N_c=\omega_c/2\pi T$ adds a fixed increment $\Delta(2S)=1$ to that size (from (X.2),
$d(2S)/d\ln N_c\to1$). Because the operator is rank one, this single geometric increment is
the *entire* content of the $\mu^*$-flow: there is no second direction for a compensating or
competing term to appear, so the flow is forced into the single first-order ODE (X.9). This is
the operator-theoretic reading of "the rank-one projection sweeping out one unit of logarithmic
phase per e-fold of frequency window."

---

## X.6 Rank-$k$ generalization (scaling remark, not a theorem)

If the Coulomb operator were instead a sum of $k$ mutually near-orthogonal rank-one pieces,
$\mathfrak K_c^{(k)}=2\sum_{i=1}^k|u_c^{(i)}\rangle\langle u_c^{(i)}|$, each with its own window
$N_c^{(i)}\sim N_c$, the trace is additive:

$$
\operatorname{tr}\mathfrak K_c^{(k)}(N_c) \;=\; \sum_{i=1}^k \|\mathfrak K_c^{(i)}(N_c)\|
\;\sim\; k\ln N_c + O(1). \tag{X.11}
$$

(X.11) is exact and analytic (trace is linear; each summand obeys (X.2)). If the same
leading-log matching of §X.3 applies term-by-term — which requires the $k$ pieces to stay
mutually near-orthogonal and each to individually saturate its own secular contribution, neither
of which is established for any concrete multi-channel construction — the coefficient in (X.6)
becomes $a=k$, giving

$$
\beta_k(\mu^*) \;=\; -k\,(\mu^*)^2. \tag{X.12}
$$

**This is illustrative, not proven.** `P4D_coulomb_uniqueness.md` (H6) already records that no
multi-atom (let alone multi-channel-Coulomb) construction has been checked against the rank-one
machinery of §3 there; the same caveat applies here. The only rank actually realized and
certified anywhere in this repo is $k=1$ (single Morel–Anderson $\mu^*$, single cutoff $\omega_c$),
which is why Gate C measures slope $-1$ and not some other integer.

---

## X.7 Relation to existing material (no duplication)

`P4D_coulomb_operator.md` §1 already records the numeric value $\|\mathfrak K_c\|=2S=\ln
N_c+1.9635$ and cites Gate C's slope test as the "main verification" locking the prefactor 2 and
the $1/\sqrt{2n+1}$ weight. This section adds: (i) the closed-form digamma derivation and exact
$O(N_c^{-2})$ rate for (X.2) (not previously written down); (ii) the explicit secular-equation
route (X.4)–(X.5) connecting the log to the rank-one resolvent machinery of
`P4D_coulomb_uniqueness.md` §3, rather than treating the log as a bare numerical fact; (iii) the
$\beta$-function formulation (X.9) and its identification with the Kondo/KT normal form,
which is new content not present in either upstream file; (iv) the honest rank-$k$ scaling
remark (§X.6), flagged as unproven, addressing the natural "why exactly $-1$" question at the
structural level beyond what Gate C's numerics alone can answer.

**Not claimed**: that (X.6)'s coefficient $a=1$ is derived from first principles without the
leading-log ansatz of §X.3 — a fully rigorous derivation of $a$ from the operator alone (without
assuming $\psi_1[n]$ tracks $u_n$'s envelope) is open, and is a natural target for upgrading
this section from "leading-log heuristic + numerical-evidence confirmation" to "analytic."
