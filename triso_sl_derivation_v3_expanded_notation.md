# Analytical Solution of Fission Product Diffusion in a TRISO Particle

**Expanded-notation working**

Sturm–Liouville analytical solution for transient fission-product diffusion in a spherical TRISO particle with a constant volumetric source, zero initial inventory, and a convective (Robin) surface condition.

Written throughout in fully expanded notation:

- every derivative in Leibniz form — no prime notation;
- every intermediate algebraic step retained;
- no dimensionless groupings — the physical constants $D$, $h$, $R$, $S_0$ appear explicitly in every expression.

Covers the source-active phase, the post-shutdown phase, and the formulation of the five-layer extension.

---

## Contents

1. [Notation](#1-notation)
2. [Problem statement](#2-problem-statement)
3. [Steady–transient split](#3-steadytransient-split)
4. [Steady solution](#4-steady-solution)
5. [Transient problem](#5-transient-problem)
6. [Separation of variables](#6-separation-of-variables)
7. [Removal of the curvature term](#7-removal-of-the-curvature-term)
8. [General solution and regularity at the origin](#8-general-solution-and-regularity-at-the-origin)
9. [Eigencondition](#9-eigencondition)
10. [Positivity of the eigenvalues](#10-positivity-of-the-eigenvalues)
11. [Orthogonality](#11-orthogonality)
12. [Normalization](#12-normalization)
13. [Projection of the initial condition](#13-projection-of-the-initial-condition)
14. [Evaluation of the integrals](#14-evaluation-of-the-integrals)
15. [Expansion coefficients](#15-expansion-coefficients)
16. [Solution with the source active](#16-solution-with-the-source-active)
17. [Solution after source switch-off](#17-solution-after-source-switch-off)
18. [Extension to the five-layer particle](#18-extension-to-the-five-layer-particle)

---

## 1. Notation

Two derivative symbols appear, and the distinction between them is carried consistently.

| Function | Depends on | Symbol | Example |
|---|---|---|---|
| $c(r,t)$ | $r$ **and** $t$ | $\partial$ | $\dfrac{\partial c}{\partial t}$, $\dfrac{\partial c}{\partial r}$ |
| $v(r,t)$ | $r$ **and** $t$ | $\partial$ | $\dfrac{\partial v}{\partial t}$, $\dfrac{\partial v}{\partial r}$ |
| $w(r)$ | $r$ only | $d$ | $\dfrac{dw}{dr}$ |
| $\phi(r)$ | $r$ only | $d$ | $\dfrac{d\phi}{dr}$, $\dfrac{d^{2}\phi}{dr^{2}}$ |
| $u(r)$ | $r$ only | $d$ | $\dfrac{du}{dr}$, $\dfrac{d^{2}u}{dr^{2}}$ |
| $T(t)$ | $t$ only | $d$ | $\dfrac{dT}{dt}$ |

The symbol $\partial$ signals that other variables are being held fixed. For a function of a single variable there is nothing to hold fixed, so the ordinary symbol $d$ is used.

The change from $\partial$ to $d$ occurs in [Section 6](#6-separation-of-variables), and is precisely what separation of variables achieves: one partial differential equation is replaced by two ordinary differential equations.

The notation

$$
\frac{d\phi}{dr}\bigg|_{r=R}
$$

means: perform the differentiation first, then substitute $r=R$ into the result.

### Symbol list

| Symbol | Meaning | Units |
|---|---|---|
| $c(r,t)$ | fission-product concentration | $\mathrm{mol\ m^{-3}}$ |
| $r$ | radial coordinate | $\mathrm{m}$ |
| $t$ | time | $\mathrm{s}$ |
| $R$ | particle radius | $\mathrm{m}$ |
| $D$ | diffusivity | $\mathrm{m^{2}\ s^{-1}}$ |
| $h$ | surface mass-transfer coefficient | $\mathrm{m\ s^{-1}}$ |
| $S_0$ | volumetric source (generation rate) | $\mathrm{mol\ m^{-3}\ s^{-1}}$ |
| $t_c$ | time at which the source is switched off | $\mathrm{s}$ |
| $w(r)$ | steady-state profile | $\mathrm{mol\ m^{-3}}$ |
| $v(r,t)$ | transient remainder, $c - w$ | $\mathrm{mol\ m^{-3}}$ |
| $\lambda_n$ | eigenvalue | $\mathrm{m^{-2}}$ |
| $\phi_n(r)$ | eigenfunction | $\mathrm{m^{-1}}$ |
| $N_n$ | normalization integral | $\mathrm{m}$ |
| $a_n$, $b_n$ | expansion coefficients | $\mathrm{mol\ m^{-2}}$ |

---

## 2. Problem statement

$$
\frac{\partial c}{\partial t}
= \frac{D}{r^{2}}\,\frac{\partial}{\partial r}\!\left(r^{2}\frac{\partial c}{\partial r}\right) + S_{0},
\qquad 0 < r < R,\quad t > 0
$$

Expanding the inner derivative by the product rule gives the equivalent form

$$
\frac{\partial c}{\partial t}
= D\left(\frac{\partial^{2} c}{\partial r^{2}} + \frac{2}{r}\,\frac{\partial c}{\partial r}\right) + S_{0}
$$

The initial and boundary data are

$$
c(r,0) = 0
$$

$$
\frac{\partial c}{\partial r}\bigg|_{r=0} = 0
\qquad \text{(spherical symmetry)}
$$

$$
-D\,\frac{\partial c}{\partial r}\bigg|_{r=R} = h\,c(R,t)
\qquad \text{(convective release)}
$$

The source is active for $0 \le t \le t_{c}$ and is switched off for $t > t_{c}$.

---

## 3. Steady–transient split

$$
c(r,t) = w(r) + v(r,t)
$$

---

## 4. Steady solution

Setting the time derivative to zero leaves a function of $r$ alone, so ordinary derivatives apply from here on.

$$
\frac{D}{r^{2}}\,\frac{d}{dr}\!\left(r^{2}\frac{dw}{dr}\right) + S_{0} = 0
$$

Multiplying through by $r^{2}/D$:

$$
\frac{d}{dr}\!\left(r^{2}\frac{dw}{dr}\right) = -\frac{S_{0}\,r^{2}}{D}
$$

Integrating once with respect to $r$, with $A$ the constant of integration:

$$
r^{2}\frac{dw}{dr} = -\frac{S_{0}\,r^{3}}{3D} + A
$$

Dividing through by $r^{2}$:

$$
\frac{dw}{dr} = -\frac{S_{0}\,r}{3D} + \frac{A}{r^{2}}
$$

The term $A/r^{2}$ diverges as $r \to 0$, which would give unbounded flux at the centre and contradict the symmetry condition. Hence

$$
A = 0
\qquad\Longrightarrow\qquad
\frac{dw}{dr} = -\frac{S_{0}\,r}{3D}
$$

Integrating a second time, with $B$ the second constant of integration:

$$
w(r) = -\frac{S_{0}\,r^{2}}{6D} + B
$$

Evaluating the gradient at the surface:

$$
\frac{dw}{dr}\bigg|_{r=R} = -\frac{S_{0}R}{3D}
$$

Substituting into the surface condition:

$$
-D\,\frac{dw}{dr}\bigg|_{r=R}
= -D\left(-\frac{S_{0}R}{3D}\right)
= \frac{S_{0}R}{3}
= h\,w(R)
= h\left(B - \frac{S_{0}R^{2}}{6D}\right)
$$

Solving for $B$:

$$
B = \frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}
$$

$$
\boxed{\;
w(r) = \frac{S_{0}R}{3h} + \frac{S_{0}}{6D}\left(R^{2} - r^{2}\right)
\;}
$$

---

## 5. Transient problem

Subtracting the steady part removes the source and homogenises both boundary conditions:

$$
\frac{\partial v}{\partial t}
= \frac{D}{r^{2}}\,\frac{\partial}{\partial r}\!\left(r^{2}\frac{\partial v}{\partial r}\right)
$$

$$
\frac{\partial v}{\partial r}\bigg|_{r=0} = 0,
\qquad
-D\,\frac{\partial v}{\partial r}\bigg|_{r=R} = h\,v(R,t),
\qquad
v(r,0) = -w(r)
$$

---

## 6. Separation of variables

Seek a solution in the form of a pure product:

$$
v(r,t) = \phi(r)\,T(t)
$$

Differentiating with respect to $t$, with $r$ held fixed — the factor $\phi(r)$ contains no $t$ and behaves as a constant multiplier:

$$
\frac{\partial v}{\partial t} = \phi(r)\,\frac{dT}{dt}
$$

Differentiating with respect to $r$, with $t$ held fixed — the factor $T(t)$ contains no $r$ and behaves as a constant multiplier:

$$
\frac{\partial v}{\partial r} = \frac{d\phi}{dr}\,T(t)
$$

Multiplying by $r^{2}$:

$$
r^{2}\frac{\partial v}{\partial r} = r^{2}\frac{d\phi}{dr}\,T(t)
$$

Differentiating once more with respect to $r$, again factoring $T(t)$ out of the $r$-derivative:

$$
\frac{\partial}{\partial r}\!\left(r^{2}\frac{\partial v}{\partial r}\right)
= \frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right) T(t)
$$

> This is where the notation changes. The quantity $r^{2}\,d\phi/dr$ contains no $t$ whatsoever, so the derivative acting on it is ordinary rather than partial.

Substituting into the transient equation:

$$
\phi(r)\,\frac{dT}{dt}
= \frac{D}{r^{2}}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right) T(t)
$$

Dividing both sides by $D\,\phi(r)\,T(t)$:

$$
\frac{1}{D\,T(t)}\,\frac{dT}{dt}
= \frac{1}{r^{2}\phi(r)}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right)
$$

The left-hand side depends only on $t$; the right-hand side depends only on $r$. Since $r$ and $t$ vary independently, both sides must equal the same constant, written $-\lambda$:

$$
\frac{1}{D\,T(t)}\,\frac{dT}{dt}
= \frac{1}{r^{2}\phi(r)}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right)
= -\lambda
$$

This yields two ordinary differential equations.

**Time equation:**

$$
\frac{dT}{dt} = -\lambda D\,T(t)
\qquad \Longrightarrow \qquad
T(t) = e^{-\lambda D t}
$$

**Spatial equation:**

$$
\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right) + \lambda\,r^{2}\phi = 0
$$

---

## 7. Removal of the curvature term

Introduce

$$
u(r) = r\,\phi(r),
\qquad\text{equivalently}\qquad
\phi(r) = \frac{u(r)}{r}
$$

Differentiating $\phi = u \cdot r^{-1}$ by the product rule:

$$
\frac{d\phi}{dr}
= \frac{1}{r}\,\frac{du}{dr} + u \cdot \left(-\frac{1}{r^{2}}\right)
= \frac{1}{r}\,\frac{du}{dr} - \frac{u}{r^{2}}
$$

Multiplying through by $r^{2}$:

$$
r^{2}\frac{d\phi}{dr} = r\,\frac{du}{dr} - u
$$

Differentiating with respect to $r$, applying the product rule to the first term:

$$
\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right)
= \frac{d}{dr}\!\left(r\,\frac{du}{dr} - u\right)
= \left(\frac{du}{dr} + r\,\frac{d^{2}u}{dr^{2}}\right) - \frac{du}{dr}
= r\,\frac{d^{2}u}{dr^{2}}
$$

The remaining term of the spatial equation becomes

$$
r^{2}\phi = r^{2}\cdot\frac{u}{r} = r\,u
$$

Substituting both results:

$$
r\,\frac{d^{2}u}{dr^{2}} + \lambda\,r\,u = 0
$$

Dividing by $r$ (valid for $r > 0$):

$$
\boxed{\;
\frac{d^{2}u}{dr^{2}} + \lambda u = 0
\;}
$$

All curvature terms have cancelled.

---

## 8. General solution and regularity at the origin

Substituting the trial form $u(r) = e^{mr}$, for which $d^{2}u/dr^{2} = m^{2}e^{mr}$:

$$
m^{2} e^{mr} + \lambda e^{mr} = 0
\qquad \Longrightarrow \qquad
m^{2} + \lambda = 0
\qquad \Longrightarrow \qquad
m = \pm\, i\sqrt{\lambda}
$$

The corresponding real general solution is

$$
u(r) = C_{1}\sin\!\left(\sqrt{\lambda}\,r\right) + C_{2}\cos\!\left(\sqrt{\lambda}\,r\right)
$$

so that

$$
\phi(r)
= C_{1}\,\frac{\sin\!\left(\sqrt{\lambda}\,r\right)}{r}
+ C_{2}\,\frac{\cos\!\left(\sqrt{\lambda}\,r\right)}{r}
$$

Expanding each branch near the origin:

$$
\frac{\sin\!\left(\sqrt{\lambda}\,r\right)}{r}
= \sqrt{\lambda} - \frac{\lambda^{3/2} r^{2}}{6} + \cdots
\;\xrightarrow[\;r \to 0\;]{}\;
\sqrt{\lambda}
\qquad \text{(finite)}
$$

$$
\frac{\cos\!\left(\sqrt{\lambda}\,r\right)}{r}
= \frac{1}{r} - \frac{\lambda r}{2} + \cdots
\;\xrightarrow[\;r \to 0\;]{}\;
\infty
\qquad \text{(unbounded)}
$$

Boundedness of the concentration at the centre therefore requires

$$
C_{2} = 0
$$

Setting $C_{1} = 1$ (any overall constant is absorbed later into the expansion coefficients):

$$
\boxed{\;
u(r) = \sin\!\left(\sqrt{\lambda}\,r\right),
\qquad
\phi(r) = \frac{\sin\!\left(\sqrt{\lambda}\,r\right)}{r}
\;}
$$

The centre condition is then satisfied automatically:

$$
\phi(r) = \sqrt{\lambda} - \frac{\lambda^{3/2} r^{2}}{6} + \cdots
\qquad \Longrightarrow \qquad
\frac{d\phi}{dr} = -\frac{\lambda^{3/2} r}{3} + \cdots
\qquad \Longrightarrow \qquad
\frac{d\phi}{dr}\bigg|_{r=0} = 0
$$

---

## 9. Eigencondition

Evaluating at the surface:

$$
\phi(R) = \frac{\sin\!\left(\sqrt{\lambda}\,R\right)}{R}
$$

Differentiating $\phi(r) = \sin\!\left(\sqrt{\lambda}\,r\right) \cdot r^{-1}$ by the product rule:

$$
\frac{d\phi}{dr}
= \sqrt{\lambda}\cos\!\left(\sqrt{\lambda}\,r\right)\cdot\frac{1}{r}
- \sin\!\left(\sqrt{\lambda}\,r\right)\cdot\frac{1}{r^{2}}
$$

Combining over a common denominator:

$$
\frac{d\phi}{dr}
= \frac{\sqrt{\lambda}\,r\cos\!\left(\sqrt{\lambda}\,r\right)
- \sin\!\left(\sqrt{\lambda}\,r\right)}{r^{2}}
$$

At the surface:

$$
\frac{d\phi}{dr}\bigg|_{r=R}
= \frac{\sqrt{\lambda}\,R\cos\!\left(\sqrt{\lambda}\,R\right)
- \sin\!\left(\sqrt{\lambda}\,R\right)}{R^{2}}
$$

Substituting into the Robin condition $-D\,(d\phi/dr)|_{r=R} = h\,\phi(R)$:

$$
-D \cdot \frac{\sqrt{\lambda}\,R\cos\!\left(\sqrt{\lambda}\,R\right)
- \sin\!\left(\sqrt{\lambda}\,R\right)}{R^{2}}
= h \cdot \frac{\sin\!\left(\sqrt{\lambda}\,R\right)}{R}
$$

Multiplying both sides by $R^{2}$:

$$
-D\left(\sqrt{\lambda}\,R\cos\!\left(\sqrt{\lambda}\,R\right)
- \sin\!\left(\sqrt{\lambda}\,R\right)\right)
= hR\sin\!\left(\sqrt{\lambda}\,R\right)
$$

Distributing the minus sign:

$$
D\sin\!\left(\sqrt{\lambda}\,R\right)
- D\sqrt{\lambda}\,R\cos\!\left(\sqrt{\lambda}\,R\right)
= hR\sin\!\left(\sqrt{\lambda}\,R\right)
$$

Collecting the sine terms:

$$
\left(D - hR\right)\sin\!\left(\sqrt{\lambda}\,R\right)
= D\,\sqrt{\lambda}\,R\cos\!\left(\sqrt{\lambda}\,R\right)
$$

Dividing by $\sin\!\left(\sqrt{\lambda}\,R\right)$ gives the eigencondition:

$$
\boxed{\;
D\,\sqrt{\lambda}\,R\cot\!\left(\sqrt{\lambda}\,R\right) = D - hR
\;}
$$

This is transcendental and admits infinitely many roots:

$$
\lambda_{1} < \lambda_{2} < \lambda_{3} < \cdots,
\qquad
(n-1)\pi < \sqrt{\lambda_{n}}\,R < n\pi
$$

with corresponding eigenfunctions

$$
\phi_{n}(r) = \frac{\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{r}
$$

**Cosine-elimination identity** — used repeatedly below:

$$
\cos\!\left(\sqrt{\lambda_{n}}\,R\right)
= \frac{\left(D - hR\right)\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}
{D\,\sqrt{\lambda_{n}}\,R}
$$

---

## 10. Positivity of the eigenvalues

Multiplying the spatial equation by $\phi$ and integrating over the sphere:

$$
\int_{0}^{R} \phi\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right) dr
+ \lambda \int_{0}^{R} r^{2}\phi^{2}\,dr = 0
$$

The first integral is treated by parts. Since

$$
\frac{d}{dr}\!\left(\phi \cdot r^{2}\frac{d\phi}{dr}\right)
= \frac{d\phi}{dr}\cdot r^{2}\frac{d\phi}{dr}
+ \phi\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right)
$$

it follows that

$$
\int_{0}^{R} \phi\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi}{dr}\right) dr
= \left[r^{2}\phi\,\frac{d\phi}{dr}\right]_{0}^{R}
- \int_{0}^{R} r^{2}\left(\frac{d\phi}{dr}\right)^{2} dr
$$

Substituting and solving for $\lambda$:

$$
\lambda =
\frac{\displaystyle -\left[r^{2}\phi\,\frac{d\phi}{dr}\right]_{0}^{R}
+ \int_{0}^{R} r^{2}\left(\frac{d\phi}{dr}\right)^{2} dr}
{\displaystyle \int_{0}^{R} r^{2}\phi^{2}\,dr}
$$

At the lower limit the factor $r^{2}$ vanishes:

$$
\left[r^{2}\phi\,\frac{d\phi}{dr}\right]_{r=0} = 0
$$

At the upper limit the Robin condition gives

$$
\frac{d\phi}{dr}\bigg|_{r=R} = -\frac{h}{D}\,\phi(R)
$$

so that

$$
-\left[r^{2}\phi\,\frac{d\phi}{dr}\right]_{r=R}
= -R^{2}\phi(R)\left(-\frac{h}{D}\phi(R)\right)
= \frac{h}{D}\,R^{2}\phi(R)^{2} \;\ge\; 0
$$

Both terms in the numerator are non-negative and the denominator is strictly positive, so $\lambda \ge 0$. The value $\lambda = 0$ would force $r^{2}\,d\phi/dr$ constant, hence zero by regularity, hence $\phi$ constant; the Robin condition with $h > 0$ then forces $\phi \equiv 0$, which is not an eigenfunction. Therefore

$$
\lambda_{n} > 0 \qquad \text{for all } n
$$

and every mode $e^{-\lambda_{n} D t}$ decays.

---

## 11. Orthogonality

Let $\phi_{n}$ and $\phi_{m}$ be eigenfunctions with $\lambda_{n} \neq \lambda_{m}$:

$$
\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{n}}{dr}\right) + \lambda_{n} r^{2}\phi_{n} = 0
$$

$$
\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{m}}{dr}\right) + \lambda_{m} r^{2}\phi_{m} = 0
$$

Multiplying the first by $\phi_{m}$, the second by $\phi_{n}$, and subtracting:

$$
\phi_{m}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{n}}{dr}\right)
- \phi_{n}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{m}}{dr}\right)
= \left(\lambda_{m} - \lambda_{n}\right) r^{2}\phi_{n}\phi_{m}
$$

The left-hand side is an exact derivative, since the cross terms cancel:

$$
\frac{d}{dr}\!\left[r^{2}\!\left(\phi_{m}\frac{d\phi_{n}}{dr}
- \phi_{n}\frac{d\phi_{m}}{dr}\right)\right]
= \phi_{m}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{n}}{dr}\right)
- \phi_{n}\,\frac{d}{dr}\!\left(r^{2}\frac{d\phi_{m}}{dr}\right)
$$

Integrating from $0$ to $R$, the boundary term vanishes at the origin because of the factor $r^{2}$:

$$
\left[r^{2}\!\left(\phi_{m}\frac{d\phi_{n}}{dr}
- \phi_{n}\frac{d\phi_{m}}{dr}\right)\right]_{r=0} = 0
$$

and vanishes at the surface because both eigenfunctions obey the same Robin relation:

$$
\left[\phi_{m}\frac{d\phi_{n}}{dr} - \phi_{n}\frac{d\phi_{m}}{dr}\right]_{r=R}
= \phi_{m}\left(-\frac{h}{D}\phi_{n}\right)
- \phi_{n}\left(-\frac{h}{D}\phi_{m}\right) = 0
$$

Hence

$$
\left(\lambda_{m} - \lambda_{n}\right)\int_{0}^{R} r^{2}\phi_{n}\phi_{m}\,dr = 0
$$

and since $\lambda_{m} \neq \lambda_{n}$:

$$
\boxed{\;
\int_{0}^{R} \phi_{n}\,\phi_{m}\, r^{2}\,dr = 0,
\qquad n \neq m
\;}
$$

---

## 12. Normalization

With $\phi_{n}^{2} r^{2} = \sin^{2}\!\left(\sqrt{\lambda_{n}}\,r\right)$:

$$
N_{n} = \int_{0}^{R} \phi_{n}^{2}\, r^{2}\,dr
= \int_{0}^{R} \sin^{2}\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

Applying the half-angle identity:

$$
N_{n}
= \int_{0}^{R} \frac{1 - \cos\!\left(2\sqrt{\lambda_{n}}\,r\right)}{2}\,dr
= \frac{R}{2} - \frac{\sin\!\left(2\sqrt{\lambda_{n}}\,R\right)}{4\sqrt{\lambda_{n}}}
$$

Using $\sin(2\theta) = 2\sin\theta\cos\theta$:

$$
N_{n} = \frac{R}{2}
- \frac{\sin\!\left(\sqrt{\lambda_{n}}\,R\right)
\cos\!\left(\sqrt{\lambda_{n}}\,R\right)}{2\sqrt{\lambda_{n}}}
$$

Eliminating the cosine with the identity from Section 9:

$$
N_{n} = \frac{R}{2}
- \frac{\left(D - hR\right)\sin^{2}\!\left(\sqrt{\lambda_{n}}\,R\right)}
{2D\,\lambda_{n} R}
$$

$$
\boxed{\;
N_{n} = \frac{R}{2}\left[1
- \frac{\left(D - hR\right)\sin^{2}\!\left(\sqrt{\lambda_{n}}\,R\right)}
{D\,\lambda_{n} R^{2}}\right]
\;}
$$

---

## 13. Projection of the initial condition

Sturm–Liouville theory guarantees that the eigenfunctions form a complete basis, so

$$
v(r,0) = \sum_{n \ge 1} a_{n}\,\phi_{n}(r)
$$

Multiplying by $\phi_{n} r^{2}$ and integrating, all cross terms vanish by orthogonality, leaving

$$
\int_{0}^{R} v(r,0)\,\phi_{n}(r)\, r^{2}\,dr = a_{n} N_{n}
\qquad\Longrightarrow\qquad
a_{n} = \frac{1}{N_{n}}\int_{0}^{R} v(r,0)\,\phi_{n}(r)\, r^{2}\,dr
$$

Using $v(r,0) = -w(r)$ and $\phi_{n}(r)\,r^{2} = r\sin\!\left(\sqrt{\lambda_{n}}\,r\right)$:

$$
a_{n} = -\frac{1}{N_{n}}
\int_{0}^{R} w(r)\, r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

Writing $w(r)$ in the form of a constant minus a quadratic:

$$
w(r) = \left(\frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}\right)
- \frac{S_{0}}{6D}\,r^{2}
$$

the required integral separates as

$$
\int_{0}^{R} w(r)\, r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= \left(\frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}\right) I_{1}
- \frac{S_{0}}{6D}\, I_{3}
$$

where

$$
I_{1} = \int_{0}^{R} r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr,
\qquad
I_{3} = \int_{0}^{R} r^{3}\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

---

## 14. Evaluation of the integrals

### 14.1 The integral $I_1$

Integrating by parts with $u = r$ and $dv = \sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr$:

$$
I_{1}
= \left[-\frac{r\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}\right]_{0}^{R}
+ \frac{1}{\sqrt{\lambda_{n}}}\int_{0}^{R}\cos\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

Evaluating both pieces:

$$
I_{1}
= -\frac{R\cos\!\left(\sqrt{\lambda_{n}}\,R\right)}{\sqrt{\lambda_{n}}}
+ \frac{1}{\sqrt{\lambda_{n}}}
\left[\frac{\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}\right]_{0}^{R}
= \frac{\sin\!\left(\sqrt{\lambda_{n}}\,R\right)
- \sqrt{\lambda_{n}}\,R\cos\!\left(\sqrt{\lambda_{n}}\,R\right)}{\lambda_{n}}
$$

From the cosine-elimination identity:

$$
\sqrt{\lambda_{n}}\,R\cos\!\left(\sqrt{\lambda_{n}}\,R\right)
= \frac{\left(D - hR\right)\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D}
$$

so that

$$
I_{1}
= \frac{\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{\lambda_{n}}
\left[1 - \frac{D - hR}{D}\right]
$$

$$
\boxed{\;
I_{1} = \frac{hR\,\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D\,\lambda_{n}}
\;}
$$

### 14.2 The integral $I_3$

Integrating by parts with $u = r^{3}$:

$$
\int r^{3}\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= -\frac{r^{3}\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}
+ \frac{3}{\sqrt{\lambda_{n}}}\int r^{2}\cos\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

The remaining integral is reduced by a further integration by parts with $u = r^{2}$:

$$
\int r^{2}\cos\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= \frac{r^{2}\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}
- \frac{2}{\sqrt{\lambda_{n}}}\int r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
$$

and the innermost integral is

$$
\int r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= -\frac{r\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}
+ \frac{\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}}
$$

Back-substituting:

$$
\int r^{2}\cos\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= \frac{r^{2}\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}
+ \frac{2r\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}}
- \frac{2\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}^{3/2}}
$$

and therefore

$$
\int r^{3}\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= -\frac{r^{3}\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\sqrt{\lambda_{n}}}
+ \frac{3r^{2}\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}}
+ \frac{6r\cos\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}^{3/2}}
- \frac{6\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{\lambda_{n}^{2}}
$$

Evaluating between $0$ and $R$, where every term vanishes at the lower limit, and grouping the sine and cosine contributions:

$$
I_{3}
= \left(\frac{3R^{2}}{\lambda_{n}} - \frac{6}{\lambda_{n}^{2}}\right)
\sin\!\left(\sqrt{\lambda_{n}}\,R\right)
- \left(\frac{R^{3}}{\sqrt{\lambda_{n}}} - \frac{6R}{\lambda_{n}^{3/2}}\right)
\cos\!\left(\sqrt{\lambda_{n}}\,R\right)
$$

Eliminating the cosine:

$$
- \left(\frac{R^{3}}{\sqrt{\lambda_{n}}} - \frac{6R}{\lambda_{n}^{3/2}}\right)
\cos\!\left(\sqrt{\lambda_{n}}\,R\right)
= -\frac{\left(D - hR\right)\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D}
\left(\frac{R^{2}}{\lambda_{n}} - \frac{6}{\lambda_{n}^{2}}\right)
$$

so that

$$
I_{3} = \sin\!\left(\sqrt{\lambda_{n}}\,R\right)
\left[\frac{R^{2}}{\lambda_{n}}\left(3 - \frac{D - hR}{D}\right)
- \frac{6}{\lambda_{n}^{2}}\left(1 - \frac{D - hR}{D}\right)\right]
$$

Simplifying the two bracketed coefficients:

$$
3 - \frac{D - hR}{D} = \frac{2D + hR}{D},
\qquad
1 - \frac{D - hR}{D} = \frac{hR}{D}
$$

gives

$$
\boxed{\;
I_{3} = \sin\!\left(\sqrt{\lambda_{n}}\,R\right)
\left[\frac{\left(2D + hR\right)R^{2}}{D\,\lambda_{n}}
- \frac{6hR}{D\,\lambda_{n}^{2}}\right]
\;}
$$

---

## 15. Expansion coefficients

Substituting $I_1$ and $I_3$ and grouping by powers of $\lambda_{n}$:

$$
\int_{0}^{R} w(r)\, r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= \frac{\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{\lambda_{n}}
\left[\left(\frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}\right)\frac{hR}{D}
- \frac{S_{0}}{6D}\cdot\frac{\left(2D + hR\right)R^{2}}{D}\right]
+ \frac{S_{0}}{6D}\cdot
\frac{6hR\,\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D\,\lambda_{n}^{2}}
$$

The bracketed term vanishes identically. Expanding its first part:

$$
\left(\frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}\right)\frac{hR}{D}
= \frac{S_{0}R^{2}}{3D} + \frac{S_{0}hR^{3}}{6D^{2}}
$$

and its second part:

$$
\frac{S_{0}}{6D}\cdot\frac{\left(2D + hR\right)R^{2}}{D}
= \frac{2S_{0}R^{2}D}{6D^{2}} + \frac{S_{0}hR^{3}}{6D^{2}}
= \frac{S_{0}R^{2}}{3D} + \frac{S_{0}hR^{3}}{6D^{2}}
$$

so that

$$
\left(\frac{S_{0}R}{3h} + \frac{S_{0}R^{2}}{6D}\right)\frac{hR}{D}
- \frac{S_{0}}{6D}\cdot\frac{\left(2D + hR\right)R^{2}}{D} = 0
$$

This cancellation occurs only when the steady state and the eigencondition are mutually consistent. What remains is

$$
\int_{0}^{R} w(r)\, r\sin\!\left(\sqrt{\lambda_{n}}\,r\right) dr
= \frac{S_{0}}{6D}\cdot
\frac{6hR\,\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D\,\lambda_{n}^{2}}
= \frac{hR\,S_{0}\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}{D^{2}\,\lambda_{n}^{2}}
$$

and therefore

$$
\boxed{\;
a_{n} = -\frac{hR\,S_{0}\,\sin\!\left(\sqrt{\lambda_{n}}\,R\right)}
{D^{2}\,\lambda_{n}^{2}\,N_{n}}
\;}
$$

---

## 16. Solution with the source active

Reassembling $c = w + v$:

$$
\boxed{\;
c(r,t) = \frac{S_{0}R}{3h} + \frac{S_{0}}{6D}\left(R^{2} - r^{2}\right)
+ \sum_{n \ge 1} a_{n}\,
\frac{\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{r}\;
e^{-\lambda_{n} D t},
\qquad 0 \le t \le t_{c}
\;}
$$

At $t = 0$ every exponential equals one and the series reproduces $-w(r)$ by construction, so the initial condition is met exactly. As $t \to \infty$ every mode decays and the profile relaxes to the steady parabola.

The instantaneous release rate through the surface is

$$
\dot{m}(t) = 4\pi R^{2}\, h\, c(R,t)
\qquad \left(\mathrm{mol\ s^{-1}}\right)
$$

and the slowest mode sets the characteristic equilibration time

$$
\tau_{1} = \frac{1}{\lambda_{1} D}
$$

---

## 17. Solution after source switch-off

For $t > t_{c}$ the governing equation is homogeneous with unchanged boundary conditions:

$$
\frac{\partial c}{\partial t}
= \frac{D}{r^{2}}\,\frac{\partial}{\partial r}\!\left(r^{2}\frac{\partial c}{\partial r}\right),
\qquad t > t_{c}
$$

so the same eigenbasis applies:

$$
c(r,t) = \sum_{n \ge 1} b_{n}\,
\frac{\sin\!\left(\sqrt{\lambda_{n}}\,r\right)}{r}\;
e^{-\lambda_{n} D\,(t - t_{c})}
$$

The new initial state is the phase-one solution evaluated at $t_{c}$. Since the projection of $w$ onto the eigenbasis is $-a_{n}$ by definition:

$$
w(r) = -\sum_{n \ge 1} a_{n}\,\phi_{n}(r)
$$

and therefore

$$
c(r,t_{c})
= w(r) + \sum_{n \ge 1} a_{n}\,\phi_{n}(r)\, e^{-\lambda_{n} D t_{c}}
= \sum_{n \ge 1}\left(-a_{n} + a_{n} e^{-\lambda_{n} D t_{c}}\right)\phi_{n}(r)
$$

giving

$$
\boxed{\;
b_{n} = a_{n}\left(e^{-\lambda_{n} D t_{c}} - 1\right)
\;}
$$

No additional integrals are required: orthogonality of the eigenbasis transfers the phase-one state directly. Since $\lambda_{n} > 0$, the inventory decays exponentially after shutdown, with the late-time behaviour governed by $\lambda_{1}$.

---

## 18. Extension to the five-layer particle

The full TRISO geometry comprises five concentric layers (kernel, buffer, IPyC, SiC, OPyC) with radii

$$
0 = r_{0} < r_{1} < r_{2} < r_{3} < r_{4} < r_{5} = R
$$

and piecewise constant diffusivities $D_{i}$. In each layer:

$$
\frac{\partial c}{\partial t}
= \frac{D_{i}}{r^{2}}\,\frac{\partial}{\partial r}\!\left(r^{2}\frac{\partial c}{\partial r}\right)
+ S_{i},
\qquad i = 1,\dots,5
$$

with the source confined to the kernel:

$$
S_{1} = S_{0},
\qquad
S_{2} = S_{3} = S_{4} = S_{5} = 0
$$

### Interface conditions

Both the concentration and the radial flux are continuous at each interface:

$$
c\!\left(r_{i}^{-},t\right) = c\!\left(r_{i}^{+},t\right),
\qquad
D_{i}\,\frac{\partial c}{\partial r}\bigg|_{r = r_{i}^{-}}
= D_{i+1}\,\frac{\partial c}{\partial r}\bigg|_{r = r_{i}^{+}},
\qquad i = 1,\dots,4
$$

### Composite spatial problem

Writing $D(r)$ for the piecewise diffusivity, the spatial problem becomes

$$
\frac{d}{dr}\!\left(D(r)\, r^{2}\frac{d\phi}{dr}\right) + \lambda\, r^{2}\phi = 0
$$

with time factor $e^{-\lambda t}$, the diffusivity now residing inside the derivative rather than in the time equation. Flux continuity is exactly continuity of $D(r)\,r^{2}\,d\phi/dr$, which is the condition under which self-adjointness, orthogonality

$$
\int_{0}^{R} \phi_{n}\,\phi_{m}\, r^{2}\,dr = 0,
\qquad n \neq m
$$

and realness of the spectrum survive the material discontinuities.

### Layer-wise solution

With $u = r\phi$ in each layer:

$$
u^{(i)}(r) = \alpha_{i}\sin\!\left(\sqrt{\frac{\lambda}{D_{i}}}\;r\right)
+ \beta_{i}\cos\!\left(\sqrt{\frac{\lambda}{D_{i}}}\;r\right),
\qquad
\beta_{1} = 0
$$

by regularity at the origin. The remaining constants satisfy eight interface equations together with the outer Robin condition, forming a homogeneous linear system conveniently assembled as a product of $2 \times 2$ transfer matrices across the layers. Nontrivial solutions exist only where the system determinant vanishes:

$$
\det M(\lambda) = 0
$$

which is the transcendental eigencondition generalizing the single-layer case. The steady-state split, the projection onto the eigenbasis, and the switch-off treatment then carry over unchanged.
