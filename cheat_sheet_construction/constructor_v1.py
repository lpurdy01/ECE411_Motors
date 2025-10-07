# Second retry with simplified mathtext expressions (avoid \ge, \Rightarrow, etc.).

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

width, height = 8.5, 11
ppath = "ECE411_DC_Machine_CheatSheet_v2.pdf"

pdf = PdfPages(ppath)
fig = plt.figure(figsize=(width, height))

left = 0.06
y = 0.96
dy_title = 0.035
dy = 0.028
small = 9
normal = 10
tiny = 8.3

def section(title):
    global y
    y -= dy_title
    fig.text(left, y, title, fontsize=12, fontweight="bold")
    y -= 0.007

def line(text, size=normal, x=None):
    global y
    if x is None:
        x = left + 0.01
    fig.text(x, y, text, fontsize=size)
    y -= dy

fig.text(0.5, 0.985, "ECE 411 — DC Machine & Control Systems Cheat Sheet", ha="center", va="top", fontsize=15, fontweight="bold")

section("1) Electrical Model")
line(r"$V_a = R_a\,I_a + L_a\,\dfrac{d I_a}{dt} + e$")
line(r"$e = K_e\,\omega_m$    (Ke = Kt in SI)")
line(r"$V_a$ [V], $R_a$ [\Omega], $L_a$ [H], $I_a$ [A], $e$ [V], $K_e$ [\mathrm{V\,s/rad}], $\omega_m$ [rad/s]", size=small)

section("2) Mechanical Model")
line(r"$J_m\,\dfrac{d\omega_m}{dt} = T_e - T_L$")
line(r"$T_e = K_t\,I_a$")
line(r"$J_m$ [kg\,m^2], $T_e$ [N\,m], $T_L$ [N\,m], $K_t$ [N\,m/A], $\omega_m$ [rad/s]", size=small)

section("3) Time Constants")
line(r"$\tau_a = \dfrac{L_a}{R_a}, \qquad \tau_m = \dfrac{J_m R_a}{K_t^2}$")
line(r"$\tau_a, \tau_m$ [s]", size=small)

section("4) Characteristic Equation (open-loop speed dynamics)")
line(r"$p^2 + \dfrac{1}{\tau_a}\, p + \dfrac{1}{\tau_a \tau_m} = 0$")

section("5) Natural Frequency & Damping Ratio")
line(r"$\omega_n = \dfrac{1}{\sqrt{\tau_a \tau_m}}, \qquad \zeta = \dfrac{1}{2}\sqrt{\dfrac{\tau_m}{\tau_a}}$")
line(r"$\omega_n$ [rad/s], $\zeta$ [1]", size=small)

section("6) Eigenvalues (roots)")
line(r"$p_{1,2} = -\dfrac{1}{2\tau_a} \pm \sqrt{\left(\dfrac{1}{2\tau_a}\right)^{2} - \dfrac{1}{\tau_a \tau_m}}$")

section("7) Minimum Inertia for Real-Valued Roots")
line("Require discriminant >= 0 for real roots:  tau_m >= 4 tau_a")
line(r"$J_{m,\min} = \dfrac{4 L_a K_t^2}{R_a^2}$")

section("8) Performance Metrics (2% criterion)")
line(r"Settling time: $\mathrm{ST} = \dfrac{4}{\omega_n \zeta}$")
line(r"Percent overshoot: $\mathrm{OS} = 100\,\exp\!\left(-\dfrac{\zeta \pi}{\sqrt{\,1-\zeta^2\,}}\right)\%$")
line(r"General solution: $x(t)=X\,e^{p t}$", size=normal)

section("9) Impedance Analogy (Electrical ↔ Rotational Mechanical)")
line(r"$V \leftrightarrow T$,  $I \leftrightarrow \omega$,  $R \leftrightarrow B$,  $L \leftrightarrow J$,  $C \leftrightarrow 1/B$", size=normal)

section("10) Variable & Unit Quick Reference")
y -= 0.01
col1 = [
    r"$V_a$ : Armature voltage [V]",
    r"$I_a$ : Armature current [A]",
    r"$R_a$ : Armature resistance [\Omega]",
    r"$L_a$ : Armature inductance [H]",
    r"$e$ : Back-EMF [V]",
    r"$K_t$ : Torque const. [N\,m/A]",
    r"$K_e$ : Voltage const. [V\,s/rad]",
]
col2 = [
    r"$T_e$ : Electromagnetic torque [N\,m]",
    r"$J_m$ : Inertia [kg\,m^2]",
    r"$\omega_m$ : Angular speed [rad/s]",
    r"$\tau_a,\tau_m$ : Time constants [s]",
    r"$\zeta$ : Damping ratio [1]",
    r"$\omega_n$ : Natural freq. [rad/s]",
    r"$p$ : Laplace variable [s^{-1}]",
]

start_y = y
for t in col1:
    line(t, size=tiny)

y_right = start_y
for t in col2:
    fig.text(0.52, y_right, t, fontsize=tiny)
    y_right -= dy

fig.text(0.5, 0.02, "Assumes viscous friction B ~= 0 in time-constant formulas; Ke = Kt in SI.", ha="center", fontsize=8)

fig.patch.set_alpha(0.0)
plt.axis('off')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()
