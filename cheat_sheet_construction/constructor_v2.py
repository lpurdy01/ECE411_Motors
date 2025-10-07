# Create an improved layout PDF (v3) and also write a LaTeX .tex source file.
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---------- Improved matplotlib PDF ----------
width, height = 8.5, 11
ppath = "ECE411_DC_Machine_CheatSheet_v3.pdf"

pdf = PdfPages(ppath)
fig = plt.figure(figsize=(width, height))

# Layout params
margin_left = 0.06
margin_right = 0.95
y = 0.965
sp_title = 0.04   # space after section title
sp_line = 0.03    # space between lines
fs_title = 12
fs_text = 10
fs_small = 9

def section(title):
    global y
    y -= sp_title
    fig.text(margin_left, y, title, fontsize=fs_title, fontweight="bold")
    y -= 0.02  # extra gap before first line

def line_math(tex, size=fs_text, x=None):
    global y
    if x is None:
        x = margin_left + 0.01
    fig.text(x, y, tex, fontsize=size)
    y -= sp_line

def line_text(txt, size=fs_text, x=None):
    global y
    if x is None:
        x = margin_left + 0.01
    fig.text(x, y, txt, fontsize=size)
    y -= sp_line

# Title
fig.text(0.5, 0.985, "ECE 411 — DC Machine & Control Systems Cheat Sheet",
         ha="center", va="top", fontsize=16, fontweight="bold")

# Sections
section("1) Electrical Model")
line_math(r"$V_a = R_a\,I_a + L_a\,\dfrac{d I_a}{dt} + e$")
line_math(r"$e = K_e\,\omega_m$   (Ke = Kt in SI)")
line_math(r"$V_a$ [V], $R_a$ [\Omega], $L_a$ [H], $I_a$ [A], $e$ [V], $K_e$ [\mathrm{V\,s/rad}], $\omega_m$ [rad/s]", fs_small)

section("2) Mechanical Model")
line_math(r"$J_m\,\dfrac{d\omega_m}{dt} = T_e - T_L$")
line_math(r"$T_e = K_t\,I_a$")
line_math(r"$J_m$ [kg\,m$^2$], $T_e$ [N\,m], $T_L$ [N\,m], $K_t$ [N\,m/A], $\omega_m$ [rad/s]", fs_small)

section("3) Time Constants")
line_math(r"$\tau_a = \dfrac{L_a}{R_a}, \qquad \tau_m = \dfrac{J_m R_a}{K_t^2}$")
line_math(r"$\tau_a, \tau_m$ [s]", fs_small)

section("4) Characteristic Equation (open-loop speed dynamics)")
line_math(r"$p^2 + \dfrac{1}{\tau_a}\, p + \dfrac{1}{\tau_a \tau_m} = 0$")

section("5) Natural Frequency & Damping Ratio")
line_math(r"$\omega_n = \dfrac{1}{\sqrt{\tau_a \tau_m}}, \qquad \zeta = \dfrac{1}{2}\sqrt{\dfrac{\tau_m}{\tau_a}}$")
line_math(r"$\omega_n$ [rad/s], $\zeta$ [1]", fs_small)

section("6) Eigenvalues (roots)")
line_math(r"$p_{1,2} = -\dfrac{1}{2\tau_a} \pm \sqrt{\left(\dfrac{1}{2\tau_a}\right)^{2} - \dfrac{1}{\tau_a \tau_m}}$")

section("7) Minimum Inertia for Real-Valued Roots")
line_text("Require discriminant >= 0 for real roots:  tau_m >= 4 tau_a", fs_text)
line_math(r"$J_{m,\min} = \dfrac{4 L_a K_t^2}{R_a^2}$")

section("8) Performance Metrics (2% criterion)")
line_math(r"Settling time: $\mathrm{ST} = \dfrac{4}{\omega_n \zeta}$")
line_math(r"Percent overshoot: $\mathrm{OS} = 100\,\exp\!\left(-\dfrac{\zeta \pi}{\sqrt{\,1-\zeta^2\,}}\right)\%$")
line_math(r"General solution: $x(t)=X\,e^{p t}$")

section("9) Impedance Analogy (Electrical ↔ Rotational Mechanical)")
line_math(r"$V \leftrightarrow T$,  $I \leftrightarrow \omega$,  $R \leftrightarrow B$,  $L \leftrightarrow J$,  $C \leftrightarrow 1/B$")

# Units table (render as plain text rows to avoid mathtext parser issues)
section("10) Variable & Unit Quick Reference")
rows_left = [
    "Va : Armature voltage [V]",
    "Ia : Armature current [A]",
    "Ra : Armature resistance [Ω]",
    "La : Armature inductance [H]",
    "e  : Back-EMF [V]",
    "Kt : Torque constant [N·m/A]",
    "Ke : Voltage constant [V·s/rad]",
]
rows_right = [
    "Te : Electromagnetic torque [N·m]",
    "Jm : Inertia [kg·m²]",
    "ωm : Angular speed [rad/s]",
    "τa, τm : Time constants [s]",
    "ζ : Damping ratio [–]",
    "ωn : Natural frequency [rad/s]",
    "p : Laplace variable [s⁻¹]",
]

# Remember where we start this block
block_y = y
for r in rows_left:
    line_text(r, fs_small)

# Right column
y2 = block_y
for r in rows_right:
    fig.text(0.52, y2, r, fontsize=fs_small)
    y2 -= sp_line

# Footer
fig.text(0.5, 0.02, "Assumes viscous friction B ≈ 0 in time-constant formulas; Ke = Kt in SI.", ha="center", fontsize=8)

plt.axis('off')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()

# ---------- LaTeX file ----------
latex_path = "ece411_dc_machine_cheatsheet.tex"
latex_src = r"""
\documentclass[11pt]{article}
\usepackage[letterpaper,margin=0.7in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{multicol}
\usepackage{siunitx}
\usepackage{titlesec}
\titleformat{\section}{\large\bfseries}{}{0pt}{}
\setlength{\parskip}{0.3em}
\setlength{\parindent}{0pt}
\begin{document}
\begin{center}
{\LARGE \textbf{ECE 411 — DC Machine \& Control Systems Cheat Sheet}}\\[2mm]
\small (Ke = Kt in SI units; viscous friction $B\approx0$ unless noted)
\end{center}

\section{Electrical Model}
\begin{align*}
V_a &= R_a I_a + L_a \frac{d I_a}{dt} + e, \qquad
e = K_e\,\omega_m
\end{align*}

\section{Mechanical Model}
\begin{align*}
J_m \frac{d\omega_m}{dt} &= T_e - T_L, \qquad
T_e = K_t I_a
\end{align*}

\section{Time Constants}
\begin{align*}
\tau_a = \frac{L_a}{R_a}, \qquad
\tau_m = \frac{J_m R_a}{K_t^2}
\end{align*}

\section{Characteristic Equation (open-loop speed)}
\begin{align*}
p^2 + \frac{1}{\tau_a}p + \frac{1}{\tau_a \tau_m} = 0
\end{align*}

\section{Natural Frequency \& Damping}
\begin{align*}
\omega_n = \frac{1}{\sqrt{\tau_a \tau_m}}, \qquad
\zeta = \frac{1}{2}\sqrt{\frac{\tau_m}{\tau_a}}
\end{align*}

\section{Eigenvalues}
\begin{align*}
p_{1,2} = -\frac{1}{2\tau_a} \pm \sqrt{\left(\frac{1}{2\tau_a}\right)^2 - \frac{1}{\tau_a\tau_m}}
\end{align*}

\section{Minimum Inertia for Real Roots}
Real roots require discriminant $\ge 0 \Rightarrow \tau_m \ge 4 \tau_a$ :
\begin{align*}
J_{m,\min} = \frac{4 L_a K_t^2}{R_a^2}
\end{align*}

\section{Performance Metrics (2\% criterion)}
\begin{align*}
\mathrm{ST} &= \frac{4}{\omega_n \zeta}, \qquad
\mathrm{OS} = 100\,e^{-\zeta \pi/\sqrt{1-\zeta^2}}\% \\
x(t) &= X e^{pt}
\end{align*}

\section{Impedance Analogy (Electrical $\leftrightarrow$ Rotational Mechanical)}
\begin{center}
$V \!\leftrightarrow\! T,\quad I \!\leftrightarrow\! \omega,\quad R \!\leftrightarrow\! B,\quad L \!\leftrightarrow\! J,\quad C \!\leftrightarrow\! 1/B$
\end{center}

\section{Variable \& Unit Quick Reference}
\begin{multicols}{2}
\begin{tabular}{ll}
$V_a$ & Armature voltage [V] \\
$I_a$ & Armature current [A] \\
$R_a$ & Armature resistance [\si{\ohm}]\\
$L_a$ & Armature inductance [H] \\
$e$   & Back-EMF [V] \\
$K_t$ & Torque constant [N\,m/A] \\
$K_e$ & Voltage constant [V\,s/rad] \\
\end{tabular}

\columnbreak

\begin{tabular}{ll}
$T_e$ & Electromagnetic torque [N\,m] \\
$J_m$ & Inertia [kg\,m$^2$] \\
$\omega_m$ & Angular speed [rad/s] \\
$\tau_a,\tau_m$ & Time constants [s] \\
$\zeta$ & Damping ratio [–] \\
$\omega_n$ & Natural frequency [rad/s] \\
$p$ & Laplace variable [s$^{-1}$] \\
\end{tabular}
\end{multicols}

\end{document}
"""
with open(latex_path, "w") as f:
    f.write(latex_src)

(ppath, latex_path)
