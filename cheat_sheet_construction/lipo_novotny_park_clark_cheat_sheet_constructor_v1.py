# Retry with mathtext-safe formulas (no LaTeX environments). Create the one-page PDF.
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

width, height = 8.5, 11
ppath = "Park_Clarke_SpaceVector_OnePage.pdf"

pdf = PdfPages(ppath)
fig = plt.figure(figsize=(width, height))

# Layout
margin_left = 0.06
y = 0.975
sp_title = 0.038
sp_line = 0.028
fs_title = 12
fs_text = 10
fs_small = 9

def section(title):
    global y
    y -= sp_title
    fig.text(margin_left, y, title, fontsize=fs_title, fontweight="bold")
    y -= 0.016

def line(tex, size=fs_text, x=None):
    global y
    if x is None:
        x = margin_left + 0.01
    fig.text(x, y, tex, fontsize=size)
    y -= sp_line

# Title
fig.text(0.5, 0.99, "Park–Clarke & Space Vector Cheat Sheet (Lipo–Novotny)",
         ha="center", va="top", fontsize=16, fontweight="bold")

# 1) Conventions
section("1) Conventions (UW Lipo–Novotny)")
line("• q-axis ∥ α at θ=0; d is −90° from q (clockwise).", fs_small)
line("• θ: electrical angle α→q, positive CCW;  ω_e = dθ/dt.", fs_small)
line("• β plotted as −Imag{·}; q = Real{·},  d = −Imag{·}.", fs_small)

# 2) Frames
section("2) Frames")
line("abc: {f_a,f_b,f_c} phase quantities (stator).", fs_small)
line("αβ0: {f_α,f_β,f_0} stationary orthogonal frame (Clarke).", fs_small)
line("qd0: {f_q,f_d,f_0} rotor-synchronous frame (Park).", fs_small)

# 3) Clarke
section("3) Clarke (abc → αβ0) — balanced 3ϕ")
line(r"f_α = (2/3)[ f_a − 0.5 f_b − 0.5 f_c ]", fs_small)
line(r"f_β = (2/3)[ (√3/2)( f_b − f_c ) ]", fs_small)
line(r"f_0 = (1/3)( f_a + f_b + f_c )", fs_small)

# 4) Park
section("4) Park (αβ0 → qd0) and reverse")
line(r"Let cθ=cosθ, sθ=sinθ.", fs_small)
line(r"f_q =  cθ f_α − sθ f_β;     f_d =  sθ f_α + cθ f_β", fs_small)
line(r"f_α =  cθ f_q + sθ f_d;     f_β = −sθ f_q + cθ f_d", fs_small)

# 5) Complex forms
section("5) Complex space vectors")
line(r"a = e^{j2π/3}", fs_small)
line(r"f_{αβ} = (2/3)( f_a + a f_b + a^2 f_c )", fs_small)
line(r"Park:  f_{qd} = f_{αβ} e^{−jθ}  ⇒  q=Real{f_{qd}}, d=−Imag{f_{qd}}", fs_small)
line(r"Back-EMF:  e_{αβ} = dλ_{αβ}/dt = j ω_e λ_{αβ}  (EMF leads flux by 90°)", fs_small)

# 6) Symbols
section("6) Symbols & meanings")
left = [
    "α, β : Clarke axes (stationary)",
    "q, d : Quadrature/Direct axes (rotating)",
    "θ : electrical angle α→q [rad]",
    "ω_e : electrical speed [rad/s]",
    "λ, Λ : flux linkage",
    "e : back-EMF",
    "i : current",
]
right = [
    "v : voltage",
    "p : pole pairs",
    "f_0 : zero-sequence component",
    "Real{·}, Imag{·} : real/imag parts",
    "j : √(−1) (90° rotation)",
    "a = e^{j2π/3} : 120° operator",
    "β plotting: often use −Imag{·}",
]
block_y = y
for r in left:
    line(r, fs_small)
# Right column
y2 = block_y
for r in right:
    fig.text(0.52, y2, r, fontsize=fs_small)
    y2 -= sp_line
y = min(y, y2)

# 7) Pitfalls
section("7) Common pitfalls")
line("• Using angle-to-d (standard) instead of angle-to-q (Lipo–Novotny).", fs_small)
line("• Forgetting minus on β or on d = −Imag{·}.", fs_small)
line("• Using mechanical speed ω_m when formulas need ω_e = p·ω_m.", fs_small)

# Footer
fig.text(0.5, 0.02, "One-page quick ref — UW Lipo–Novotny convention.  e = dλ/dt ⇒ EMF ⟂ flux; q∥α at θ=0.",
         ha="center", fontsize=8)

plt.axis('off')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()
