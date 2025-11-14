import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Page setup
width, height = 8.5, 11
ppath = "ECE411_HW3_HW4_CheatSheet.pdf"

pdf = PdfPages(ppath)
fig = plt.figure(figsize=(width, height))

# Layout parameters
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
fig.text(0.5, 0.99, "ECE 411 – HW3 & HW4 Exam Cheat Sheet (Lipo–Novotny Conventions)",
         ha="center", va="top", fontsize=16, fontweight="bold")

# 1) Reference Frames & Conventions
section("1) Frames & Conventions")
line("• q-axis ∥ α at θ=0; d is −90° from q (clockwise).", fs_small)
line("• θ_e: electrical angle (α→q, CCW); ω_e = dθ_e/dt.", fs_small)
line("• q = Real{·}, d = −Imag{·}; β = −Imag{αβ} for plotting.", fs_small)

# 2) Clarke & Park Transforms
section("2) Clarke & Park Transforms")
line(r"Clarke: α = (2/3)[ Va − 0.5Vb − 0.5Vc ]", fs_small)
line(r"β = (2/3)(√3/2)(Vc − Vb)", fs_small)
line(r"Park: q =  α cosθ − β sinθ", fs_small)
line(r"       d =  α sinθ + β cosθ", fs_small)
line(r"Inverse: α =  q cosθ + d sinθ;   β = −q sinθ + d cosθ", fs_small)

# 3) Phasor Definitions
section("3) Phasor Definitions")
line(r"a = e^{j2π/3};   x_{αβ} = (2/3)(xa + a xb + a² xc)", fs_small)
line(r"Park: x_{qd} = x_{αβ} e^{−jθ_e};   q = Re{x_{qd}}, d = −Im{x_{qd}}", fs_small)
line(r"E_0 = ω_e Λ_0;   Λ_0 = E_0 / ω_e", fs_small)

# 4) Base Quantities (per-unit)
section("4) Base Quantities (per-unit)")
line(r"V_B = V_{dc}/√3,   I_B = I_{max},   Z_B = V_B/I_B", fs_small)
line(r"P_B = 1.5 V_B I_B,   ω_B = 2πf_B,   T_B = P_B/ω_{mB}", fs_small)

# 5) Voltage Equations
section("5) Voltage Equations")
line(r"Round-rotor:  V_q = R_s I_q + ω_e L_s I_d + E_q", fs_small)
line(r"               V_d = R_s I_d − ω_e L_s I_q + E_d", fs_small)
line(r"Salient-pole: V_q = R_s I_q + ω_e L_d I_d + E_q", fs_small)
line(r"               V_d = R_s I_d − ω_e L_q I_q + E_d", fs_small)

# 6) Torque & Power
section("6) Torque & Power")
line(r"P = 1.5(V_q I_q + V_d I_d)", fs_small)
line(r"T_e = P / ω_e", fs_small)
line(r"T_e = (E_q I_q)/ω_e + ((L_d−L_q) I_d I_q)/ω_e", fs_small)
line(r"Salient-pole:", fs_small)
line(r"   T_field = −(E_0 V_0)/(ω_e L_d) sinδ", fs_small)
line(r"   T_rel = −(V_0²(L_d−L_q))/(2ω_e L_d L_q) sin(2δ)", fs_small)
line(r"   T_e = T_field + T_rel", fs_small)

# 7) Field Weakening
section("7) Field Weakening")
line(r"Voltage limit (R_s=0):  V_0² = (ω_e Λ_0)² − (ω_e L_s I_0)²", fs_small)
line(r"Λ_0 = √[(V_0/ω_e)² + (L_s I_0)²]", fs_small)
line(r"T_max = V_0 I_0 / ω_e;   P_max = V_0 I_0", fs_small)
line(r"Flux ratio = Λ_new / Λ_rated", fs_small)

# 8) Rs ≠ 0 Effects
section("8) Rs ≠ 0 Effects")
line(r"Z_s = √(R_s² + (ω_e L_s)²)", fs_small)
line(r"T_max = (V_0E_0)/(ω_e Z_s) − (E_0²R_s)/(ω_e Z_s²)", fs_small)
line(r"δ_max = sin⁻¹(R_s/Z_s) − π/2", fs_small)

# 9) Saliency & Ratios
section("9) Saliency & Ratios")
line(r"Saliency ratio ξ = L_d / L_q", fs_small)
line(r"T_rel,max = T_field,max → L_q = (V_0 L_d)/(2E_0 + V_0 R_s)", fs_small)

# 10) Quick Reference
section("10) Quick Reference (exam tips)")
line("• Open-circuit: vq ≈ ω_e Λ_0, vd ≈ 0 (confirm q–α alignment).", fs_small)
line("• Generator: P_mech = −P_elec; motoring → +T_e.", fs_small)
line("• δ increases with load torque; small δ ≈ linear region.", fs_small)
line("• Keep per-unit bases consistent (V_B, I_B, L in pu).", fs_small)
line("• q = torque-producing axis; d = flux (field) axis.", fs_small)

# Footer
fig.text(0.5, 0.02, "One-page review — HW3 & HW4 core equations under Lipo–Novotny convention.",
         ha="center", fontsize=8)

plt.axis('off')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()
