import fitz
doc = fitz.open('drozdov2015.pdf')
page = doc[8]
drawings = page.get_drawings()

# axis calibration (panel c), derived from vector tick-mark rects
X0, X0_P = 434.25, 150.0   # major x-tick pos (pdf pt) -> Pressure GPa
X1, X1_P = 541.90, 250.0
Y0, Y0_T = 312.33, 0.0     # major y-tick pos (pdf pt) -> Tc K
Y1, Y1_T = 158.58, 200.0

def to_pressure(x):
    return X0_P + (x - X0) * (X1_P - X0_P) / (X1 - X0)

def to_tc(y):
    return Y0_T + (y - Y0) * (Y1_T - Y0_T) / (Y1 - Y0)

markers = []
for i, d in enumerate(drawings):
    if d.get('type') != 'fs':
        continue
    r = d.get('rect')
    if r.x0 < 400:
        continue
    fill = d.get('fill')
    color = d.get('color')
    cx = r.x0 + r.width/2
    cy = r.y0 + r.height/2
    P = to_pressure(cx)
    Tc = to_tc(cy)
    markers.append((i, fill, color, round(P,1), round(Tc,1)))

markers.sort(key=lambda m: m[3])
print(f"{'idx':>4} {'fill':>18} {'color':>18} {'P(GPa)':>7} {'Tc(K)':>7}")
for m in markers:
    print(f"{m[0]:>4} {str(m[1]):>18} {str(m[2]):>18} {m[3]:>7} {m[4]:>7}")
