import fitz
doc = fitz.open('drozdov2015.pdf')
page = doc[8]
drawings = page.get_drawings()

X0, X0_P = 434.25, 150.0
X1, X1_P = 541.90, 250.0
Y0, Y0_T = 312.33, 0.0
Y1, Y1_T = 158.58, 200.0

def to_pressure(x):
    return X0_P + (x - X0) * (X1_P - X0_P) / (X1 - X0)

def to_tc(y):
    return Y0_T + (y - Y0) * (Y1_T - Y0_T) / (Y1 - Y0)

blacks = []
for d in drawings:
    r = d.get('rect')
    if r is None or r.x0 < 400:
        continue
    color = d.get('color')
    fill = d.get('fill')
    if color is None:
        continue
    # black filled squares (H3S): color near (0,0,0), fill also near (0,0,0)
    if max(color) < 0.25 and fill is not None and max(fill) < 0.25:
        cx = r.x0 + r.width / 2
        cy = r.y0 + r.height / 2
        P = to_pressure(cx)
        Tc = to_tc(cy)
        blacks.append((round(P, 2), round(Tc, 2), d.get('type'), round(r.width, 2)))

blacks.sort()
print(f"{'P(GPa)':>8} {'Tc(K)':>8} {'type':>5} {'w':>6}")
for m in blacks:
    print(f"{m[0]:>8} {m[1]:>8} {m[2]:>5} {m[3]:>6}")
print("n=", len(blacks))

# look for near-duplicate pressures (within 1 GPa) among H3S markers
print("\nnear-duplicate P pairs (<1.5 GPa apart):")
for i in range(len(blacks)):
    for j in range(i+1, len(blacks)):
        if abs(blacks[i][0]-blacks[j][0]) < 1.5:
            print(blacks[i], blacks[j])
