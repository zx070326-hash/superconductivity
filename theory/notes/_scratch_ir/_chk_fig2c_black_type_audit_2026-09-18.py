import fitz
doc = fitz.open('theory/notes/_scratch_ir/drozdov2015.pdf')
page = doc[8]
drawings = page.get_drawings()

X0, X0_P = 434.25, 150.0
X1, X1_P = 541.90, 250.0
Y0, Y0_T = 312.33, 0.0
Y1, Y1_T = 158.58, 200.0
def to_pressure(x): return X0_P + (x - X0) * (X1_P - X0_P) / (X1 - X0)
def to_tc(y): return Y0_T + (y - Y0) * (Y1_T - Y0_T) / (Y1 - Y0)

blacks = []
for d in drawings:
    r = d.get('rect')
    if r is None or r.x0 < 400:
        continue
    color = d.get('color')
    fill = d.get('fill')
    if color is None:
        continue
    if max(color) < 0.25 and fill is not None and max(fill) < 0.25:
        cmds = tuple(it[0] for it in d.get('items', []))
        cx, cy = r.x0 + r.width/2, r.y0 + r.height/2
        blacks.append((round(to_pressure(cx),2), round(to_tc(cy),2), d.get('type'), cmds, round(r.width,2), round(r.height,2)))

for m in sorted(blacks):
    print(m)
print("n=", len(blacks))
