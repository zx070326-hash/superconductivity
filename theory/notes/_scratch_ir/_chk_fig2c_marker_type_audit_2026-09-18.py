import fitz
doc = fitz.open('../../../drozdov2015.pdf') if False else None
import os
# locate pdf
cands = []
for root, dirs, files in os.walk('D:/physics/superconductivity'):
    for f in files:
        if f == 'drozdov2015.pdf':
            cands.append(os.path.join(root, f))
print("pdf candidates:", cands)
path = cands[0]
doc = fitz.open(path)
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

def item_cmds(d):
    return [it[0] for it in d.get('items', [])]

print("=== RED (stroke-color) candidates: item command signatures ===")
reds = []
for d in drawings:
    r = d.get('rect')
    if r is None or r.x0 < 400:
        continue
    color = d.get('color')
    if color is None:
        continue
    if color[0] > 0.6 and color[1] < 0.4 and color[2] < 0.4:
        cmds = item_cmds(d)
        cx, cy = r.x0 + r.width/2, r.y0 + r.height/2
        reds.append((round(to_pressure(cx),2), round(to_tc(cy),2), d.get('type'), tuple(cmds), round(r.width,2), round(r.height,2), d.get('fill')))

for m in sorted(reds):
    print(m)
print("n_red=", len(reds))
