import re

doc = open("/home/ys/feedback/uidreczz28dlJndpQsN/artifacts/report.html").read()
svgs = re.findall(r"<svg.*?</svg>", doc, re.S)
print("svg count:", len(svgs))
issues = []
for k, s in enumerate(svgs):
    W, H = map(float, re.search(r'viewBox="0 0 (\d+) (\d+)"', s).groups())
    for m in re.finditer(r'<(circle|rect)\b[^>]*?\b(c?x|x)="([-\d.]+)"[^>]*?\b(c?y|y)="([-\d.]+)"', s):
        x, y = float(m.group(3)), float(m.group(5))
        if not (-1 <= x <= W + 1 and -1 <= y <= H + 1):
            issues.append((k, m.group(1), x, y))
    for m in re.finditer(r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*>(.*?)</text>', s, re.S):
        x, y = float(m.group(1)), float(m.group(2))
        if not (-1 <= x <= W + 1 and -1 <= y <= H + 1):
            issues.append((k, "text-out", x, y, m.group(3)[:20]))
    for m in re.finditer(r'points="([^"]+)"', s):
        for pair in m.group(1).split():
            x, y = map(float, pair.split(","))
            if not (-1 <= x <= W + 1 and -1 <= y <= H + 1):
                issues.append((k, "poly-pt", x, y))
print("out-of-bounds issues:", issues if issues else "none")
for k, s in enumerate(svgs):
    labels = re.findall(r'class="end-label">(.*?)</text>', s)
    ys = sorted(float(y) for _, y, _ in
                re.findall(r'<text x="(\d+\.\d+)" y="(\d+\.\d+)" class="end-label">(.*?)</text>', s))
    collisions = [ys[i] for i in range(1, len(ys)) if ys[i] - ys[i - 1] < 14]
    print("chart", k, "end-labels=", len(labels),
          "y-collisions=", collisions if collisions else "none")
print("table rows (<tr>):", doc.count("<tr>"), "tables:", doc.count("<table>"))
print("tooltip listeners:", doc.count("addEventListener"))
