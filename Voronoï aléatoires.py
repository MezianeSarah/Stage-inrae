import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from collections import defaultdict
import os


def voronoi_finite_polygons_2d(vor, radius=None):
    if vor.points.shape[1] != 2:
        raise ValueError("Supports 2D diagrams only.")
    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max() * 2
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))
    for p1, region_idx in enumerate(vor.point_region):
        vertices = vor.regions[region_idx]
        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]
        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue
            t = vor.points[p2] - vor.points[p1]
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])
            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius
            new_vertices.append(far_point.tolist())
            new_region.append(len(new_vertices) - 1)
        new_regions.append(new_region)
    return new_regions, np.array(new_vertices)

def is_segment_inside(p1, p2, x_min, x_max, y_min, y_max):
    return (
        x_min <= p1[0] <= x_max and y_min <= p1[1] <= y_max and
        x_min <= p2[0] <= x_max and y_min <= p2[1] <= y_max
    )

def is_horizontal(p1, p2):
    return abs(p1[1] - p2[1]) < abs(p1[0] - p2[0])

def segments_intersect(x1, y1, x2, y2, x3, y3, x4, y4):
    def ccw(A, B, C):
        return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
    A, B = (x1, y1), (x2, y2)
    C, D = (x3, y3), (x4, y4)
    intersect = ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)
    if not intersect:
        return False, None
    denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
    if denom == 0:
        return False, None
    px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / denom
    return True, (px, py)

def point_outside_box(x1, y1, x2, y2, xmin, ymin, xmax, ymax):
    for (x, y) in [(x1, y1), (x2, y2)]:
        if not (xmin <= x <= xmax and ymin <= y <= ymax):
            return (x, y)
    return None

def peridynamic_output3(points, bounds, output_file="output_peridynamic.txt"):
    vor = Voronoi(points)
    xmin, xmax, ymin, ymax = bounds
    xmin += 0.2 * xmax
    ymin += 0.2 * ymax
    xmax -= 0.2 * xmax
    ymax -= 0.2 * ymax

    with open(output_file, 'w') as f:
        compt = 0
        for ridge in vor.ridge_vertices:
            if -1 in ridge or len(ridge) != 2:
                continue
            v0, v1 = vor.vertices[ridge]
            epaisseur = 0.002
            if (xmin <= v0[0] <= xmax and ymin <= v0[1] <= ymax and
                xmin <= v1[0] <= xmax and ymin <= v1[1] <= ymax):
                f.write(f"{v0[0]}\t{v0[1]}\t! coordonnées du point 1 du segment {compt}\n")
                f.write(f"{v1[0]}\t{v1[1]}\t! coordonnées du point 2 du segment {compt}\n")
                f.write(f"{epaisseur}\t! épaisseur du segment {compt}\n\n")
                compt += 1
            else:
                for segment_boite in [[xmin,ymin,xmax,ymin],[xmax,ymin,xmax,ymax],[xmin,ymax,xmax,ymax],[xmin,ymin,xmin,ymax]]:
                    x3, y3, x4, y4 = segment_boite
                    intersect, point = segments_intersect(v0[0], v0[1], v1[0], v1[1], x3, y3, x4, y4)
                    if not intersect:
                        continue
                    out = point_outside_box(v0[0], v0[1], v1[0], v1[1], xmin, ymin, xmax, ymax)
                    if out == (v0[0], v0[1]):
                        v0 = np.array(point)
                    elif out == (v1[0], v1[1]):
                        v1 = np.array(point)
                    if (xmin <= v0[0] <= xmax and ymin <= v0[1] <= ymax and
                        xmin <= v1[0] <= xmax and ymin <= v1[1] <= ymax):
                        f.write(f"{v0[0]}\t{v0[1]}\t! coordonnées du point 1 du segment {compt}\n")
                        f.write(f"{v1[0]}\t{v1[1]}\t! coordonnées du point 2 du segment {compt}\n")
                        f.write(f"{epaisseur}\t! épaisseur du segment {compt}\n\n")
                        compt += 1

        for (x0, y0), (x1, y1) in [[(xmin, ymax), (xmax, ymax)], [(xmax, ymin), (xmin, ymin)]]:
            epaisseur = 0.004
            f.write(f"{x0}\t{y0}\t! coordonnées du point 1 du segment {compt}\n")
            f.write(f"{x1}\t{y1}\t! coordonnées du point 2 du segment {compt}\n")
            f.write(f"{epaisseur}\t! épaisseur du segment {compt}\n\n")
            compt += 1

    with open(output_file, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f"{compt}\t! nombre de segments\n\n" + content)


base_x_range = (0, 45)
base_y_range = (0, 60)
base_cols_range = (50, 60)
y_step_range = (0.6, 1.2)
threshold = 0.3
num_models = 20
output_dir = "voronoi_models_scaled"
os.makedirs(output_dir, exist_ok=True)

for model_id in range(1, num_models + 1):
    scale_factor = np.random.uniform(0.6, 1.4)
    x_min, x_max = [val * scale_factor for val in base_x_range]
    y_min, y_max = [val * scale_factor for val in base_y_range]
    x_range = (x_min, x_max)
    cols_range = (
        int(base_cols_range[0] * scale_factor),
        int(base_cols_range[1] * scale_factor)
    )

    points = []
    y = y_min
    while y < y_max:
        n_cols = np.random.randint(*cols_range)
        offset_x = np.random.uniform(-0.2, 0.9)
        xs = np.linspace(x_range[0], x_range[1], n_cols) + offset_x
        for x in xs:
            jitter_x = np.random.uniform(-0.1, 0.1)
            jitter_y = np.random.uniform(-0.01, 0.01)
            points.append([x + jitter_x, y + jitter_y])
        y += np.random.uniform(*y_step_range)

    points = np.array(points)
    min_y = np.min(points[:, 1])
    max_y = np.max(points[:, 1])

    vor = Voronoi(points)
    regions, vertices = voronoi_finite_polygons_2d(vor)

    segments = []
    highlight_segments = []

    for i, region in enumerate(regions):
        point = vor.points[i]
        is_top = np.abs(point[1] - max_y) < threshold
        is_bottom = np.abs(point[1] - min_y) < threshold

        for j in range(len(region)):
            v1 = vertices[region[j]]
            v2 = vertices[(region[(j + 1) % len(region)])]
            if is_segment_inside(v1, v2, x_min, x_max, y_min, y_max):
                segments.append((v1, v2))
                if (is_top or is_bottom) and is_horizontal(v1, v2):
                    highlight_segments.append((tuple(v1), tuple(v2)))

    segments.extend(highlight_segments)

    img_filename = os.path.join(output_dir, f"model_{model_id:03d}.png")
    fig, ax = plt.subplots(figsize=(10, 10))
    for p1, p2 in segments:
        if (tuple(p1), tuple(p2)) in highlight_segments or (tuple(p2), tuple(p1)) in highlight_segments:
            continue
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-', linewidth=1.0)

    line_groups = defaultdict(list)
    for p1, p2 in highlight_segments:
        y_avg = round((p1[1] + p2[1]) / 2, 3)
        line_groups[y_avg].extend([p1, p2])

    for y_val, pts in line_groups.items():
        pts = np.array(pts)
        pts = pts[np.argsort(pts[:, 0])]
        y_corrected = np.full_like(pts[:, 1], y_val)
        ax.plot(pts[:, 0], y_corrected, 'k-', linewidth=12.0)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(img_filename)
    plt.close()

  
    peridynamic_filename = os.path.join(output_dir, f"model_{model_id:03d}_peridynamic.txt")
    peridynamic_output3(points, (x_min, x_max, y_min, y_max), output_file=peridynamic_filename)

print(f"✅ {num_models} modèles générés avec images et fichiers Peridynamic uniquement dans : '{output_dir}'")
