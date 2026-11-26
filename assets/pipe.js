// Conforming inner and outer cylindrical lofts, in physical millimeters.
const TAU = 2 * Math.PI;
const STEPS = 24;

function quad(triangles, a, b, c, d) {
    triangles.push([a, b, c], [a, c, d]);
}

function skin(sections, thickness) {
    const main = 2 + thickness;
    const [floor, left, right] = sections[0];
    const center = (left + right) / 2;
    const aperture = (right - left) / 2 + thickness;
    if (aperture > main + 1e-9) throw Error("Branch aperture exceeds main bore");
    const alpha = Math.asin(Math.min(1, aperture / main));
    const angles = [], edgeZ = [];

    for (let i = 0; i < STEPS; i++) angles.push(-Math.PI + (Math.PI - alpha) * i / STEPS);
    for (let i = 0; i <= STEPS; i++) {
        const z = aperture * Math.sin(-Math.PI / 2 + Math.PI * i / STEPS);
        edgeZ.push(z);
        angles.push(Math.asin(z / main));
    }
    for (let i = 1; i < STEPS; i++) angles.push(alpha + (Math.PI - alpha) * i / STEPS);

    const sectionsX = angles.map((theta, index) => {
        const y = 2 + main * Math.cos(theta), z = main * Math.sin(theta);
        const i = index - STEPS;
        const dx = i > 0 && i < STEPS ? Math.sqrt(aperture ** 2 - edgeZ[i] ** 2) : 0;
        return [[0, y, z], [center - dx, y, z], [center + dx, y, z], [20, y, z]];
    });
    const mainFaces = [];
    sectionsX.forEach((row, i) => {
        const next = sectionsX[(i + 1) % sectionsX.length];
        quad(mainFaces, row[0], next[0], next[1], row[1]);
        quad(mainFaces, row[2], next[2], next[3], row[3]);
    });

    const apertureRing = [
        ...sectionsX.slice(STEPS, 2 * STEPS + 1).map(row => row[2]),
        ...sectionsX.slice(STEPS + 1, 2 * STEPS).reverse().map(row => row[1]),
    ];
    const angular = apertureRing.map(point => [(point[0] - center) / aperture, point[2] / aperture]);
    const rings = [apertureRing];
    for (const [y, lo, hi] of sections) {
        const x = (lo + hi) / 2, radius = (hi - lo) / 2 + thickness;
        const displayY = y + thickness * Math.max(0, (6 - y) / (6 - floor));
        rings.push(angular.map(([u, v]) => [x + radius * u, displayY, radius * v]));
    }
    const branchFaces = [];
    for (let j = 0; j < rings.length - 1; j++) {
        const a = rings[j], b = rings[j + 1];
        a.forEach((point, i) => {
            const k = (i + 1) % a.length;
            quad(branchFaces, point, a[k], b[k], b[i]);
        });
    }

    return {
        mainFaces, branchFaces,
        ports: [sectionsX.map(row => row[0]), sectionsX.map(row => row[3]), rings.at(-1)],
    };
}

const key = p => p.map(x => Math.round(x * 1e8)).join("/");
const normal = (a, b, c) => [
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
    (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]),
];

export function buildPipe(sections, thickness, audit = false) {
    const inner = skin(sections, 0), outer = skin(sections, thickness);
    const rims = [];
    outer.ports.forEach((ring, j) => {
        const hole = inner.ports[j];
        ring.forEach((point, i) => {
            const k = (i + 1) % ring.length;
            quad(rims, point, ring[k], hole[k], hole[i]);
        });
    });
    const data = [], allTriangles = [];
    function append(faces, invert, material) {
        const normals = new Map();
        const triangles = [];
        for (const face of faces) {
            const n = normal(...face);
            const len = Math.hypot(...n);
            if (len < 1e-10) continue;
            triangles.push(face);
            for (const point of face) {
                const id = key(point);
                const existing = normals.get(id) ?? [0, 0, 0];
                n.forEach((value, i) => { existing[i] += invert ? -value : value; });
                normals.set(id, existing);
            }
        }
        for (const face of triangles) {
            for (const point of face) {
                const n = normals.get(key(point));
                const len = Math.hypot(...n) || 1;
                data.push(point[0] - 10, point[1] - 3.5, point[2], ...n.map(x => x / len), material, 0);
            }
        }
        if (audit) allTriangles.push(...triangles);
    }

    append(outer.mainFaces, false, 0);
    append(outer.branchFaces, true, 0);
    append(inner.mainFaces, true, 1);
    append(inner.branchFaces, false, 1);
    append(rims, false, 2);
    if (audit) {
        const edges = new Map(), adjacency = new Map();
        for (const face of allTriangles) {
            const ids = face.map(key);
            ids.forEach((a, i) => {
                const b = ids[(i + 1) % 3];
                const edge = [a, b].sort().join(":");
                edges.set(edge, (edges.get(edge) ?? 0) + 1);
                if (!adjacency.has(a)) adjacency.set(a, []);
                adjacency.get(a).push(b);
            });
        }
        const openEdges = [...edges.values()].filter(n => n !== 2).length;
        const pending = [adjacency.keys().next().value], seen = new Set(pending);
        while (pending.length) for (const b of adjacency.get(pending.pop())) {
            if (!seen.has(b)) { seen.add(b); pending.push(b); }
        }
        return { openEdges, connected: seen.size === adjacency.size, triangles: allTriangles.length, ports: inner.ports };
    }
    return new Float32Array(data);
}
