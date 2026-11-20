"""Geometry for the 3D finite-volume solver. All coordinates are in meters."""

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from sol_rans_3d.repair_hex_mesh import read_list


@dataclass
class Mesh:
    centers: np.ndarray
    volumes: np.ndarray
    faces: np.ndarray
    vectors: np.ndarray
    owner: np.ndarray
    neighbour: np.ndarray
    patches: np.ndarray


def totals(cells, values, size):
    if values.ndim == 1:
        return np.bincount(cells, values, minlength=size)
    flat = values.reshape(len(cells), -1)
    return np.stack([np.bincount(cells, column, minlength=size) for column in flat.T],
                    axis=1).reshape((size, *values.shape[1:]))


def from_polymesh(directory):
    """Read binary compact faces and little-endian 32-bit labels from the Sol mesh."""
    directory = Path(directory)
    raw = (directory / "points").read_bytes()
    if b'LSB;label=32;scalar=64' not in raw[:1024] or not re.search(rb"format\s+binary;", raw[:1024]):
        raise ValueError("Expected a binary LSB;label=32;scalar=64 polyMesh")
    points = read_list(raw, "<f8", width=3)[0].reshape(-1, 3)
    raw = (directory / "faces").read_bytes()
    if b'faceCompactList' not in raw[:1024]:
        raise ValueError("Expected binary compact faces")
    offsets, end = read_list(raw, "<i4")
    verts, _ = read_list(raw, "<i4", end)
    owner, _ = read_list((directory / "owner").read_bytes(), "<i4")
    neighbour, _ = read_list((directory / "neighbour").read_bytes(), "<i4")
    patches = np.full(len(owner), "", dtype="U11")
    boundary = (directory / "boundary").read_text(encoding="utf-8")
    recognized = {"inlet", "main_outlet", "side_outlet", "walls"}
    for name, body in re.findall(r"(?m)^\s*(\w+)\s*\{([^{}]*)\}", boundary):
        if name == "FoamFile":
            continue
        start = int(re.search(r"\bstartFace\s+(\d+)", body)[1])
        count = int(re.search(r"\bnFaces\s+(\d+)", body)[1])
        if count and name not in recognized:
            raise ValueError(f"Unsupported patch {name}")
        patches[start:start+count] = name
    if not np.all(patches[len(neighbour):] != "") or np.any(patches[:len(neighbour)] != ""):
        raise ValueError("The boundary patches do not match the face addressing")

    # Decompose each polygon into triangles around its average vertex position.
    widths = np.diff(offsets)
    slots = np.arange(widths.max())
    valid = slots < widths[:, None]
    indices = verts[offsets[:-1, None] + np.where(valid, slots, 0)]
    vertices = points[indices]
    center = (vertices * valid[..., None]).sum(1) / widths[:, None]
    following = np.roll(vertices, -1, axis=1)
    triangles = np.cross(vertices - center[:, None], following - center[:, None]) / 2
    vectors = triangles.sum(1)
    weight = np.einsum("fij,fj->fi", triangles, vectors)
    faces = (weight[..., None]*(vertices + following + center[:, None])).sum(1) / (3*weight.sum(1)[:, None])

    # Sum signed tetrahedron volumes and first moments separately for each side.
    n = int(owner.max()) + 1
    internal = len(neighbour)
    cells = np.r_[owner, neighbour]
    reference = totals(cells, np.concatenate((faces, faces[:internal])), n)
    reference /= np.bincount(cells, minlength=n)[:, None]

    def pyramids(c, sign, count):
        base = vertices[:count]
        origin = reference[c]
        tetra = sign / 3 * np.einsum("fij,fij->fi", triangles[:count], base - origin[:, None])
        moment = (base + following[:count] + center[:count, None] + origin[:, None]) / 4
        return tetra.sum(1), np.einsum("fi,fij->fj", tetra, moment)

    v, m = pyramids(owner, 1, len(owner))
    vn, mn = pyramids(neighbour, -1, internal)
    volumes = totals(cells, np.r_[v, vn], n)
    if not np.isfinite(volumes).all() or np.any(volumes <= 0):
        raise ValueError("Invalid oriented mesh volumes")
    centers = totals(cells, np.concatenate((m, mn)), n) / volumes[:, None]
    graph = coo_matrix((np.ones(internal, bool), (owner[:internal], neighbour)), shape=(n, n))
    if connected_components(graph, directed=False, return_labels=False) != 1:
        raise ValueError("Expected one connected fluid region")
    return Mesh(centers, volumes, faces, vectors, owner, neighbour, patches)


def duct(nx, ny, length=.020, width=.0035):
    """A uniform square duct with inlet and outlet normals along the x axis."""
    dims = (nx, ny, ny)
    spacing = np.array([length/nx, width/ny, width/ny])
    index = np.arange(np.prod(dims)).reshape(dims)
    centers = (np.indices(dims).reshape(3, -1).T + .5) * spacing
    volume = np.prod(spacing)
    internal, boundary = [], []
    for axis in range(3):
        line = np.moveaxis(index, axis, 0)
        area = np.eye(3)[axis] * volume / spacing[axis]
        a, b = line[:-1].ravel(), line[1:].ravel()
        internal.append((a, b, (centers[a]+centers[b])/2, len(a)*[area], [""]*len(a)))
        for cells, sign in ((line[0].ravel(), -1), (line[-1].ravel(), 1)):
            face = centers[cells] + sign * np.eye(3)[axis] * spacing[axis]/2
            name = ("inlet" if sign < 0 else "main_outlet") if axis == 0 else "walls"
            boundary.append((cells, [], face, len(cells)*[sign*area], [name]*len(cells)))
    a, b, faces, vectors, patches = (np.concatenate(values) for values in zip(*(internal+boundary)))
    return Mesh(centers, np.full(len(centers), volume), faces, vectors, a, b.astype(int), patches)
