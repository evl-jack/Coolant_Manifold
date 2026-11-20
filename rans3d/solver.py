"""Steady, incompressible, collocated finite volumes with SIMPLE pressure correction."""

from collections import deque
from dataclasses import dataclass
from hashlib import sha256

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import bicgstab, cg, splu
from scipy.spatial import cKDTree

from .mesh import totals


@dataclass(frozen=True)
class Settings:
    density: float = 998.2
    viscosity: float = 1.002e-3
    intensity: float = .05
    length_scale: float = .00056
    pseudo_cfl: float = 3.
    tolerance: float = 1e-5
    order: int = 1
    turbulence: bool = True


def linear_solve(matrix, rhs, symmetric=False):
    if len(rhs) <= 8000:
        return splu(matrix.tocsc()).solve(rhs)
    preconditioner = diags(1 / matrix.diagonal())
    method = cg if symmetric else bicgstab

    def solve_component(source):
        value, status = method(matrix, source, M=preconditioner, rtol=1e-7, atol=0, maxiter=3000)
        if status != 0:
            raise RuntimeError(f"{method.__name__} failed to converge (status {status})")
        return value

    return solve_component(rhs) if rhs.ndim == 1 else np.column_stack([solve_component(v) for v in rhs.T])


class FiniteVolume:
    def __init__(self, mesh):
        self.m = mesh
        self.n = len(mesh.volumes)
        self.b = mesh.neighbour
        self.i = len(self.b)
        self.a = mesh.owner[:self.i]
        self.patches = mesh.patches[self.i:]
        self.boundary = mesh.owner[self.i:]
        self.inlet = self.patches == "inlet"
        self.wall = self.patches == "walls"
        self.outlet = ~(self.inlet | self.wall)
        if not (self.inlet.any() and self.wall.any() and self.outlet.any()):
            raise ValueError("The mesh needs an inlet, pressure outlet and walls")
        self.area = np.linalg.norm(mesh.vectors, axis=1)
        self.normal = mesh.vectors / self.area[:, None]
        self.r = mesh.faces - mesh.centers[mesh.owner]
        self.delta = mesh.centers[self.b] - mesh.centers[self.a]
        self.distance = np.einsum("ij,ij->i", np.vstack((self.delta, self.r[self.i:])), self.normal)
        if np.any(self.distance <= 0) or not np.isfinite(self.distance).all():
            raise ValueError("Face normals must point from the owner toward the other cell or boundary")
        self.weight = np.clip(np.einsum("ij,ij->i", self.r[:self.i], self.normal[:self.i]) /
                              self.distance[:self.i], 0., 1.)
        self.nonorth = self.normal[:self.i] - self.delta / self.distance[:self.i, None]
        self.area_over_distance = self.area / self.distance
        self.rows = np.r_[np.arange(self.n), self.a, self.b]
        self.cols = np.r_[np.arange(self.n), self.b, self.a]
        self.state_gradient = self.gradient_operator(~self.outlet)
        self.pressure_gradient = self.gradient_operator(self.outlet)

    def interp(self, value):
        w = self.weight.reshape((-1,) + (1,)*(value.ndim-1))
        return (1-w)*value[self.a] + w*value[self.b]

    def sum_faces(self, internal, boundary):
        return (totals(self.a, internal, self.n) - totals(self.b, internal, self.n)
                + totals(self.boundary, boundary, self.n))

    def matrix(self, diag, upper, lower):
        return coo_matrix((np.r_[diag, upper, lower], (self.rows, self.cols)),
                          shape=(self.n, self.n)).tocsr()

    def gradient_operator(self, mask):
        cells = self.boundary[mask]
        delta, boundary_delta = self.delta, self.r[self.i:][mask]
        unit = delta / np.linalg.norm(delta, axis=1)[:, None]
        boundary_unit = boundary_delta / np.linalg.norm(boundary_delta, axis=1)[:, None]
        internal_outer = np.einsum("fi,fj->fij", unit, unit)
        lhs = (totals(self.a, internal_outer, self.n) + totals(self.b, internal_outer, self.n)
               + totals(cells, np.einsum("fi,fj->fij", boundary_unit, boundary_unit), self.n))
        inverse = np.linalg.inv(lhs)
        stencil = unit / np.linalg.norm(delta, axis=1)[:, None]
        boundary_stencil = boundary_unit / np.linalg.norm(boundary_delta, axis=1)[:, None]

        def grad(q, bc):
            values = q.reshape(self.n, -1)
            change = values[self.b] - values[self.a]
            rhs = stencil[:, :, None]*change[:, None]
            bchange = np.asarray(bc).reshape(len(mask), -1)[mask] - values[cells]
            accum = (totals(self.a, rhs, self.n) + totals(self.b, rhs, self.n)
                     + totals(cells, boundary_stencil[:, :, None]*bchange[:, None], self.n))
            value = inverse @ accum
            return value.reshape(self.n, 3, *q.shape[1:])

        return grad

    def transport(self, q, diffusion, phi, bc, dt, grad, order, source=0., sink=0.):
        volume = self.m.volumes
        def column(value):
            return value[:, None] if q.ndim == 2 else value

        positive, negative = np.maximum(phi[:self.i], 0), np.maximum(-phi[:self.i], 0)
        fb = phi[self.i:]
        face_diffusion = self.interp(diffusion)
        df = face_diffusion * self.area_over_distance[:self.i]
        db = diffusion[self.boundary] * self.area_over_distance[self.i:] * ~self.outlet
        diag = (totals(self.a, df+positive, self.n) + totals(self.b, df+negative, self.n)
                + totals(self.boundary, db+np.maximum(fb, 0), self.n) + volume/dt + sink*volume)
        rhs = column(volume)*source + column(volume/dt)*q
        rhs += totals(self.boundary, column(db+np.maximum(-fb, 0))*bc, self.n)
        face_grad = self.interp(grad)
        correction = column(face_diffusion*self.area[:self.i])*np.einsum("fi...,fi->f...", face_grad, self.nonorth)
        if order == 2:
            upwind = np.where(phi[:self.i] >= 0, self.a, self.b)
            offset = self.m.faces[:self.i] - self.m.centers[upwind]
            extrap = q[upwind] + np.einsum("fi...,fi->f...", grad[upwind], offset)
            limited = np.clip(extrap, np.minimum(q[self.a], q[self.b]), np.maximum(q[self.a], q[self.b]))
            correction -= column(phi[:self.i])*(limited-q[upwind])
        rhs += totals(self.a, correction, self.n) - totals(self.b, correction, self.n)
        return linear_solve(self.matrix(diag, -df-negative, -df-positive), rhs), diag

    def correct(self, velocity, pressure, gp, inv_diagonal, ubc, pbc):
        rf = self.interp(inv_diagonal)
        jump = pressure[self.b] - pressure[self.a] - np.einsum("ij,ij->i", self.interp(gp), self.delta)
        di = rf*self.area_over_distance[:self.i]
        db = inv_diagonal[self.boundary]*self.area_over_distance[self.i:]*self.outlet
        internal = np.einsum("ij,ij->i", self.interp(velocity), self.m.vectors[:self.i]) - di*jump
        vb = np.where(self.outlet[:, None], velocity[self.boundary], ubc)
        boundary = np.einsum("ij,ij->i", vb, self.m.vectors[self.i:])
        boundary += db*(pressure[self.boundary] - pbc + np.einsum("ij,ij->i", gp[self.boundary], self.r[self.i:]))
        diag = (totals(self.a, di, self.n) + totals(self.b, di, self.n)
                + totals(self.boundary, db, self.n))
        pc = linear_solve(self.matrix(diag, -di, -di), -self.sum_faces(internal, boundary), symmetric=True)
        internal += di*(pc[self.a]-pc[self.b])
        boundary += db*pc[self.boundary]
        gc = self.pressure_gradient(pc, np.zeros(len(boundary)))
        return velocity - inv_diagonal[:, None]*gc, pressure + .3*pc, np.r_[internal, boundary]


def sst(k, omega, grad_u, grad_k, grad_omega, wall, nu):
    strain_tensor = (grad_u + grad_u.transpose(0, 2, 1))/2
    strain = np.sqrt(2*np.sum(strain_tensor**2, axis=(1, 2)))
    viscous = 500*nu/(wall**2*omega)
    turbulent = np.sqrt(k)/(.09*omega*wall)
    cross = 2*.856*np.einsum("ij,ij->i", grad_k, grad_omega)/omega
    arg1 = np.minimum(np.maximum(turbulent, viscous), 4*.856*k/(np.maximum(cross, 1e-10)*wall**2))
    blend = np.tanh(np.minimum(arg1, 10)**4)
    f2 = np.tanh(np.minimum(np.maximum(2*turbulent, viscous), 10)**2)
    nut = .31*k/np.maximum(.31*omega, strain*f2)
    production = np.minimum(nut*strain**2, .9*k*omega)
    return nut, production, blend, (1-blend)*cross


def solve(mesh, inlet_velocity, steps, settings=Settings(), monitor=None, restart=None):
    if steps < 1 or settings.order not in (1, 2):
        raise ValueError("Use a positive iteration count and convection order 1 or 2")
    fv = FiniteVolume(mesh)
    nu = settings.viscosity/settings.density
    size = len(mesh.volumes)
    signature = sha256()
    for values in (mesh.centers, mesh.volumes, mesh.faces, mesh.vectors, mesh.owner, mesh.neighbour, mesh.patches):
        signature.update(np.ascontiguousarray(values).view(np.uint8))
    mesh_key = signature.hexdigest()
    ubc = np.zeros((len(fv.boundary), 3))
    ubc[fv.inlet] = inlet_velocity
    speed = np.linalg.norm(ubc, axis=1).max()
    if speed <= 0:
        raise ValueError("Inlet velocity must be nonzero")
    k_in = 1.5*(speed*settings.intensity)**2
    omega_in = np.sqrt(k_in)/(.09**.25*settings.length_scale)
    u, p, nut = np.zeros((size, 3)), np.zeros(size), np.zeros(size)
    k, omega = np.full(size, k_in), np.full(size, omega_in)
    if not settings.turbulence:
        k[:] = 0
        omega[:] = 0
    pbc = np.zeros(len(fv.boundary))
    kbc = np.where(fv.wall, 0., k_in)
    obc = np.full(len(fv.boundary), omega_in)
    if settings.turbulence:
        wall = cKDTree(mesh.faces[fv.i:][fv.wall]).query(mesh.centers)[0]
        np.minimum.at(wall, fv.boundary[fv.wall], fv.distance[fv.i:][fv.wall])
        obc[fv.wall] = 6*nu/(.075*fv.distance[fv.i:][fv.wall]**2)
    phi = np.r_[np.zeros(fv.i), np.einsum("ij,ij->i", ubc, mesh.vectors[fv.i:])]
    qin = -phi[fv.i:][fv.inlet].sum()
    if qin <= 0:
        raise ValueError("The specified inlet must have net inward flux")
    if restart is not None:
        if str(restart["mesh_key"]) != mesh_key:
            raise ValueError("Restart fields must use the same cell and face geometry and ordering")
        shapes = {"velocity": u.shape, "pressure_pa": p.shape, "k": k.shape,
                  "omega": omega.shape, "nut": nut.shape, "face_flux_m3_s": phi.shape}
        if any(restart[name].shape != shape or not np.isfinite(restart[name]).all()
               for name, shape in shapes.items()):
            raise ValueError("Restart fields have invalid dimensions or values")
        if settings.turbulence and (np.any(restart["k"] <= 0) or np.any(restart["omega"] <= 0)):
            raise ValueError("An SST restart requires positive k and omega from an SST solution")
        if not settings.turbulence and np.any(restart["nut"] > 0):
            raise ValueError("A laminar restart requires a laminar solution")
        u, k, omega, nut = (restart[name] for name in ("velocity", "k", "omega", "nut"))
        p = restart["pressure_pa"]/settings.density + 2*k/3
        phi = restart["face_flux_m3_s"].copy()
        fixed_flux = fv.inlet | fv.wall
        phi[fv.i:][fixed_flux] = np.einsum("ij,ij->i", ubc, mesh.vectors[fv.i:])[fixed_flux]
    dt = settings.pseudo_cfl * np.cbrt(mesh.volumes)/speed
    converged = False
    window = deque(maxlen=20)
    for iteration in range(1, steps+1):
        old_u, old_p, old_k, old_omega = u, p, k, omega
        pbc[fv.outlet] = 2*k[fv.boundary[fv.outlet]]/3
        gu = fv.state_gradient(u, ubc)
        gp = fv.pressure_gradient(p, pbc)
        if settings.turbulence:
            gk, go = fv.state_gradient(k, kbc), fv.state_gradient(omega, obc)
            nut, production, f1, cross = sst(k, omega, gu, gk, go, wall, nu)
        # Transposed viscous stress is explicit. Isotropic turbulence enters the pressure.
        grad_transpose = fv.interp(gu).transpose(0, 2, 1)
        stress_flux = fv.interp(nu+nut)[:, None]*np.einsum("fij,fi->fj", grad_transpose, mesh.vectors[:fv.i])
        stress = fv.sum_faces(stress_flux, np.zeros_like(ubc))/mesh.volumes[:, None]
        predicted, diagonal = fv.transport(u, nu+nut, phi, ubc, dt, gu, settings.order, source=stress-gp)
        u, p, phi = fv.correct(predicted, p, gp, mesh.volumes/diagonal, ubc, pbc)
        if settings.turbulence:
            k = fv.transport(k, nu+(.85*f1+1-f1)*nut, phi, kbc, dt, gk, settings.order,
                             source=production, sink=.09*omega)[0].clip(1e-12)
            beta = .075*f1 + .0828*(1-f1)
            gamma = (5/9)*f1 + .44*(1-f1)
            omega = fv.transport(omega, nu+(.5*f1+.856*(1-f1))*nut, phi, obc, dt, go, settings.order,
                                 source=gamma*production/np.maximum(nut, 1e-20)+np.maximum(cross, 0),
                                 sink=beta*omega + np.maximum(-cross, 0)/omega)[0].clip(1e-8)
        change = max(np.max(abs(u-old_u))/speed, np.max(abs(p-old_p))/speed**2)
        if settings.turbulence:
            change = max(change, np.max(abs(k-old_k))/max(k_in, 1e-12),
                         np.max(abs(omega-old_omega))/max(omega_in, np.max(abs(omega))))
        continuity = np.sum(abs(fv.sum_faces(phi[:fv.i], phi[fv.i:])))/qin
        flows = {name: float(phi[fv.i:][fv.patches == name].sum()) for name in np.unique(fv.patches)}
        balance = abs(sum(flows.values())) / qin
        pressure_pa = settings.density*(p-2*k/3)
        pressure_drop = np.average(pressure_pa[fv.boundary[fv.inlet]], weights=fv.area[fv.i:][fv.inlet])
        window.append([pressure_drop, *(value for name, value in flows.items() if name.endswith("outlet"))])
        recent = np.asarray(window)
        pressure_drift = np.ptp(recent[:, 0])/max(abs(np.mean(recent[:, 0])), 1e-12)
        flow_drift = np.ptp(recent[:, 1:], axis=0).max()/qin
        metrics = dict(iteration=iteration, field_change=float(change), continuity=float(continuity),
                       mass_imbalance=float(balance), pressure_drop_pa=float(pressure_drop),
                       monitor_samples=len(window), outlet_flow_drift=float(flow_drift),
                       pressure_drop_drift=float(pressure_drift),
                       flow_m3_s=flows, reverse_flow_m3_s={name: float(-np.minimum(
                           phi[fv.i:][fv.patches == name], 0).sum()) for name in flows if name.endswith("outlet")})
        if not np.isfinite(change+continuity+balance):
            raise FloatingPointError(f"Non-finite solver state at iteration {iteration}")
        if monitor:
            monitor(metrics)
        if (max(change, continuity, balance) < settings.tolerance
                and len(window) == window.maxlen and max(pressure_drift, flow_drift) <= 1e-3):
            converged = True
            break
    if settings.turbulence:
        nut = sst(k, omega, fv.state_gradient(u, ubc), fv.state_gradient(k, kbc),
                  fv.state_gradient(omega, obc), wall, nu)[0]
    return dict(metrics, converged=converged), dict(centers=mesh.centers, velocity=u, pressure_pa=pressure_pa,
                                                    k=k, omega=omega, nut=nut, face_flux_m3_s=phi, mesh_key=mesh_key)
