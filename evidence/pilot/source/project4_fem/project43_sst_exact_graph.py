"""Sparse derivative of the archived flat-plate graph residual."""

import numpy as np
from scipy.sparse import block_diag, coo_matrix, csr_matrix

from .project43_sst_bounded import BETA


def transport_volume(model, coefficients):
    """Off-diagonal transport entries, before global assembly."""
    _, gradients, _, weights, h, _ = model.geometry
    velocity, *_, diffusivity, _, _ = model.fields(coefficients, model.geometry)
    along = np.einsum("eqi,eai->eqa", velocity, gradients)
    convection = np.einsum("eq,eqa,qb->eab", weights, along, model.p)
    speed_term = 4*np.sum(velocity*velocity, axis=2)/h**2
    return [
        np.einsum("eq,eai,ebi->eab", weights*d, gradients, gradients) - convection
        + np.einsum("eq,eqa,eqb->eab", weights/np.sqrt(speed_term+(4*d/h**2)**2),
                    along, along)
        for d in diffusivity
    ]


class GraphEdges:
    def __init__(self, iteration):
        m = iteration.model
        pairs = m.triangles[:, [[0, 1], [0, 2], [1, 2]]].reshape(-1, 2)
        self.i, self.j = np.unique(np.sort(pairs, axis=1), axis=0).T
        self.model = m

    def branches(self, matrix):
        a = np.asarray(matrix[self.i, self.j]).ravel()
        b = np.asarray(matrix[self.j, self.i]).ravel()
        winner = np.argmax(np.stack((np.zeros_like(a), a, b)), axis=0)
        tied = ((a == b) & (a >= 0)) | ((a == 0) & (b <= 0)) | ((b == 0) & (a <= 0))
        return winner, tied

    def jacobian(self, state, operators):
        m = self.model
        windows = m.blocks[2:]
        graph = block_diag([csr_matrix((m.fem.size, m.fem.size)),
                            *[item[1] for item in operators]], format="csr")
        winners = [self.branches(item[0])[0] for item in operators]

        def weights(nodes, j):
            i, k = self.i, self.j
            active = winners[j] != 0
            forward = winners[j] == 1
            row, col = np.where(forward, i, k)[active], np.where(forward, k, i)[active]
            block = windows[j]
            values = state[block][row] - state[block][col]
            selector = csr_matrix((values, (row, col)), shape=(m.nb, m.nb))
            a, b = nodes[:, :, None], nodes[:, None, :]
            aa, bb = np.broadcast_arrays(a, b)
            return np.asarray(selector[aa.ravel(), bb.ravel()]).reshape(aa.shape)

        volume_weights = [weights(m.triangles, j) for j in range(2)]
        edges = m.edges[:, :2]
        edge_weights = [weights(edges, j) for j in range(2)]
        row_chunks, col_chunks, data_chunks = [], [], []

        def scatter(nodes, columns, derivative, selector, block, column):
            weighted = derivative * selector
            data = weighted.sum(axis=2) - weighted.sum(axis=1)
            keep = data != 0
            row_chunks.append(nodes[keep] + block.start)
            col_chunks.append(np.broadcast_to(columns[:, column, None], nodes.shape)[keep])
            data_chunks.append(data[keep])

        eps = 1e-28
        coefficients = state[m.dofs].astype(complex)
        for column in (*range(12), *range(15, 21)):
            coefficients[:, column] += 1j*eps
            for j, volume in enumerate(transport_volume(m, coefficients)):
                scatter(m.triangles, m.dofs, volume.imag/eps, volume_weights[j],
                        windows[j], column)
            coefficients[:, column] -= 1j*eps
        normals = m.geometry[-1]
        for column in range(6):
            flux = normals[:, column//3, None]*m.edge_velocity[:, column % 3][None, :]
            derivative = np.einsum("eq,q,qa,qb->eab", flux, m.edge_w,
                                   m.edge_scalar, m.edge_scalar)
            for j, block in enumerate(windows):
                scatter(edges, m.edge_local_dofs, derivative, edge_weights[j], block, column)
        coefficients_jac = coo_matrix((np.concatenate(data_chunks),
                                       (np.concatenate(row_chunks), np.concatenate(col_chunks))),
                                      shape=graph.shape).tocsr()
        return BETA*(graph + coefficients_jac)
