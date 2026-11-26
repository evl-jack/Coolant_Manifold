"""Analytical checks of scaled optimization without a CFD solve."""

import json
import sys
import unittest
from pathlib import Path

import numpy as np

from .project43_optimizer import DesignEvaluator, optimize, to_normalized, to_physical
from .project43_optimality import kkt_report


def analytical_response(h_mm, s_mm):
    x1 = h_mm / 4
    x2 = s_mm
    return {"converged": True,
            "power": 1 + (x1 - .2)**2 + (x2 - 1)**2,
            "side_share": .5 + .004*(x1 - x2),
            "power_gradient_mm": np.array([2*(x1-.2)/4, 2*(x2-1)]),
            "side_share_gradient_mm": np.array([.001, -.004])}


class OptimizerChecks(unittest.TestCase):
    def setUp(self):
        self.evaluator = DesignEvaluator(analytical_response, 1.)

    def test_conversion(self):
        np.testing.assert_allclose(to_physical([.375, .625]), [1.5, .625])
        np.testing.assert_allclose(to_normalized([1.5, .625]), [.375, .625])

    def test_derivatives(self):
        x = np.array([.3, .3])
        v = self.evaluator.evaluate(x)
        step = 1e-6
        for j in range(2):
            d = np.eye(2)[j]*step
            a = self.evaluator.evaluate(x+d)
            b = self.evaluator.evaluate(x-d)
            self.assertAlmostEqual(v.gradient[j], (a.objective-b.objective)/(2*step), places=8)
            np.testing.assert_allclose(v.jacobian[:,j], (a.constraints-b.constraints)/(2*step), atol=1e-8)
        np.testing.assert_allclose(v.jacobian, [[4,-4],[-4,4],[-1,-1]])

    def test_cache(self):
        first = self.evaluator.evaluate([.3, .3])
        self.assertIs(first, self.evaluator.evaluate(np.array([.3, .3])))
        self.assertEqual(self.evaluator.unique_evaluations, 1)
        self.evaluator.evaluate([.3, .30000000000000004])
        self.assertEqual(self.evaluator.unique_evaluations, 2)

    def test_invalid_response(self):
        for response in (lambda h,s: {"converged": False},
                         lambda h,s: dict(analytical_response(h,s), power=float('nan'))):
            with self.subTest(response=response), self.assertRaises(ValueError):
                DesignEvaluator(response, 1.).evaluate([.3,.3])

    def test_solution_and_kkt(self):
        result = optimize(self.evaluator)
        self.assertTrue(result.success, result.message)
        np.testing.assert_allclose(result.x, [.375, .625], atol=1e-6, rtol=0)
        value = self.evaluator.evaluate(result.x)
        self.assertAlmostEqual(value.objective, 1.17125, delta=1e-6)
        self.assertAlmostEqual(value.side_share, .499, delta=1e-8)
        report = kkt_report(self.evaluator, result.x, result)
        self.assertTrue(report["kkt_verified"], report)
        self.assertAlmostEqual(report["nominal_area_mm2"], 83., places=6)
        np.testing.assert_allclose(list(report["multipliers"].values()),
                                   [.1375, 0, .2, 0, 0, 0, 0], atol=1e-5)
        self.assertEqual(report["active_constraints"], ["side_share_lower", "nominal_area"])
        if len(sys.argv) > 1 and sys.argv[-1].endswith('.json'):
            Path(sys.argv[-1]).write_text(json.dumps(report, indent=2))

    def test_nonoptimal(self):
        self.assertFalse(kkt_report(self.evaluator, [.3,.3])["kkt_verified"])


if __name__ == "__main__":
    output = next((arg for arg in sys.argv[1:] if arg.endswith('.json')), None)
    sys.argv = sys.argv[:1]
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(OptimizerChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if output and result.wasSuccessful():
        evaluator = DesignEvaluator(analytical_response, 1.)
        solved = optimize(evaluator)
        Path(output).write_text(json.dumps({"tests_run": result.testsRun,
            "failures": len(result.failures), "errors": len(result.errors),
            "kkt": kkt_report(evaluator, solved.x, solved)}, indent=2))
    sys.exit(0 if result.wasSuccessful() else 1)
