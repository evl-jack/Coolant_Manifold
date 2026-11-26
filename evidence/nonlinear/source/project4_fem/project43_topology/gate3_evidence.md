# Gate 3 fitted-wall reference verification

Sol job 2395603 solved the exact sharp and obstructed 83 mm² geometries at three fitted resolution levels. The obstructed fixture is the fixed geometry described in [the diagnostic protocol](gate3_protocol.md).

| Mesh spacing (mm) | Sharp triangles | Sharp side flow | Sharp power (µW/m) | Obstructed triangles | Obstructed side flow | Obstructed power (µW/m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 2,656 | 4.039488% | 9.280431 | 4,114 | 6.785281% | 11.487428 |
| 0.125 | 10,624 | 4.031835% | 9.276688 | 16,456 | 6.775616% | 11.481570 |
| 0.0625 | 42,496 | 4.026955% | 9.275452 | 65,918 | 6.771887% | 11.482819 |

| Case | First split change (pp) | Second split change (pp) | First relative power change | Second relative power change |
|---|---:|---:|---:|---:|
| Sharp | 0.007653 | 0.004879 | 0.04035% | 0.01332% |
| Obstructed | 0.009666 | 0.003728 | 0.05102% | 0.01087% |

Both changes in each fixture meet the prospective 0.02 pp and 0.5% limits. The obstructed power changed nonmonotonically. The largest relative mass imbalance was 2.91e-14 and the largest nonlinear residual was 4.80e-12. The mesh changes provide convergence evidence for the fitted reference ladder. The immersed comparison and layer resolution checks determine whether a Brinkman configuration qualifies.

The three-resistance, three-mesh binary screen was submitted as Sol array 2395875 after the fitted reference check. The fitted and screen reports are under `/scratch/jtfoste6/p43-topology-study-20261125/`.

## Uniform binary Brinkman screen

All nine mesh and resistance combinations completed. Every combination failed at least one prospective physical or leakage requirement on each fixture. The finest results below use axis spacing 0.0625 mm. All leakage percentages refer to inlet flux.

| Fixture | Maximum resistance kg/(m³·s) | Fitted split error (pp) | Relative power error | Gross ingress | Gross egress | Largest absolute cut flux | Solid RMS / inlet mean |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sharp | 100,000 | 1.85593 | 4.49995% | 3.00370% | 2.98534% | 1.35684% | 2.78741% |
| Obstructed | 100,000 | 1.08048 | 11.38402% | 10.14853% | 10.12089% | 1.42619% | 5.15600% |
| Sharp | 1,000,000 | 0.58041 | 1.46082% | 0.49114% | 0.48404% | 0.18878% | 0.52307% |
| Obstructed | 1,000,000 | 0.42812 | 3.78116% | 2.40194% | 2.38849% | 0.18879% | 1.04923% |
| Sharp | 10,000,000 | 0.15227 | 0.39495% | 0.06865% | 0.06708% | 0.02147% | 0.09725% |
| Obstructed | 10,000,000 | 0.11730 | 1.08358% | 0.41385% | 0.41067% | 0.02191% | 0.19182% |
| Limit | | 0.02 | 0.5% | 0.10% | 0.10% | 0.02% | 1.0% |

At the largest resistance, successive immersed split changes were 0.02989 and 0.02733 pp for the sharp fixture and 0.02247 and 0.02214 pp for the obstructed fixture. All four exceed 0.02 pp. The finest run at this resistance used peak RSS 2.95 million KiB and 7 minutes 7 seconds for both fixtures. The fine axis spacing is about 6.24 times the nominal penetration length at this resistance.

Sol array 2396361 was submitted for two fixed-fixture, interface-refined runs at the largest screened resistance. It applies the prospective half-penetration-length spacing criterion and a 140,000-triangle resource cap. The two cases are a diagnostic confirmation effort. There is no selected production resistance or refinement envelope.

The sharp and obstructed refined meshes contain 112,104 and 134,040 triangles, respectively, or 506,326 and 605,118 state DOFs. Both resolve the fixed material-interface normal and tangential spacings to 0.00390625 mm, below the 0.005005 mm requirement. The certified bands extend 0.02002 mm on each side of the fixed interfaces. The separate port-band maximum diameter is 0.04420 mm against 0.0625 mm. Uniform axis spacing of 0.005005 mm over the entire 60 mm² design region alone would imply roughly 4.8 million structured triangles before accounting for the remaining domain or grading. The local meshes do not establish costs or resolution for moving production interfaces.

| Locally resolved fixture | Fitted split error (pp) | Relative power error | Gross ingress | Gross egress | Largest absolute cut flux | Solid RMS / inlet mean | Wall solve time (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sharp | 0.168779 | 0.46204% | 0.071786% | 0.071731% | 0.020395% | 0.09155% | 287 |
| Obstructed | 0.129576 | 1.22345% | 0.435522% | 0.435419% | 0.021200% | 0.18442% | 381 |

Both completed solves fail the approved physical and through-solid criteria. The obstructed fixture also fails the power and gross exchange criteria. This confirms a discrepancy persists after resolving the tested fixed interfaces at 10,000,000 kg/(m³·s). Further work requires a reviewed production configuration or a separately reviewed physical-method change.
