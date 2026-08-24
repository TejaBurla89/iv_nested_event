import tempfile
import unittest
from pathlib import Path

from weak_iv_finite_strata import project, read_rows


class FiniteStrataProjectionTests(unittest.TestCase):
    def test_positive_first_stage_intervals_pass_orientation_diagnostic(self):
        rows = [
            {
                "stratum": "x0",
                "reference_weight": 1.0,
                "rho_pre": 0.40,
                "se_rho_pre": 0.04,
                "pi_pre": 0.20,
                "se_pi_pre": 0.01,
                "rho_post": 0.60,
                "se_rho_post": 0.05,
                "pi_post": 0.20,
                "se_pi_post": 0.01,
            }
        ]
        result = project(rows, alpha=0.05)
        self.assertEqual(
            result["first_stage_orientation"]["point_estimate_status"],
            "positive",
        )
        self.assertEqual(
            result["first_stage_orientation"]["simultaneous_interval_assessment"],
            "consistent_with_declared_positive_orientation",
        )
        self.assertFalse(
            result["first_stage_orientation"]["causal_interpretation_stop"]
        )

    def test_weak_required_first_stage_returns_all_real(self):
        rows = [
            {
                "stratum": "x0",
                "reference_weight": 1.0,
                "rho_pre": 0.40,
                "se_rho_pre": 0.04,
                "pi_pre": 0.20,
                "se_pi_pre": 0.015,
                "rho_post": 0.03,
                "se_rho_post": 0.05,
                "pi_post": 0.01,
                "se_pi_post": 0.02,
            }
        ]
        result = project(rows, alpha=0.05)
        self.assertTrue(result["confidence_set"]["all_real"])
        self.assertFalse(
            result["first_stage_orientation"]["causal_interpretation_stop"]
        )
        self.assertEqual(
            result["first_stage_orientation"]["simultaneous_interval_assessment"],
            "inconclusive_because_an_interval_spans_zero",
        )

    def test_mixed_first_stage_signs_trigger_causal_stop(self):
        rows = [
            {
                "stratum": "x0",
                "reference_weight": 1.0,
                "rho_pre": 0.40,
                "se_rho_pre": 0.04,
                "pi_pre": 0.20,
                "se_pi_pre": 0.015,
                "rho_post": -0.60,
                "se_rho_post": 0.05,
                "pi_post": -0.20,
                "se_pi_post": 0.015,
            }
        ]
        result = project(rows, alpha=0.05)
        self.assertEqual(
            result["first_stage_orientation"]["point_estimate_status"],
            "mixed_or_zero",
        )
        self.assertTrue(
            result["first_stage_orientation"]["causal_interpretation_stop"]
        )
        self.assertEqual(
            result["first_stage_orientation"]["simultaneous_interval_assessment"],
            "contradicts_declared_positive_orientation",
        )

    def test_noisy_negative_point_sign_is_warning_not_causal_stop(self):
        rows = [
            {
                "stratum": "x0",
                "reference_weight": 1.0,
                "rho_pre": 0.40,
                "se_rho_pre": 0.04,
                "pi_pre": 0.20,
                "se_pi_pre": 0.015,
                "rho_post": -0.03,
                "se_rho_post": 0.05,
                "pi_post": -0.01,
                "se_pi_post": 0.10,
            }
        ]
        result = project(rows, alpha=0.05)
        self.assertEqual(
            result["first_stage_orientation"]["point_estimate_status"],
            "mixed_or_zero",
        )
        self.assertEqual(
            result["first_stage_orientation"]["simultaneous_interval_assessment"],
            "inconclusive_because_an_interval_spans_zero",
        )
        self.assertFalse(
            result["first_stage_orientation"]["causal_interpretation_stop"]
        )

    def test_zero_weight_row_is_rejected(self):
        contents = (
            "stratum,reference_weight,rho_pre,se_rho_pre,pi_pre,se_pi_pre,"
            "rho_post,se_rho_post,pi_post,se_pi_post\n"
            "x0,1,0.4,0.04,0.2,0.015,0.6,0.05,0.2,0.015\n"
            "unused,0,0,0,0,0,0,0,0,0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            path.write_text(contents, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "strictly positive"):
                read_rows(path)

        zero_weight_rows = [
            {
                "stratum": "unused",
                "reference_weight": 0.0,
                "rho_pre": 0.0,
                "se_rho_pre": 0.0,
                "pi_pre": 0.0,
                "se_pi_pre": 0.0,
                "rho_post": 0.0,
                "se_rho_post": 0.0,
                "pi_post": 0.0,
                "se_pi_post": 0.0,
            }
        ]
        with self.assertRaisesRegex(ValueError, "strictly positive"):
            project(zero_weight_rows, alpha=0.05)


if __name__ == "__main__":
    unittest.main()
