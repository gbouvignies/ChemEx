"""Tests for B1 distribution configuration schemas, validation, and migration."""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest
from pydantic import ValidationError

from chemex.configuration.b1_config import (
    CustomDistributionConfig,
    DephasingConfig,
    GaussianDistributionConfig,
    HomogeneousDistributionConfig,
    SkewedDistributionConfig,
    parse_distribution_config,
)
from chemex.configuration.experiment import B1InhomogeneityMixin
from chemex.configuration.methods import Selection
from chemex.containers.experiment import DeprecatedB1InhomogeneityNotice
from chemex.experiments import experiment_types
from chemex.experiments.catalog import (
    cest_1hn_ap,
    cest_1hn_ip_ap,
    cest_13c,
    cest_15n,
    cest_15n_cw,
    cest_15n_tr,
    cest_ch3_1h_ip_ap,
    coscest_1hn_ip_ap,
    coscest_13c,
    dcest_13c,
    dcest_15n,
)
from chemex.experiments.catalog.wip import (
    cest_15n_test,
    relaxation_15n_r1rho,
    relaxation_15n_r1rho_eig,
)
from chemex.experiments.experiment_types import ExperimentConfigurationError
from chemex.experiments.loader import register_experiments
from chemex.models.model import ModelSpec
from chemex.nmr.distributions._quadrature import MAX_SKEWNESS
from chemex.runtime import AnalysisSession, ensure_plugins_registered

ROOT = Path(__file__).parents[2]
EXAMPLES_DIR = ROOT / "examples"

_CW_CATALOG_MODULES = [
    (
        cest_13c,
        {"name": "cest_13c", "time_t1": 0.3, "carrier": 18.0, "b1_frq": 15.0},
    ),
    (
        cest_15n,
        {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 25.0},
    ),
    (
        cest_15n_cw,
        {
            "name": "cest_15n_cw",
            "time_t1": 0.5,
            "carrier": 118.0,
            "carrier_dec": 8.0,
            "b1_frq": 25.0,
            "b1_frq_dec": 200.0,
        },
    ),
    (
        cest_15n_tr,
        {"name": "cest_15n_tr", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 25.0},
    ),
    (
        cest_1hn_ap,
        {"name": "cest_1hn_ap", "time_t1": 0.3, "carrier": 8.0, "b1_frq": 25.0},
    ),
    (
        cest_1hn_ip_ap,
        {
            "name": "cest_1hn_ip_ap",
            "time_t1": 0.3,
            "carrier": 8.0,
            "b1_frq": 25.0,
            "d1": 0.5,
        },
    ),
    (
        cest_ch3_1h_ip_ap,
        {
            "name": "cest_ch3_1h_ip_ap",
            "time_t1": 0.3,
            "carrier": 0.8,
            "b1_frq": 25.0,
            "d1": 0.5,
        },
    ),
]

_NON_CW_CATALOG_MODULES = [
    (
        coscest_13c,
        {
            "name": "coscest_13c",
            "time_t1": 0.5,
            "carrier": 55.0,
            "b1_frq": 25.0,
            "sw": 800.0,
            "cos_n": 3,
        },
    ),
    (
        coscest_1hn_ip_ap,
        {
            "name": "coscest_1hn_ip_ap",
            "time_t1": 0.5,
            "carrier": 8.0,
            "b1_frq": 25.0,
            "d1": 0.5,
            "sw": 800.0,
            "cos_n": 3,
        },
    ),
    (
        dcest_13c,
        {
            "name": "dcest_13c",
            "time_t1": 0.2,
            "carrier": 55.0,
            "b1_eff": 20.0,
            "pw90": 40.0e-6,
            "sw": 700.0,
        },
    ),
    (
        dcest_15n,
        {
            "name": "dcest_15n",
            "time_t1": 0.2,
            "carrier": 118.0,
            "b1_eff": 20.0,
            "pw90": 40.0e-6,
            "sw": 700.0,
        },
    ),
]

ALL_CEST_CATALOG_MODULES = _CW_CATALOG_MODULES + _NON_CW_CATALOG_MODULES


class DummyB1Settings(B1InhomogeneityMixin):
    pass


class DummyDefaultGaussianB1Settings(B1InhomogeneityMixin):
    legacy_b1_inh_scale_default: ClassVar[float | None] = 0.2
    legacy_b1_inh_res_default: ClassVar[int | None] = 7


class TestStaticDistributionUnion:
    """Verify the five-model discriminated union and error diagnostics."""

    def test_parses_all_five_distribution_types_and_case_insensitivity(self) -> None:
        h = parse_distribution_config({"TYPE": "HOMOGENEOUS"})
        assert isinstance(h, HomogeneousDistributionConfig)

        d = parse_distribution_config({"type": "Dephasing"})
        assert isinstance(d, DephasingConfig)

        g = parse_distribution_config({"type": "Gaussian", "scale": 0.12})
        assert isinstance(g, GaussianDistributionConfig)
        assert g.scale == 0.12
        assert g.res is None

        s = parse_distribution_config(
            {
                "TYPE": "SKEWED",
                "scale": 0.1,
                "skewness": -0.5,
                "CENTER": "MODE",
                "res": 21,
            }
        )
        assert isinstance(s, SkewedDistributionConfig)
        assert s.center == "mode"
        assert s.res == 21

        c = parse_distribution_config(
            {"type": "custom", "scales": [0.9, 1.0], "weights": [0.4, 0.6]}
        )
        assert isinstance(c, CustomDistributionConfig)
        assert c.scales == [0.9, 1.0]

    def test_skewed_requires_explicit_center_with_actionable_error(self) -> None:
        with pytest.raises(ValidationError, match="requires an explicit 'center'"):
            parse_distribution_config(
                {"type": "skewed", "scale": 0.1, "skewness": -0.5}
            )

        with pytest.raises(ValidationError, match=r'center = "mode"'):
            parse_distribution_config(
                {"type": "skewed", "scale": 0.1, "skewness": -0.5, "center": None}
            )

    @pytest.mark.parametrize(
        "bad_skewness",
        [-MAX_SKEWNESS, MAX_SKEWNESS, -1.0, 1.0, 2.0],
    )
    def test_skewed_rejects_out_of_range_skewness(self, bad_skewness: float) -> None:
        with pytest.raises(ValidationError, match="skew-normal range"):
            parse_distribution_config(
                {
                    "type": "skewed",
                    "scale": 0.1,
                    "skewness": bad_skewness,
                    "center": "mode",
                }
            )

    def test_untyped_dict_and_retired_types_emit_migration_errors(self) -> None:
        with pytest.raises((ValidationError, ValueError), match="must specify 'type'"):
            parse_distribution_config({})

        with pytest.raises(ValidationError, match="must specify 'type'"):
            cest_15n.Cest15NSettings.model_validate(
                {
                    "name": "cest_15n",
                    "time_t1": 0.5,
                    "carrier": 118.0,
                    "b1_frq": 25.0,
                    "b1_distribution": {},
                }
            )

        with pytest.raises(
            (ValidationError, ValueError), match="Unknown B1 distribution type"
        ):
            parse_distribution_config({"type": "nonexistent"})

        with pytest.raises(
            (ValidationError, ValueError), match="consolidated into 'gaussian'"
        ):
            parse_distribution_config({"type": "hermite", "scale": 0.1})

        for retired in ("beta", "truncated_skewed"):
            with pytest.raises((ValidationError, ValueError), match="unphysical"):
                parse_distribution_config({"type": retired, "scale": 0.1})
            with pytest.raises(ValidationError, match="unphysical"):
                cest_15n.Cest15NSettings.model_validate(
                    {
                        "name": "cest_15n",
                        "time_t1": 0.5,
                        "carrier": 118.0,
                        "b1_frq": 25.0,
                        "b1_distribution": {"type": retired, "scale": 0.1},
                    }
                )

    def test_extra_keys_forbidden(self) -> None:
        with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
            parse_distribution_config({"type": "gaussian", "sclae": 0.1})


class TestCatalogCestB1Contracts:
    """Validate mandatory b1_distribution, dephasing applicability, CW XOR, and D-CEST guards."""

    @staticmethod
    def _build_from_dict(
        tmp_path: Path,
        exp_dict: dict[str, object],
        *,
        spin_name: str = "10N-HN",
    ) -> experiment_types.ExperimentBuildResult:
        register_experiments()
        data_file = tmp_path / "profile.out"
        data_file.write_text("-20000 1.0 0.02\n0 0.5 0.02\n", encoding="utf-8")

        def format_toml_val(val: object) -> str:
            if isinstance(val, str):
                return f'"{val}"'
            if isinstance(val, dict):
                inner = ", ".join(f"{k} = {format_toml_val(v)}" for k, v in val.items())
                return f"{{ {inner} }}"
            if isinstance(val, list):
                inner = ", ".join(format_toml_val(v) for v in val)
                return f"[{inner}]"
            return repr(val)

        exp_lines = "\n".join(
            f"{k} = {format_toml_val(v)}" for k, v in exp_dict.items()
        )
        exp_file = tmp_path / f"{exp_dict['name']}.toml"
        exp_file.write_text(
            f"""[experiment]
{exp_lines}

[conditions]
h_larmor_frq = 600.0

[data]
path = "."

[data.profiles]
{spin_name} = "profile.out"
""",
            encoding="utf-8",
        )
        session = AnalysisSession.create()
        session.set_model("2st")
        source = experiment_types.open(exp_file)
        return experiment_types.build(
            source,
            selection=Selection(include=None, exclude=None),
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

    @pytest.mark.parametrize(("module", "exp_base"), ALL_CEST_CATALOG_MODULES)
    def test_all_11_cest_modules_reject_omitted_b1_distribution(
        self,
        tmp_path: Path,
        module,
        exp_base: dict[str, object],
    ) -> None:
        with pytest.raises(ExperimentConfigurationError) as err_info:
            self._build_from_dict(tmp_path, exp_base)
        message = str(err_info.value)
        assert "require an explicit 'b1_distribution'" in message
        allows_dephasing = getattr(module, "__name__", "")
        if any(
            non_cw in allows_dephasing
            for non_cw in ("coscest_13c", "coscest_1hn_ip_ap", "dcest_13c", "dcest_15n")
        ):
            assert 'type = "dephasing"' not in message
        else:
            assert 'type = "dephasing"' in message

    @pytest.mark.parametrize(("module", "exp_base"), _CW_CATALOG_MODULES)
    def test_cw_cest_modules_accept_dephasing(
        self,
        module,
        exp_base: dict[str, object],
    ) -> None:
        settings_cls = module.EXPERIMENT_TYPE.config_type.model_fields[
            "experiment"
        ].annotation
        settings = settings_cls.model_validate(
            {**exp_base, "b1_distribution": {"type": "dephasing"}}
        )
        assert isinstance(settings.b1_distribution, DephasingConfig)

    @pytest.mark.parametrize(("module", "exp_base"), _NON_CW_CATALOG_MODULES)
    def test_non_cw_modules_reject_dephasing(
        self,
        tmp_path: Path,
        module,
        exp_base: dict[str, object],
    ) -> None:
        with pytest.raises(
            ExperimentConfigurationError, match="cannot be approximated"
        ):
            self._build_from_dict(
                tmp_path,
                {**exp_base, "b1_distribution": {"type": "dephasing"}},
            )

    def test_cest_1hn_ip_ap_dephasing_allowed_only_when_eta_block_is_zero(
        self,
        tmp_path: Path,
    ) -> None:
        base = {
            "name": "cest_1hn_ip_ap",
            "time_t1": 0.3,
            "carrier": 8.0,
            "b1_frq": 25.0,
            "d1": 0.5,
            "b1_distribution": {"type": "dephasing"},
        }
        allowed = cest_1hn_ip_ap.Cest1HnIpApSettings.model_validate(
            {**base, "eta_block": 0}
        )
        assert isinstance(allowed.b1_distribution, DephasingConfig)

        with pytest.raises(ExperimentConfigurationError, match="eta_block"):
            self._build_from_dict(tmp_path, {**base, "eta_block": 4})

    def test_cw_and_coscest_enforce_b1_frq_xor_pw90_and_nominal_from_pw90(
        self,
        tmp_path: Path,
    ) -> None:
        # Both b1_frq and pw90 omitted -> ExperimentConfigurationError
        with pytest.raises(
            ExperimentConfigurationError,
            match="Either 'b1_frq' .* or 'pw90'",
        ):
            self._build_from_dict(
                tmp_path,
                {
                    "name": "cest_15n",
                    "time_t1": 0.5,
                    "carrier": 118.0,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                },
            )

        # Both b1_frq and pw90 provided -> ExperimentConfigurationError
        with pytest.raises(
            ExperimentConfigurationError,
            match="Use either 'b1_frq' or 'pw90'",
        ):
            self._build_from_dict(
                tmp_path,
                {
                    "name": "coscest_13c",
                    "time_t1": 0.5,
                    "carrier": 55.0,
                    "b1_frq": 25.0,
                    "pw90": 0.01,
                    "sw": 800.0,
                    "cos_n": 3,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                },
            )

        # pw90 only works cleanly and yields nominal = 1 / (4 * pw90)
        settings = cest_15n.Cest15NSettings.model_validate(
            {
                "name": "cest_15n",
                "time_t1": 0.5,
                "carrier": 118.0,
                "pw90": 0.01,
                "b1_distribution": {"type": "homogeneous"},
            }
        )
        assert settings.get_b1_nominal() == pytest.approx(25.0)

    def test_cli_reports_cw_b1_frq_pw90_conflict_as_known_error(
        self,
        tmp_path: Path,
    ) -> None:
        data_file = tmp_path / "10N-HN.out"
        data_file.write_text("-20000 1.0 0.02\n0 0.5 0.02\n", encoding="utf-8")
        exp_file = tmp_path / "bad_xor.toml"
        exp_file.write_text(
            """[experiment]
name = "cest_15n"
time_t1 = 0.5
carrier = 118.0
b1_frq = 25.0
pw90 = 0.01
b1_distribution = { type = "homogeneous" }

[conditions]
h_larmor_frq = 600.0

[data]
path = "."

[data.profiles]
10N-HN = "10N-HN.out"
""",
            encoding="utf-8",
        )
        param_file = tmp_path / "params.toml"
        param_file.write_text("[GLOBAL]\nPB = 0.1\nKEX_AB = 100.0\n", encoding="utf-8")

        completed = subprocess.run(  # noqa: S603 - fixed local python entrypoint check
            [
                sys.executable,
                "-m",
                "chemex",
                "simulate",
                "-e",
                str(exp_file),
                "-p",
                str(param_file),
                "-o",
                str(tmp_path / "out"),
                "--plot",
                "nothing",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        combined = completed.stdout + completed.stderr
        assert completed.returncode == 1
        assert "Experiment configuration is invalid" in completed.stderr
        assert "Use either 'b1_frq' or 'pw90'" in completed.stderr
        assert "unexpected internal error" not in combined
        assert "Traceback" not in combined

    def test_dcest_rejects_non_positive_sw_and_duty_cycle_overflow(self) -> None:
        with pytest.raises(ValidationError, match="greater than 0"):
            dcest_15n.DCest15NSettings.model_validate(
                {
                    "name": "dcest_15n",
                    "time_t1": 0.2,
                    "carrier": 118.0,
                    "b1_eff": 20.0,
                    "pw90": 40.0e-6,
                    "sw": 0.0,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                }
            )

        with pytest.raises(ValidationError, match="duty cycle 4 \\* pw90 \\* b1_eff"):
            dcest_15n.DCest15NSettings.model_validate(
                {
                    "name": "dcest_15n",
                    "time_t1": 0.2,
                    "carrier": 118.0,
                    "b1_eff": 300.0,
                    "pw90": 1.0e-3,
                    "sw": 500.0,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                }
            )


class TestLegacyB1CompatibilityAndNotices:
    """Verify one-release deprecation mapping, single notice per file, and WIP behavior."""

    def test_legacy_mappings_and_res_preservation(self) -> None:
        g_auto = cest_15n.Cest15NSettings.model_validate(
            {
                "name": "cest_15n",
                "time_t1": 0.5,
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": 0.1,
            }
        )
        assert isinstance(g_auto.b1_distribution, GaussianDistributionConfig)
        assert g_auto.b1_distribution.scale == 0.1
        assert g_auto.b1_distribution.res is None
        assert g_auto.legacy_b1_deprecation == (
            ("b1_inh_scale",),
            '{ type = "gaussian", scale = 0.1 }',
        )

        g_res_only = dcest_15n.DCest15NSettings.model_validate(
            {
                "name": "dcest_15n",
                "time_t1": 0.2,
                "carrier": 118.0,
                "b1_eff": 20.0,
                "pw90": 40.0e-6,
                "sw": 700.0,
                "b1_inh_res": 9,
            }
        )
        assert isinstance(g_res_only.b1_distribution, GaussianDistributionConfig)
        assert g_res_only.b1_distribution.scale == 0.1
        assert g_res_only.b1_distribution.res == 9
        assert g_res_only.legacy_b1_deprecation == (
            ("b1_inh_res",),
            '{ type = "gaussian", scale = 0.1, res = 9 }',
        )

        h_zero = DummyB1Settings.model_validate({"b1_frq": 15.0, "b1_inh_scale": 0.0})
        assert isinstance(h_zero.b1_distribution, HomogeneousDistributionConfig)

        h_res1 = DummyB1Settings.model_validate({"b1_frq": 15.0, "b1_inh_res": 1})
        assert isinstance(h_res1.b1_distribution, HomogeneousDistributionConfig)

        d_cw = cest_15n.Cest15NSettings.model_validate(
            {
                "name": "cest_15n",
                "time_t1": 0.5,
                "carrier": 118.0,
                "b1_frq": 25.0,
                "b1_inh_scale": np.inf,
            }
        )
        assert isinstance(d_cw.b1_distribution, DephasingConfig)

        with pytest.raises(ValidationError, match="cannot be approximated"):
            dcest_15n.DCest15NSettings.model_validate(
                {
                    "name": "dcest_15n",
                    "time_t1": 0.2,
                    "carrier": 118.0,
                    "b1_eff": 20.0,
                    "pw90": 40.0e-6,
                    "sw": 700.0,
                    "b1_inh_scale": np.inf,
                }
            )

    def test_rejects_mixing_legacy_keys_and_b1_distribution_including_instances(
        self,
    ) -> None:
        with pytest.raises(
            ValidationError, match="Use either 'b1_distribution' or the legacy"
        ):
            DummyB1Settings.model_validate(
                {
                    "b1_frq": 15.0,
                    "b1_distribution": {"type": "gaussian", "scale": 0.1},
                    "b1_inh_scale": 0.2,
                }
            )

        with pytest.raises(
            ValidationError, match="Use either 'b1_distribution' or the legacy"
        ):
            DummyB1Settings.model_validate(
                {
                    "b1_frq": 15.0,
                    "b1_distribution": GaussianDistributionConfig(scale=0.1, res=5),
                    "b1_inh_res": 7,
                }
            )

    def test_build_emits_exactly_one_deprecation_notice_per_multi_profile_file(
        self,
        tmp_path: Path,
    ) -> None:
        register_experiments()
        for name in ("10N-HN.out", "11N-HN.out", "12N-HN.out"):
            (tmp_path / name).write_text(
                "-20000 1.0 0.02\n0 0.5 0.02\n", encoding="utf-8"
            )
        exp_file = tmp_path / "legacy.toml"
        exp_file.write_text(
            """[experiment]
name = "cest_15n"
time_t1 = 0.5
carrier = 118.0
b1_frq = 25.0
b1_inh_scale = 0.1
b1_inh_res = 11

[conditions]
h_larmor_frq = 600.0

[data]
path = "."

[data.profiles]
10N-HN = "10N-HN.out"
11N-HN = "11N-HN.out"
12N-HN = "12N-HN.out"
""",
            encoding="utf-8",
        )
        session = AnalysisSession.create()
        session.set_model("2st")
        source = experiment_types.open(exp_file)
        result = experiment_types.build(
            source,
            selection=Selection(include=None, exclude=None),
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

        # Verify model_name assignment with validate_assignment=True succeeds cleanly.
        assert len(result.experiment) == 3
        b1_notices = [
            n for n in result.notices if isinstance(n, DeprecatedB1InhomogeneityNotice)
        ]
        assert len(b1_notices) == 1
        assert b1_notices[0].filename == exp_file
        assert b1_notices[0].legacy_fields == ("b1_inh_scale", "b1_inh_res")
        assert (
            b1_notices[0].replacement == '{ type = "gaussian", scale = 0.1, res = 11 }'
        )

    def test_wip_default_legacy_classes_emit_zero_deprecation_notices(self) -> None:
        s1 = cest_15n_test.Cest15NSettings.model_validate(
            {
                "name": "wip.cest_15n_test",
                "time_t1": 0.2,
                "carrier": 118.0,
                "b1_frq": 25.0,
            }
        )
        assert isinstance(s1.b1_distribution, GaussianDistributionConfig)
        assert s1.legacy_b1_deprecation is None

        s2 = relaxation_15n_r1rho.Relaxation15NR1RhoSettings.model_validate(
            {
                "name": "wip.relaxation_15n_r1rho",
                "carrier": 118.0,
                "b1_frq": 25.0,
            }
        )
        assert isinstance(s2.b1_distribution, DephasingConfig)
        assert s2.legacy_b1_deprecation is None

        s3 = relaxation_15n_r1rho_eig.Relaxation15NR1RhoSettings.model_validate(
            {
                "name": "wip.relaxation_15n_r1rho_eig",
                "carrier": 118.0,
                "b1_frq": 25.0,
            }
        )
        assert isinstance(s3.b1_distribution, DephasingConfig)
        assert s3.legacy_b1_deprecation is None

    def test_dcest_b1_frq_alias_preserved(self) -> None:
        settings = dcest_15n.DCest15NSettings.model_validate(
            {
                "name": "dcest_15n",
                "time_t1": 0.2,
                "carrier": 118.0,
                "b1_frq": 20.0,
                "pw90": 40.0e-6,
                "sw": 700.0,
                "b1_distribution": {"type": "gaussian", "scale": 0.1},
            }
        )
        assert settings.b1_eff == 20.0
        assert settings.b1_frq is None
        settings.model_name = "2st"
        assert settings.b1_eff == 20.0

    def test_all_shipped_cest_example_tomls_build_without_deprecation_notices(
        self,
    ) -> None:
        ensure_plugins_registered()
        exp_types = {
            mod.EXPERIMENT_NAME: mod.EXPERIMENT_TYPE
            for mod, _ in ALL_CEST_CATALOG_MODULES
        }
        cest_prefixes = ("cest_", "coscest_", "dcest_")
        shipped_tomls: list[tuple[Path, dict[str, object]]] = []
        for toml_path in sorted(EXAMPLES_DIR.rglob("*.toml")):
            rel_parts = toml_path.relative_to(EXAMPLES_DIR).parts
            if any(
                part.startswith("Output") or part == "run_info" for part in rel_parts
            ):
                continue
            payload = tomllib.loads(toml_path.read_text(encoding="utf-8"))
            exp_table = payload.get("experiment")
            if isinstance(exp_table, dict) and str(
                exp_table.get("name", "")
            ).startswith(cest_prefixes):
                shipped_tomls.append((toml_path, payload))

        assert len(shipped_tomls) == 39
        for toml_path, payload in shipped_tomls:
            source = experiment_types.open(toml_path)
            exp_type = exp_types[source.experiment_type_name]
            model_name = "3st" if "3States" in str(toml_path) else "2st"
            config = exp_type.config_type.model_validate(
                payload, context={"model": ModelSpec.from_name(model_name)}
            )
            assert config.experiment.legacy_b1_deprecation is None, str(toml_path)
            assert config.experiment.b1_distribution is not None, str(toml_path)
