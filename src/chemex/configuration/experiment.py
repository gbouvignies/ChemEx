from __future__ import annotations

import math
from functools import cached_property
from typing import Literal, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from chemex.configuration.utils import key_to_lower

T = TypeVar("T")


class ExperimentSettings(BaseModel):
    observed_state: Literal["a", "b", "c", "d"] = "a"

    model_config = ConfigDict(str_to_lower=True)

    _key_to_lower = model_validator(mode="before")(key_to_lower)


class RelaxationSettings(ExperimentSettings):
    """Base settings for relaxation-based experiments (CEST, DCEST, CPMG, etc).

    Attributes:
        cs_evolution_prior (bool):
            Controls the initial condition (equilibrium vs non-equilibrium).
            If False (default), the initial magnetization assumes thermal
            equilibrium with all states populated according to their equilibrium
            populations. If True, only the observed state is initially populated
            (non-equilibrium condition), which is appropriate when chemical shift
            evolution occurs before the CEST/CPMG element in the pulse sequence,
            breaking the equilibrium assumption. This can affect the accuracy of
            extracted kinetic parameters, especially for slow exchange rates
            (see Yuwen et al., J. Biomol. NMR 2016, 65:143-156).
        detect_all_states (bool):
            If True, detection will use all states (e.g., '[iz]'),
            otherwise only the observed state (e.g., '[iz_a]').

    Properties:
        suffix_start: Suffix for starting terms, respects cs_evolution_prior.
        suffix_detect: Suffix for detection, respects detect_all_states.
            Returns '' if detect_all_states is True, else '_{observed_state}'.

    """

    cs_evolution_prior: bool = False
    detect_all_states: bool = False

    @cached_property
    def suffix_start(self) -> str:
        """Suffix for starting terms.

        Returns '_{observed_state}' if cs_evolution_prior is True (non-equilibrium
        initial condition with only the observed state populated), else ''
        (equilibrium initial condition with all states populated).
        """
        return f"_{self.observed_state}" if self.cs_evolution_prior else ""

    @cached_property
    def suffix_detect(self) -> str:
        """Suffix for detection.

        Returns '' if detect_all_states is True, else '_{observed_state}'.
        """
        return "" if self.detect_all_states else f"_{self.observed_state}"


class CpmgSettings(RelaxationSettings):
    even_ncycs: bool = False


class CpmgSettingsEvenNcycs(RelaxationSettings):
    even_ncycs: bool = True


class CestSettings(RelaxationSettings):
    sw: float = math.inf


class MFCestSettings(RelaxationSettings):
    sw: float
