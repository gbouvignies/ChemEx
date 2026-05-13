from chemex.parameters.spin_system import SpinSystem


def test_spin_system_name_parser_does_not_mutate_input_dict() -> None:
    raw = {"name": "G23N-HN"}

    spin_system = SpinSystem.model_validate(raw)

    assert spin_system.name == "G23N-H"
    assert raw == {"name": "G23N-HN"}
