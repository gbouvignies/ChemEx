from chemex import helper as ch

_SCHEMA_METHODS_PARAM = {
    "type": "object",
    "additionalProperties": {
        "type": "object",
        "properties": {
            "include": {
                "oneOf": [
                    {
                        "type": "array",
                        "items": {"anyOf": [{"type": "integer"}, {"type": "string"}]},
                    },
                    {"type": "string"},
                ]
            },
            "exclude": {
                "oneOf": [
                    {
                        "type": "array",
                        "items": {"anyOf": [{"type": "integer"}, {"type": "string"}]},
                    },
                    {"type": "string"},
                ]
            },
            "fit": {"type": "array", "items": {"type": "string"}},
            "fix": {"type": "array", "items": {"type": "string"}},
            "constraints": {"type": "array", "items": {"type": "string"}},
            "grid": {"type": "array", "items": {"type": "string"}},
            "statistics": {
                "type": "object",
                "properties": {
                    "mc": {"type": "integer", "minimum": 0},
                    "bs": {"type": "integer", "minimum": 0},
                    "bsn": {"type": "integer", "minimum": 0},
                },
            },
        },
    },
}


def read_methods(filenames):
    print("\nReading methods...")
    return ch.read_toml_multi(filenames, _SCHEMA_METHODS_PARAM)
