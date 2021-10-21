FILTER_PLANES_SCHEMA = {
    "type": "array",
    "items": {"type": "integer"},
    "default": [],
}

PATH_SCHEMA = {"type": "string", "default": "./"}

PROFILES_LIST_SCHEMA = {
    "type": "array",
    "items": {
        "type": "array",
        "minItems": 2,
        "maxItems": 2,
        "items": {"type": "string"},
    },
}

PROFILES_DICT_SCHEMA = {
    "type": "object",
    "additionalProperties": {
        "oneOf": [
            {"type": "string"},
            {"type": "array", "items": {"type": "string"}},
        ]
    },
}

CEST_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "scatter"],
                    "default": "file",
                },
                "filter_offsets": {
                    "type": "array",
                    "items": {
                        "type": "array",
                        "minItems": 2,
                        "maxItems": 2,
                        "items": {"type": "number"},
                    },
                    "default": [[0.0, 0.0]],
                },
                "filter_planes": FILTER_PLANES_SCHEMA,
                "filter_ref_planes": {"type": "boolean", "default": False},
                "path": PATH_SCHEMA,
                "profiles": {"oneOf": [PROFILES_LIST_SCHEMA, PROFILES_DICT_SCHEMA]},
            },
            "required": ["profiles"],
        }
    },
}

CPMG_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "duplicates"],
                    "default": "file",
                },
                "filter_planes": FILTER_PLANES_SCHEMA,
                "path": PATH_SCHEMA,
                "profiles": {"oneOf": [PROFILES_LIST_SCHEMA, PROFILES_DICT_SCHEMA]},
            },
            "required": ["profiles"],
        }
    },
}

RELAXATION_SCHEMA = CPMG_SCHEMA

SHIFT_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "duplicates"],
                    "default": "file",
                },
                "path": PATH_SCHEMA,
            },
            "required": ["path"],
        }
    },
}
