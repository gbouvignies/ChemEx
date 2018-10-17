import inspect
import pkgutil
from importlib import import_module

from chemex.experiments.base.base_profile import BaseProfile


def grab(exp_name):

    try:
        imported_module = import_module(f"{__name__}.{exp_name}")

    except ModuleNotFoundError:
        raise ImportError(
            "{} is not part of our experiment collection!".format(exp_name)
        )

    for dir_name in dir(imported_module):

        attribute = getattr(imported_module, dir_name)

        if (
            inspect.isclass(attribute)
            and issubclass(attribute, BaseProfile)
            and not inspect.isabstract(attribute)
        ):
            profile_class = attribute
            break

    else:
        raise ImportError(
            f"We currently don't have {exp_name}, but you are welcome to send in the "
            f"request for it!"
        )

    return profile_class


def get_experiment_docs():

    docs = {}

    for (_, name, _) in pkgutil.walk_packages(__path__, __name__ + "."):

        imported_module = import_module(f"{name}")

        for dir_name in vars(imported_module):
            attribute = getattr(imported_module, dir_name)

            if (
                inspect.isclass(attribute)
                and issubclass(attribute, BaseProfile)
                and not inspect.isabstract(attribute)
            ):
                exp_name = name.replace(__name__ + ".", "")
                docs[exp_name] = imported_module.__doc__

    return docs
