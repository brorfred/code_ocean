"""Configuration management for seagrid using Dynaconf.

Loads settings from ``settings.toml`` and ``.secrets.toml`` files with
support for multiple environments. The default environment is ``no_regions``.
"""

import pathlib

from dynaconf import Dynaconf

preload_file = pathlib.Path(__file__).parent/"settings.toml"
settings = Dynaconf(
    environments=True,
    envvar_prefix="DYNACONF",
    env="no_regions",
    preload=[preload_file,],
    settings_files=['settings.toml', '.secrets.toml'],
)

# `envvar_prefix` = export envvars with `export DYNACONF_FOO=bar`.
# `settings_files` = Load these files in the order.
