"""Shared access to the quality-filter implementation."""
import importlib.util
from pathlib import Path
_spec=importlib.util.spec_from_file_location("quality_filter",Path(__file__).with_name("04_quality_cadd.py"))
_mod=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
load_controls=_mod.load_controls
to_num=_mod.to_num
