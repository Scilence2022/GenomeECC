"""Make the flat top-level modules importable from tests/ without installing."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
