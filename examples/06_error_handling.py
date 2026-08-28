"""What py4dggs raises, and when.

Run:  python examples/06_error_handling.py
"""
from py4dggs import IGEO7, ISEA3H, InvalidZoneError

def show(label, fn):
    try:
        result = fn()
    except Exception as e:
        print(f"  {label:38} {type(e).__name__}: {e}")
    else:
        print(f"  {label:38} -> {result}")

print("resolution bounds (IGEO7 max is 19, ISEA3H max is 33):")
show("IGEO7 res 19", lambda: IGEO7.zone_from_geo(0.0, 0.0, 19).text_id)
show("IGEO7 res 20", lambda: IGEO7.zone_from_geo(0.0, 0.0, 20).text_id)
show("IGEO7 res -1", lambda: IGEO7.zone_from_geo(0.0, 0.0, -1).text_id)

print("\nunparseable text ids:")
show("IGEO7 'not-a-zone'", lambda: IGEO7.zone_from_text("not-a-zone"))
show("IGEO7 '0999'", lambda: IGEO7.zone_from_text("0999"))
show("ISEA3H 'nonsense'", lambda: ISEA3H.zone_from_text("nonsense"))

print("\nend of the hierarchy is empty, not an error:")
base = IGEO7.zone_from_geo(0.0, 0.0, 0)
show("parents of a res-0 zone", lambda: base.parents)
show("children at max resolution", lambda: IGEO7.zone_from_geo(0.0, 0.0, 19).children)

print("\nInvalidZoneError is the one to catch:")
try:
    IGEO7.zone_from_geo(0.0, 0.0, 99)
except InvalidZoneError as e:
    print(f"  caught InvalidZoneError: {e}")
