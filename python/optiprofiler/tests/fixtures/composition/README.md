# Legacy single-feature goldens

`legacy_single_feature_goldens.json` pins the single-feature behaviour of the
`python` baseline `4f821d8658700b8f3ec0bd480962310a74e47667` (before Feature
composition): 176 scenarios (problem kinds u/b/l/n, n = 3, seeds 0 and 12345,
`MAX_EVAL = 4`) recording structure, oracle outputs, reference histories,
budget behaviour and the exact root-callback call sequence of every feature.

Provenance recorded inside the fixture:

- `source_sha`: the base commit the scenarios were run against;
- `scenario_module_sha256`
  `2ba635ad3225f9205d4854d8f9f007ffeed8f9b4d74e73e0dc1e4a373be0ef1d`: the
  generator `composition_goldens.py` as of commit `269358a`, which produced
  the fixture. The shipped generator was rewritten afterwards (it imports
  `optiprofiler.experiment.resolve_plan` and reads `Feature.stages`, which the
  base commit does not have), so it cannot regenerate the fixture at the base
  commit. `test_composition_goldens.py` asserts the recorded hash so that an
  edited fixture is detected; a regenerated fixture must record the hash of
  the generator that wrote it and the commit it ran against.

The fixture is a regression net, not the only evidence: independent
base-versus-candidate differentials (522 end-to-end cases and 3,510
low-level scenarios in the 2026-09 audits) confirmed the pinned behaviour.
