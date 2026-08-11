# Synchronizing SED functionality from MDtrace

pySED and MDtrace are maintained as separate projects. Their current Git
histories do not share a merge base, so SED updates are ported by behavior;
MDtrace commits are not merged or cherry-picked into pySED.

## Recorded source snapshot

The pySED 3.0 synchronization uses the following read-only source:

- Repository: `https://github.com/Tingliangstu/MDtrace.git`
- Branch: `main`
- Commit: `a61eae3a1a3c9b9ba75bafee52f4ca4b220d6178`
- Commit date: `2026-08-04T21:16:15+08:00`
- Capture date: `2026-08-11T10:22:28+08:00`

That MDtrace worktree also contained uncommitted SED changes for
`dump_group_supercell` and scale-independent fitting of small-amplitude SED
peaks. The exact working-file hashes are stored in
`sync/mdtrace-sed-snapshot.json` and
`sync/mdtrace-sed-dirty-files.sha256`. A hash entry describes the complete
working file, but only the SED behavior identified by its `scope` field is in
the pySED synchronization scope.

In particular, `mdtrace/parameters.py` and `mdtrace/pipeline.py` contain both
SED and DSF working-tree changes. Do not copy those files wholesale. Port only
the `dump_group_supercell` parameter, its SED validation, and selection of the
local cell for SED Brillouin-zone construction.

The manifest records integrity and provenance; it does not vendor MDtrace
source files into pySED. Preserve copyright and attribution from both projects
when porting implementation content, and adapt namespaces, messages, cache
suffixes, and documentation to pySED.

## Deliberate exclusions

This synchronization is SED-only. It excludes:

- all DSF and EELS modules, parameters, documentation, examples, tests, and
  roadmap items;
- the pySED `add_dynamic_structure_factor` branch;
- MDtrace's uncommitted mode-projected SED roadmap in `TODO.md`;
- generated `tmp/` content and the untracked
  `ni_diamond_111_10nm_24x24_sed/` analysis tree.

Names such as `method = dsf`, `method = eels`, `atom_types`, `dsf_qpoints`,
and `experiment` must therefore not become pySED 3.0 public interfaces.

## Future synchronization procedure

1. Record the new MDtrace branch and full commit hash before inspecting local
   modifications. Keep the MDtrace worktree read-only.
2. Compare SED behavior and tests against the previous snapshot. Ignore the
   excluded subsystems above, including changes mixed into shared files.
3. Port the selected behavior into pySED namespaces and add pySED regression
   tests; do not merge or cherry-pick unrelated history.
4. Update the JSON and SHA256 manifests together. Hash files as raw bytes with
   SHA-256 and record whether each source file was tracked or untracked.
5. Run pySED's CPU test matrix, documentation build, package build, and CLI
   smoke checks. Confirm the MDtrace `git status` and source hashes are
   unchanged after the port.

The SHA256 list can be checked from an MDtrace checkout with either
`sha256sum --check <path-to-pySED>/sync/mdtrace-sed-dirty-files.sha256` or an
equivalent platform-native SHA-256 tool. The paths in that file are relative
to the MDtrace repository root.
