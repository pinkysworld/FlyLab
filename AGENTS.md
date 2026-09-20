# FlyLab agent instructions

These rules apply to any automated coding agent, LLM assistant, bot, or human using agent-generated patches in this repository.

## Commit attribution policy

AI-assisted work is allowed.

However, **Git commits must be attributed to the repository owner or the human committer, not to the AI system that helped create the change.**

Do not create commits with:

- author name `Claude`
- committer name `Claude`
- any author/committer email at `@anthropic.com`
- `Co-Authored-By: Claude ...` trailers
- `Claude-Session: ...` trailers

The same principle applies to any future automated agent identity that would incorrectly appear as a human project contributor unless the repository owner explicitly requests otherwise.

Use the repository owner's configured Git identity for commits created on the owner's behalf.

This policy changes attribution only. It does **not** require removing, squashing, or concealing AI-assisted code, prose, analyses, or other substantive contributions.

CI enforces the current Claude/Anthropic attribution rule with:

```bash
python scripts/check_commit_attribution.py
```

If the check fails, rewrite the offending commit metadata while preserving the commit's file tree/content.

## Scientific integrity

Before changing scientific claims or generated paper artifacts, read:

- `HANDOFF.md`
- `docs/ARCHITECTURE.md`
- `docs/NOVELTY.md`
- `papers/STATUS.md`

Do not hand-edit generated numerical paper results. Use the reproduction pipeline.

## Minimum checks

For ordinary changes:

```bash
python scripts/check_commit_attribution.py
python -m pytest -q
```

For release-level changes, also run the slow/reproduction checks documented in the repository.
