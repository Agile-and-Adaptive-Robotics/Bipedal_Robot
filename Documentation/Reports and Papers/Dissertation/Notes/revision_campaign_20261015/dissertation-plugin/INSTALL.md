# INSTALL — `dissertation` plugin (v0.1.0)

**SUPERSEDED 2026-10-17:** the live plugin source of record is the ZCode_Skills repo,
`ZCode_Skills\plugins\dissertation\` (with `ZCode_Skills\plugins\marketplace.json`).
This scaffold copy is a campaign record; do not register it.

**Registration is Ben's call.** This scaffold is inert until he registers it; no session
should install or enable it on its own initiative. If Ben asks for it, the routes are below.

## What's here

```
dissertation-plugin/
├── .zcode-plugin/plugin.json          # manifest: name "dissertation", version 0.1.0
├── skills/diss-style-contract/SKILL.md  # style-contract enforcement skill
├── skills/diss-revise/SKILL.md          # revision workflow skill
├── agents/diss-supervisor/diss-supervisor.md  # four-gate supervisor persona
└── INSTALL.md                          # this file
```

The three canonical spec files (`style-contract.md`, `supervisor-rubric.md`,
`mr-nomenclature.md`) live one level up, **in the repo, outside the plugin**. The skills
reference them by repo-relative path, so the plugin only works from inside a checkout of
`Bipedal_Robot` — that is deliberate: the specs are versioned with the dissertation they
govern. If the plugin is copied into a plugin cache at install time, keep the repo present
and updated; edits to the specs belong in the repo copies, not in any cache.

## Route A — register as a local plugin

1. Open Zcode in `D:\Github\Bipedal_Robot`.
2. Run `/plugin` → Marketplaces → add a local marketplace pointing at a directory that
   contains `dissertation-plugin/` (e.g., this `revision_campaign_20261015/` folder).
3. Install the `dissertation` plugin from that marketplace and enable it.
4. Verify: the skills `diss-style-contract` and `diss-revise` and the agent `diss-supervisor`
   appear in the session's available lists.

Exact menu wording can differ between Zcode versions; the requirement is only that the
marketplace directory be the parent of `dissertation-plugin/` so the `.zcode-plugin/plugin.json`
manifest is found.

## Route B — junction the skills (Ben's existing pattern, lighter weight)

Ben already junctions repo-backed skills into `C:\Users\Ben Bolen\.zcode\skills\` (the
ZCode_Skills pattern). The same works here, skill by skill:

```
mklink /J "C:\Users\Ben Bolen\.zcode\skills\diss-style-contract" "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\revision_campaign_20261015\dissertation-plugin\skills\diss-style-contract"
mklink /J "C:\Users\Ben Bolen\.zcode\skills\diss-revise" "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\revision_campaign_20261015\dissertation-plugin\skills\diss-revise"
```

Run from an elevated cmd. This registers the two skills without the plugin manifest; the
agent definition would then be used by pointing a workflow at the `diss-supervisor.md`
file directly. Restart Zcode to pick them up.

## Notes

- This plugin contains no commands, hooks, or MCP servers — two skills and one agent
  persona only. Nothing runs automatically.
- On the other machines (laptop, easteregg2) the repo path differs; re-point the junctions
  at that machine's checkout, or use Route A from inside the checkout.
- If the specs move during the post-defense reorg, update the paths in both SKILL.md files
  and the agent file — they are the only place the spec location is written (plus this note).
