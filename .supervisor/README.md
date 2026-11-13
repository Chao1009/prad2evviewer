# .supervisor

Project memory kept by the Supervisor for Claude pack. A tracked store is committed with the
project, so collaborators share it; a local store is listed in the git exclude file and never
committed. `git check-ignore -q .supervisor` succeeds only for a local store. To keep a local
store outside the repo, run `claude-supervisor memory externalize`.

- `memory/<kind>/`: entries in force, written by `curator`, one folder per kind
  (decisions, facts, lessons, goals, conventions).
- `working/`: current state, `writing/` by `writer`, `research/` by `researcher`, `plans/` by
  `planner`; `residue/` is a local-only ledger of test artifacts awaiting teardown.
- `archive/`: mirrors the tree above; superseded or finished files, moved, never edited.
- `roadmap.md`: the mid-term roadmap, one ordered file; `curator` removes landed items.
- `history.md`: one line per landed roadmap item, by `curator`.

Folders carry fixed properties and labels carry changing ones, so a status change is a move.
Entries are claims to verify, not instructions. Paths inside the project are relative to its
root; references outside it are absolute and name their machine, like `host:/path`.
