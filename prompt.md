# Prompt For Implementation Agent

You are working in `/home/simone/projects/reticolo`.

Before doing any code work, follow this exact git workflow:

1. Start from `rodney`:
   - `git checkout rodney`
2. Sync remotes:
   - `git fetch origin --prune`
3. Switch to a local branch that tracks `origin/translate_to_python`:
   - If local branch exists: `git checkout translate_to_python`
   - If not: `git checkout -b translate_to_python --track origin/translate_to_python`
4. Confirm you are on `translate_to_python`:
   - `git branch --show-current`

Then execute the full task:
- Translate the whole RETICOLO V9 implementation (`V9/reticolo_allege_v9`) into native Python.
- Make all necessary code changes and add/adjust tests.
- Run the required test suite and parity checks before finalizing.
- Commit your changes with clear commit messages.
- Push your branch to origin:
  - `git push origin translate_to_python`

After all work is complete and pushed, switch back to `rodney`:
- `git checkout rodney`

In your final report include:
- Files changed
- Commands run
- Test/parity results
- Confirmation that you returned to `rodney`
