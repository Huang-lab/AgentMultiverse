#!/usr/bin/env bash
# Keep agents in this project blind to completed analyses and to earlier session records.
# PreToolUse hook for Bash: refuse any command that names a protected location.
# The framework docs (README.md, framework/) describe the experiment itself, so agents may not read them either.
# Finished studies go into the archive only through tools/archive_study.sh STUDY, which prints sizes and counts, never contents.
# The shared result cache (.fmcache/) is reached only through tools/fmcache.{R,py}: its metadata would reveal other runs' choices.
# Relaxed while .claude/hooks/UNBLIND exists (delete it to restore blinding before launching new cohorts).
[ -e "$(dirname "$0")/UNBLIND" ] && exit 0
cmd=$(jq -r '.tool_input.command // empty')
case "$cmd" in
  *GWASAgentLite_archive*|*AgentMultiverse_archive*|*46fa2140-2c99-45cb-8d31-04cf397027a3*|*f614b217-f199-44f5-a9d8-797e197ce4fa*|*.claude/projects*.jsonl*|*.fmcache*|*[[:space:]]README.md*|README.md*|./README.md*|*AgentMultiverse/README.md*|*framework/AMA.md*|*framework/|*framework/\**)
    jq -n '{hookSpecificOutput: {hookEventName: "PreToolUse", permissionDecision: "deny",
      permissionDecisionReason: "Blocked: earlier analyses and session records are off-limits to agents in this project. Work from CLAUDE.md and the resources it lists."}}'
    ;;
esac
exit 0
