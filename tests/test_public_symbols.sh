#!/bin/sh
set -eu

lib="${1:-}"
if [ -z "$lib" ]; then
  echo "usage: $0 /path/to/libifxc.so" >&2
  exit 2
fi

if ! command -v nm >/dev/null 2>&1; then
  exit 0
fi

bad="$(
  nm -D --defined-only "$lib" \
    | awk '{print $3}' \
    | grep -Ev '^(ifxc_|IFXC_|__bss_start$|_edata$|_end$)' \
    || true
)"

if [ -n "$bad" ]; then
  echo "unexpected exported symbols:" >&2
  printf '%s\n' "$bad" >&2
  exit 1
fi
