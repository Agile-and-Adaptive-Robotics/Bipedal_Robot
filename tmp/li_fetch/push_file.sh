#!/bin/bash
# Push files to easteregg2 over ssh stdin as base64 (scp rejects the
# spaced username; copy con hangs on EOF). Usage: push_file.sh <local> <remote>
set -e
KEY=~/.ssh/id_ed25519
HOST="ben bolen@easteregg2.mme.pdx.edu"
LOCAL="$1"
REMOTE="$2"
B64=$(base64 -w0 "$LOCAL")
B64FILE="${REMOTE}.b64"
# write the base64 payload remotely via PowerShell stdin (ASCII-safe)
printf '%s' "$B64" | ssh -i $KEY -o BatchMode=yes "$HOST" \
  "powershell -NoProfile -Command \"[IO.File]::WriteAllText('$B64FILE', [Console]::In.ReadToEnd())\""
# decode in place, clean up
ssh -i $KEY -o BatchMode=yes "$HOST" \
  "certutil -decode $B64FILE $REMOTE >nul 2>&1 && del $B64FILE && echo PUSHED $REMOTE"
