#!/bin/bash

# Smart rebuild script for GWAS Alfred workflow
# - Checks if rebuild is already running
# - If running: reports status
# - Checks if current version is already built
# - If not: starts rebuild in background

# match the PATH the query script filters use, so the rebuild runs under the
# same interpreter the user installed pandas into
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WF_DATA="${alfred_workflow_data:-.}"
LOG_FILE="$WF_DATA/rebuild.log"
DB_FILE="$WF_DATA/index.db"
PID_FILE="$WF_DATA/rebuild.pid"

# Check if a build we started is still running. A pid file is used rather than
# pgrep, which also matches editors and log tails holding the script name.
RUNNING_PID=""
if [ -f "$PID_FILE" ]; then
    CANDIDATE=$(cat "$PID_FILE" 2>/dev/null)
    if [ -n "$CANDIDATE" ] && kill -0 "$CANDIDATE" 2>/dev/null &&
       ps -p "$CANDIDATE" -o command= 2>/dev/null | grep -q "build-GWAS-index.py"; then
        RUNNING_PID="$CANDIDATE"
    else
        rm -f "$PID_FILE"   # stale entry from a build that has since exited
    fi
fi

if [ -n "$RUNNING_PID" ]; then
    ELAPSED=$(ps -o etime= -p "$RUNNING_PID" 2>/dev/null | xargs)
    cat <<EOF
{
  "rerun": 3,
  "items": [{
    "title": "⏳️ Rebuild in progress...",
    "subtitle": "PID: $RUNNING_PID | Running for: $ELAPSED | You can close Alfred, it will continue",
    "valid": false,
    "icon": {"path": "icon.png"},
    "text": {
      "largetype": "GWAS Catalog Rebuild In Progress\n\nProcess ID: $RUNNING_PID\nRunning for: $ELAPSED\n\nYou can close Alfred - the rebuild will continue in the background"
    }
  }]
}
EOF
    exit 0
fi

# Pick an interpreter that can actually run the build
PYTHON=""
for CANDIDATE in $(command -v python3) /opt/homebrew/bin/python3 /usr/local/bin/python3 /usr/bin/python3; do
    if [ -x "$CANDIDATE" ] && "$CANDIDATE" -c "import pandas" 2>/dev/null; then
        PYTHON="$CANDIDATE"
        break
    fi
done

if [ -z "$PYTHON" ]; then
    cat <<EOF
{
  "items": [{
    "title": "⚠️ pandas is not installed",
    "subtitle": "The rebuild needs pandas: run 'pip3 install pandas' and try again",
    "valid": false,
    "icon": {"path": "icons/Warning.png"}
  }]
}
EOF
    exit 0
fi

# Check server for current release filename
SERVER_FILENAME=$(curl -sI "https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false" \
    | grep -i 'Content-Disposition' \
    | sed 's/.*filename=//;s/\r//')

if [ -n "$SERVER_FILENAME" ]; then
    TSV_NAME="${SERVER_FILENAME%.zip}.tsv"
    TSV_PATH="$WF_DATA/$TSV_NAME"

    # Extract colophon from filename (e.g. e115_r2026-02-16_full)
    COLOPHON=$(echo "$TSV_NAME" | sed 's/.*associations_//;s/\.tsv//')

    # Check if TSV exists and DB is complete (has GeneTrait table)
    if [ -f "$TSV_PATH" ]; then
        HAS_GENETRAIT=$(sqlite3 "$DB_FILE" "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='GeneTrait';" 2>/dev/null)
        if [ "$HAS_GENETRAIT" = "1" ]; then
            cat <<EOF
{
  "items": [{
    "title": "You already have the most current version of the database",
    "subtitle": "Version: $COLOPHON",
    "valid": false,
    "icon": {"path": "icon.png"}
  }]
}
EOF
            exit 0
        fi
    fi
else
    COLOPHON="latest"
fi

# Launch the build in background, fully detached from Alfred
cd "$SCRIPT_DIR"
nohup "$PYTHON" build-GWAS-index.py > "$LOG_FILE" 2>&1 &
PID=$!
echo "$PID" > "$PID_FILE"
disown $PID 2>/dev/null

cat <<EOF
{
  "items": [{
    "title": "🚀 GWAS catalog rebuild started",
    "subtitle": "Building version $COLOPHON in background (PID: $PID). Run ::rebuild to check status.",
    "valid": false,
    "icon": {"path": "icon.png"},
    "text": {
      "largetype": "GWAS Catalog Rebuild Started\n\nProcess ID: $PID\nVersion: $COLOPHON\n\n✓ You can close Alfred\n✓ Rebuild continues in background\n✓ Run ::rebuild again to check status"
    }
  }]
}
EOF
