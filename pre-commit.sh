#!/usr/bin/env bash

# Load only staged changes
STASH_NAME=pre-commit-$(date +%s)
echo $STASH_NAME > $GIT_DIR/commit_info/stash
git stash save --include-untracked --keep-index --quiet $STASH_NAME

# Ensure documentation is up-to-date
R --slave -e "devtools::document()"

# Ensure documentation can be built
R --slave -e "devtools::build_manual(path = '/tmp/')"

git add man/ NAMESPACE
